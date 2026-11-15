"""Random-access reader for remote, uncompressed multi-extension FITS files.

SPHEREx Level-2 spectral-image MEFs are ~65 MB each but are stored as
uncompressed IEEE images.  Every HDU therefore has a deterministic byte layout
(FITS Standard 4.0, Sect. 3-7): a header made of 2880-byte blocks followed by a
data unit padded to a multiple of 2880 bytes.  Knowing the header lets us fetch
only the rows of an image that intersect a sky cutout with an HTTP Range
request, instead of downloading the entire file.

This module is deliberately small and dependency-light: astropy parses headers,
numpy interprets the big-endian payload.
"""
from __future__ import annotations

import hashlib
import logging
import os
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path

import httpx
import numpy as np
from astropy.io import fits

log = logging.getLogger("parallax.fitsrange")

BLOCK = 2880
CARD = 80

_client_lock = threading.Lock()
_client: httpx.Client | None = None


def client() -> httpx.Client:
    """Shared keep-alive HTTP client (thread-safe for concurrent requests)."""
    global _client
    with _client_lock:
        if _client is None:
            _client = httpx.Client(
                timeout=httpx.Timeout(60.0, connect=20.0),
                limits=httpx.Limits(max_connections=64, max_keepalive_connections=32),
                headers={"User-Agent": "parallax-research/0.1 (SPHEREx range reader)"},
                follow_redirects=True,
            )
        return _client


class RemoteUnavailable(RuntimeError):
    """Raised when the remote archive cannot serve a request after retries."""


def fetch_range(url: str, start: int, stop: int, retries: int = 4) -> bytes:
    """Fetch bytes [start, stop) from ``url`` with retry and exponential backoff."""
    if stop <= start:
        return b""
    hdr = {"Range": f"bytes={start}-{stop - 1}"}
    delay = 0.5
    last: Exception | None = None
    for attempt in range(retries):
        try:
            r = client().get(url, headers=hdr)
            if r.status_code == 206:
                return r.content
            if r.status_code == 200:  # server ignored Range: slice locally
                return r.content[start:stop]
            if r.status_code in (403, 404):
                raise RemoteUnavailable(f"{r.status_code} for {url}")
            last = RemoteUnavailable(f"HTTP {r.status_code} for {url}")
        except RemoteUnavailable:
            raise
        except Exception as exc:  # network error: retry
            last = exc
        time.sleep(delay)
        delay *= 2
    raise RemoteUnavailable(f"failed after {retries} attempts: {url}: {last}")


@dataclass
class HDUInfo:
    index: int
    name: str
    header: fits.Header
    header_offset: int
    data_offset: int
    data_nbytes: int  # unpadded

    @property
    def shape(self) -> tuple[int, ...]:
        n = self.header.get("NAXIS", 0)
        # FITS axis order is fastest-first; numpy wants slowest-first.
        return tuple(int(self.header[f"NAXIS{i}"]) for i in range(n, 0, -1))

    @property
    def dtype(self) -> np.dtype:
        return np.dtype({8: ">u1", 16: ">i2", 32: ">i4", 64: ">i8", -32: ">f4", -64: ">f8"}[
            int(self.header["BITPIX"])])


def _data_nbytes(h: fits.Header) -> int:
    naxis = int(h.get("NAXIS", 0))
    if naxis == 0:
        return 0
    n = 1
    for i in range(1, naxis + 1):
        n *= int(h[f"NAXIS{i}"])
    gcount = int(h.get("GCOUNT", 1))
    pcount = int(h.get("PCOUNT", 0))
    return abs(int(h["BITPIX"])) // 8 * gcount * (pcount + n)


def _padded(n: int) -> int:
    return ((n + BLOCK - 1) // BLOCK) * BLOCK


@dataclass
class RemoteFITS:
    """Lazily-indexed remote MEF.  Header blocks are cached on disk."""

    url: str
    cache_dir: Path | None = None
    hdus: list[HDUInfo] = field(default_factory=list)
    _complete: bool = False

    # ------------------------------------------------------------------ headers
    def _read_header_at(self, offset: int, guess_blocks: int = 16) -> tuple[fits.Header, int]:
        """Return (header, n_header_bytes) for the HDU starting at ``offset``."""
        buf = b""
        nblocks = guess_blocks
        while True:
            chunk = fetch_range(self.url, offset + len(buf), offset + nblocks * BLOCK)
            if not chunk:
                raise RemoteUnavailable(f"truncated FITS at {offset} in {self.url}")
            buf += chunk
            # look for END card on 80-char boundaries
            for i in range(0, len(buf), CARD):
                if buf[i:i + 8] == b"END     " or buf[i:i + 80].rstrip() == b"END":
                    hbytes = _padded(i + CARD)
                    return fits.Header.fromstring(buf[:i + CARD].decode("ascii", "replace")), hbytes
            nblocks *= 2

    def _cache_path(self) -> Path | None:
        if self.cache_dir is None:
            return None
        key = hashlib.sha1(self.url.encode()).hexdigest()[:20]
        return self.cache_dir / f"hdr_{key}.txt"

    def index(self, upto: int | None = None) -> list[HDUInfo]:
        """Index HDUs (all, or the first ``upto``+1).  Uses on-disk cache."""
        if self._complete or (upto is not None and len(self.hdus) > upto):
            return self.hdus
        cp = self._cache_path()
        if cp is not None and cp.exists() and not self.hdus:
            self._load_cache(cp)
            if self._complete or (upto is not None and len(self.hdus) > upto):
                return self.hdus
        offset = 0 if not self.hdus else self.hdus[-1].data_offset + _padded(self.hdus[-1].data_nbytes)
        i = len(self.hdus)
        while upto is None or i <= upto:
            try:
                h, hb = self._read_header_at(offset, guess_blocks=1 if i == 0 else 16)
            except RemoteUnavailable:
                if i > 0:  # reached EOF
                    self._complete = True
                    break
                raise
            nb = _data_nbytes(h)
            name = str(h.get("EXTNAME", "PRIMARY" if i == 0 else f"HDU{i}")).strip()
            self.hdus.append(HDUInfo(i, name, h, offset, offset + hb, nb))
            offset = offset + hb + _padded(nb)
            i += 1
            if name == "WCS-WAVE":  # last HDU in SPHEREx L2 MEFs
                self._complete = True
                break
        if cp is not None:
            self._save_cache(cp)
        return self.hdus

    def _save_cache(self, cp: Path) -> None:
        cp.parent.mkdir(parents=True, exist_ok=True)
        lines = [f"#complete {int(self._complete)}"]
        for h in self.hdus:
            lines.append(f"#hdu {h.index} {h.name} {h.header_offset} {h.data_offset} {h.data_nbytes}")
            lines.append(h.header.tostring(sep="\n"))
        tmp = cp.with_suffix(".tmp")
        tmp.write_text("\n".join(lines), encoding="ascii", errors="replace")
        os.replace(tmp, cp)

    def _load_cache(self, cp: Path) -> None:
        text = cp.read_text(encoding="ascii", errors="replace").split("\n")
        self._complete = text[0].strip().endswith("1")
        cur: list[str] = []
        meta = None
        def flush():
            if meta is not None:
                idx, name, ho, do, nb = meta
                self.hdus.append(HDUInfo(int(idx), name, fits.Header.fromstring("\n".join(cur), sep="\n"),
                                         int(ho), int(do), int(nb)))
        for line in text[1:]:
            if line.startswith("#hdu "):
                flush()
                meta = line.split()[1:]
                cur = []
            else:
                cur.append(line)
        flush()

    def hdu(self, name: str) -> HDUInfo:
        for h in self.index():
            if h.name == name:
                return h
        raise KeyError(name)

    # --------------------------------------------------------------------- data
    def read_rows(self, name: str, y0: int, y1: int) -> np.ndarray:
        """Read full rows [y0, y1) of a 2-D image HDU (0-based numpy indexing)."""
        h = self.hdu(name)
        ny, nx = h.shape
        y0, y1 = max(0, y0), min(ny, y1)
        if y1 <= y0:
            return np.zeros((0, nx), dtype=h.dtype.newbyteorder("="))
        rowb = nx * h.dtype.itemsize
        raw = fetch_range(self.url, h.data_offset + y0 * rowb, h.data_offset + y1 * rowb)
        return np.frombuffer(raw, dtype=h.dtype).reshape(y1 - y0, nx).astype(h.dtype.newbyteorder("="))

    def read_plane(self, name: str, k: int) -> np.ndarray:
        """Read plane ``k`` of a 3-D cube HDU (e.g. one ePSF zone)."""
        h = self.hdu(name)
        nz, ny, nx = h.shape
        pb = ny * nx * h.dtype.itemsize
        raw = fetch_range(self.url, h.data_offset + k * pb, h.data_offset + (k + 1) * pb)
        return np.frombuffer(raw, dtype=h.dtype).reshape(ny, nx).astype(h.dtype.newbyteorder("="))

    def read_all(self, name: str) -> bytes:
        h = self.hdu(name)
        return fetch_range(self.url, h.data_offset, h.data_offset + h.data_nbytes)
