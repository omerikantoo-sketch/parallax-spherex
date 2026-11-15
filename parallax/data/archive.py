"""SPHEREx archive access: discovery of exposures covering a sky position and
science-ready cutouts with IMAGE/VARIANCE/FLAGS/ZODI/ePSF/wavelength.

Discovery strategy (in order):

1. IRSA SIA2 (``https://irsa.ipac.caltech.edu/SIA``) -- the official service.
2. A local footprint index built from the public AWS Open Data mirror
   (``s3://nasa-irsa-spherex``).  The mirror has no spatial index, so we use a
   physics prior: SPHEREx scans great circles near ecliptic longitude
   lambda_sun +/- ~90-120 deg, so only planning periods whose solar longitude
   places the target near those scan longitudes can contain it.  Candidate
   periods are first *sampled* sparsely (every Nth header), then scanned densely
   where a sampled pointing lands near the target, and finally every detector
   header is tested for exact containment with its full TAN-SIP WCS.

   This is itself a coarse-to-fine search with an explicit (conservative)
   rejection radius; it never fabricates coverage -- a position is reported as
   covered only after an exact WCS test on the real header.

Cutouts are fetched with HTTP Range requests (see ``fitsrange``).
"""
from __future__ import annotations

import concurrent.futures as cf
import datetime as dt
import hashlib
import io
import json
import logging
import re
import sqlite3
import threading
import time
import warnings
import xml.etree.ElementTree as ET
from dataclasses import asdict
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import GeocentricTrueEcliptic, SkyCoord, get_sun
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS, FITSFixedWarning

from .fitsrange import RemoteFITS, RemoteUnavailable, client, fetch_range
from .model import Cutout, ObservationMeta

warnings.simplefilter("ignore", FITSFixedWarning)
log = logging.getLogger("parallax.archive")

S3_BASE = "https://nasa-irsa-spherex.s3.amazonaws.com/"
IRSA_BASE = "https://irsa.ipac.caltech.edu/ibe/data/spherex/"
SIA_URL = "https://irsa.ipac.caltech.edu/SIA"
NPIX = 2040
PIXSCALE_ARCSEC = 6.15
ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "cache"


# =============================================================================
# S3 listing (with on-disk cache)
# =============================================================================
class Catalog:
    def __init__(self, path: Path | None = None):
        self.path = path or (CACHE / "index.sqlite")
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.Lock()
        self.db = sqlite3.connect(self.path, check_same_thread=False)
        self.db.executescript(
            """
            CREATE TABLE IF NOT EXISTS listing(prefix TEXT PRIMARY KEY, fetched REAL, body TEXT);
            CREATE TABLE IF NOT EXISTS pointing(
                key TEXT PRIMARY KEY, release TEXT, period TEXT, obsid TEXT, detector INT,
                mjd_avg REAL, crval1 REAL, crval2 REAL, elon REAL, elat REAL, header TEXT);
            CREATE INDEX IF NOT EXISTS pointing_period ON pointing(period, detector);
            CREATE TABLE IF NOT EXISTS coverage(qkey TEXT PRIMARY KEY, fetched REAL, body TEXT);
            CREATE TABLE IF NOT EXISTS sampled(period_dir TEXT PRIMARY KEY, step INT, dense INT);
            """
        )
        self.db.commit()

    def get(self, sql: str, args=()):
        with self._lock:
            return self.db.execute(sql, args).fetchall()

    def put(self, sql: str, args=()):
        with self._lock:
            self.db.execute(sql, args)
            self.db.commit()

    def putmany(self, sql: str, rows):
        with self._lock:
            self.db.executemany(sql, rows)
            self.db.commit()


_catalog: Catalog | None = None


def catalog() -> Catalog:
    global _catalog
    if _catalog is None:
        _catalog = Catalog()
    return _catalog


def list_prefix(prefix: str, delimiter: bool = True, max_age_s: float = 7 * 86400) -> tuple[list[str], list[tuple[str, int]]]:
    """List an S3 prefix -> (sub-prefixes, [(key, size)]).  Cached."""
    cat = catalog()
    row = cat.get("SELECT fetched, body FROM listing WHERE prefix=?", (prefix + ("/" if delimiter else "#"),))
    if row and time.time() - row[0][0] < max_age_s:
        body = json.loads(row[0][1])
        return body["p"], [tuple(k) for k in body["k"]]
    prefixes: list[str] = []
    keys: list[tuple[str, int]] = []
    token = None
    ns = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}
    while True:
        params = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if delimiter:
            params["delimiter"] = "/"
        if token:
            params["continuation-token"] = token
        for attempt in range(4):
            try:
                r = client().get(S3_BASE, params=params)
                r.raise_for_status()
                break
            except Exception:
                if attempt == 3:
                    raise RemoteUnavailable(f"S3 listing failed for {prefix}")
                time.sleep(1 + attempt)
        root = ET.fromstring(r.content)
        prefixes += [e.text for e in root.findall("s:CommonPrefixes/s:Prefix", ns)]
        keys += [(c.find("s:Key", ns).text, int(c.find("s:Size", ns).text)) for c in root.findall("s:Contents", ns)]
        nt = root.find("s:NextContinuationToken", ns)
        if nt is None:
            break
        token = nt.text
    cat.put("INSERT OR REPLACE INTO listing VALUES (?,?,?)",
            (prefix + ("/" if delimiter else "#"), time.time(), json.dumps({"p": prefixes, "k": keys})))
    return prefixes, keys


# =============================================================================
# Survey geometry prior
# =============================================================================
_PERIOD_RE = re.compile(r"(\d{4})W(\d{2})_(\d)([A-Z])")


def period_midtime(period: str) -> Time:
    m = _PERIOD_RE.match(period)
    if not m:
        raise ValueError(period)
    y, w, half = int(m.group(1)), int(m.group(2)), int(m.group(3))
    d = dt.datetime.combine(dt.date.fromisocalendar(y, w, 1), dt.time(0)) + dt.timedelta(days=0.875 + 3.5 * (half - 1))
    return Time(d, scale="utc")


def sun_lon(t: Time) -> float:
    s = get_sun(t).transform_to(GeocentricTrueEcliptic(equinox=t))
    return float(s.lon.deg)


def ecliptic(ra: float, dec: float) -> tuple[float, float]:
    c = SkyCoord(ra * u.deg, dec * u.deg).transform_to(GeocentricTrueEcliptic())
    return float(c.lon.deg), float(c.lat.deg)


def _wrap180(x):
    return (np.asarray(x) + 180.0) % 360.0 - 180.0


def period_could_cover(period: str, ra: float, dec: float, slack_deg: float = 22.0) -> bool:
    """Conservative prior on whether a planning period can observe (ra, dec).

    Empirically (verified on DR1 headers) SPHEREx boresight ecliptic longitudes
    in a period lie within ~lambda_sun + (+85..+120) or lambda_sun + (-100..-85)
    deg (great circles through the ecliptic poles).  At |beta| > 60 deg all
    longitudes converge, so the prior is disabled there.
    """
    lam, beta = ecliptic(ra, dec)
    if abs(beta) > 70:
        return True
    elong = solar_elongation(period, lam, beta)
    # Boresight elongations measured from DR1 headers span ~91-116 deg; the
    # 3-detector field extends ~6 deg from boresight and a period spans 3.5 d
    # (~3.5 deg of solar motion).  ``slack_deg`` widens this conservatively.
    return bool(91.0 - 6.0 - 2.0 - slack_deg / 4 < elong < 116.0 + 6.0 + 2.0 + slack_deg / 4)


def solar_elongation(period: str, lam: float, beta: float) -> float:
    ls = sun_lon(period_midtime(period))
    c = np.cos(np.radians(beta)) * np.cos(np.radians(lam - ls))
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


# =============================================================================
# Header scanning
# =============================================================================
def _read_image_header(key: str) -> fits.Header | None:
    """Read the IMAGE (HDU 1) header with a single range request."""
    url = S3_BASE + key
    try:
        b = fetch_range(url, 2880, 2880 + 2880 * 13)
    except RemoteUnavailable:
        return None
    s = b.decode("ascii", "replace")
    for i in range(0, len(s), 80):
        if s[i:i + 80].rstrip() == "END":
            return fits.Header.fromstring(s[:i + 80])
    return None


_KEEP = ["NAXIS1", "NAXIS2", "CRPIX1", "CRPIX2", "PC1_1", "PC1_2", "PC2_1", "PC2_2", "CDELT1", "CDELT2",
         "CUNIT1", "CUNIT2", "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "LONPOLE", "LATPOLE", "RADESYS",
         "A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER", "DETECTOR", "OBSID", "MJD-BEG", "MJD-AVG", "MJD-END",
         "XPOSURE", "DATE-AVG", "PSF_FWHM", "OMEGA_MEDIAN", "X_SC", "Y_SC", "Z_SC", "SPS_ELON", "SPS_ELAT",
         "L2DQAFLG", "BUNIT"]


def _slim_header(h: fits.Header) -> str:
    out = fits.Header()
    for k, v in h.items():
        if k in _KEEP or re.match(r"^(A|B|AP|BP)_\d_\d$", k):
            out[k] = v
    # astropy strips HIERARCH prefix on read; preserve a few useful values
    for hk, nk in [("OMEGA_MEDIAN", "OMEGAMED")]:
        if hk in h:
            out[nk] = h[hk]
    return out.tostring()


def _store_headers(rows: list[tuple[str, str, str, fits.Header]]) -> None:
    recs = []
    for key, release, period, h in rows:
        recs.append((key, release, period, str(h.get("OBSID", "")).strip(), int(h.get("DETECTOR", 0)),
                     float(h.get("MJD-AVG", np.nan)), float(h["CRVAL1"]), float(h["CRVAL2"]),
                     float(h.get("SPS_ELON", np.nan)), float(h.get("SPS_ELAT", np.nan)), _slim_header(h)))
    catalog().putmany("INSERT OR REPLACE INTO pointing VALUES (?,?,?,?,?,?,?,?,?,?,?)", recs)


def scan_keys(keys: list[str], release: str, period: str, workers: int = 32) -> int:
    have = {r[0] for r in catalog().get(
        f"SELECT key FROM pointing WHERE key IN ({','.join('?' * len(keys))})", keys)} if keys else set()
    todo = [k for k in keys if k not in have]
    out = []
    with cf.ThreadPoolExecutor(workers) as ex:
        for k, h in zip(todo, ex.map(_read_image_header, todo)):
            if h is not None and "CRVAL1" in h:
                out.append((k, release, period, h))
    if out:
        _store_headers(out)
    return len(out)


def period_dirs(release: str) -> list[str]:
    periods, _ = list_prefix(f"{release}/level2/")
    return [p for p in periods if _PERIOD_RE.search(p)]


def detector_keys(period_dir: str, detector: int) -> list[str]:
    """All MEF keys for one detector in a planning-period directory (all versions)."""
    versions, _ = list_prefix(period_dir)
    keys = []
    for v in versions:
        _, ks = list_prefix(f"{v}{detector}/", delimiter=False)
        keys += [k for k, _ in ks if k.endswith(".fits")]
    # prefer the newest processing of an OBSID (retry dirs re-process failures)
    by_obs: dict[str, str] = {}
    for k in sorted(keys):
        m = re.search(r"level2_(.+?)D\d_spx", k)
        if m:
            by_obs[m.group(1)] = k
    return [by_obs[o] for o in sorted(by_obs)]


def angsep(ra1, dec1, ra2, dec2):
    ra1, dec1, ra2, dec2 = map(np.radians, (ra1, dec1, ra2, dec2))
    c = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


def header_contains(hdr_text: str, ra: float, dec: float, margin_px: float = 0.0) -> tuple[bool, float, float]:
    h = fits.Header.fromstring(hdr_text)
    w = WCS(h)
    # cheap TAN-only pre-test (SIP terms are < ~2 px here); the iterative SIP
    # inverse can diverge for points far outside the detector.
    x, y = w.wcs_world2pix([[ra, dec]], 0)[0]
    if not (np.isfinite(x) and np.isfinite(y)) or not (-60 - margin_px <= x <= NPIX + 60 + margin_px
                                                      and -60 - margin_px <= y <= NPIX + 60 + margin_px):
        return False, float(x), float(y)
    try:
        x, y = w.all_world2pix([[ra, dec]], 0, quiet=True)[0]
    except Exception:
        return False, float(x), float(y)
    ok = (-margin_px <= x <= NPIX - 1 + margin_px) and (-margin_px <= y <= NPIX - 1 + margin_px)
    return bool(ok), float(x), float(y)


def find_coverage_s3(ra: float, dec: float, release: str = "dr1", progress=None,
                     sample_step: int = 12, near_deg: float = 14.0, max_periods: int | None = None) -> list[dict]:
    """Find all detector exposures in ``release`` whose footprint contains (ra, dec)."""
    qkey = f"{release}:{ra:.4f}:{dec:.4f}"
    row = catalog().get("SELECT body FROM coverage WHERE qkey=?", (qkey,))
    if row:
        return json.loads(row[0][0])
    pdirs = [p for p in period_dirs(release) if period_could_cover(_PERIOD_RE.search(p).group(0), ra, dec)]
    if max_periods:
        pdirs = pdirs[:max_periods]
    log.info("candidate periods for (%.4f, %.4f): %d", ra, dec, len(pdirs))
    hits: list[dict] = []
    for ip, pdir in enumerate(pdirs):
        period = _PERIOD_RE.search(pdir).group(0)
        if progress:
            progress(f"sampling {period}", ip / max(1, len(pdirs)))
        k1 = detector_keys(pdir, 1)
        if not k1:
            continue
        sample = k1[::sample_step]
        scan_keys(sample, release, period)
        rows = catalog().get(
            f"SELECT key, crval1, crval2 FROM pointing WHERE key IN ({','.join('?' * len(sample))})", sample)
        if not rows:
            continue
        d = angsep(np.array([r[1] for r in rows]), np.array([r[2] for r in rows]), ra, dec)
        if d.min() > near_deg + 12.0:  # consecutive pointings move <~ 25 deg
            continue
        # dense scan of detector-1 headers in this period
        if progress:
            progress(f"scanning {period}", ip / max(1, len(pdirs)))
        scan_keys(k1, release, period)
        rows = catalog().get(
            f"SELECT key, obsid, crval1, crval2 FROM pointing WHERE period=? AND detector=1", (period,))
        near = [r for r in rows if angsep(r[2], r[3], ra, dec) < near_deg]
        # all six detectors of the nearby exposures
        cand_keys = []
        for r in near:
            for det in range(1, 7):
                cand_keys.append(r[0].replace(f"/1/level2_", f"/{det}/level2_").replace("D1_spx", f"D{det}_spx"))
        scan_keys(cand_keys, release, period)
        if cand_keys:
            for key, hdr, det, mjd, obsid in catalog().get(
                    f"SELECT key, header, detector, mjd_avg, obsid FROM pointing WHERE key IN "
                    f"({','.join('?' * len(cand_keys))})", cand_keys):
                ok, x, y = header_contains(hdr, ra, dec, margin_px=0)
                if ok:
                    hits.append({"key": key, "release": release, "period": period, "obsid": obsid,
                                 "detector": det, "mjd": mjd, "x": x, "y": y})
    hits.sort(key=lambda h: (h["mjd"], h["detector"]))
    catalog().put("INSERT OR REPLACE INTO coverage VALUES (?,?,?)", (qkey, time.time(), json.dumps(hits)))
    return hits


# =============================================================================
# SIA (official) discovery -- used when IRSA is reachable
# =============================================================================
def find_coverage_sia(ra: float, dec: float, collection: str = "spherex_dr1", timeout: float = 45.0) -> list[dict] | None:
    try:
        r = client().get(SIA_URL, params={"COLLECTION": collection, "POS": f"CIRCLE {ra} {dec} 0.0005"},
                         timeout=timeout)
    except Exception as exc:
        log.warning("SIA unreachable: %s", exc)
        return None
    if r.status_code != 200 or b"VOTABLE" not in r.content[:2000].upper():
        log.warning("SIA returned HTTP %s", r.status_code)
        return None
    from astropy.io.votable import parse_single_table
    t = parse_single_table(io.BytesIO(r.content)).to_table()
    hits = []
    for row in t:
        url = str(row["access_url"])
        if "/spherex/" not in url or not url.endswith(".fits"):
            continue
        key = url.split("/ibe/data/spherex/")[-1]
        m = re.search(r"D(\d)_spx", key)
        hits.append({"key": key, "release": key.split("/")[0], "period": _PERIOD_RE.search(key).group(0),
                     "obsid": re.search(r"level2_(.+?)D\d_spx", key).group(1),
                     "detector": int(m.group(1)) if m else 0,
                     "mjd": float(row["t_min"]) if "t_min" in t.colnames else np.nan})
    hits.sort(key=lambda h: (h["mjd"], h["detector"]))
    return hits


def find_coverage(ra: float, dec: float, release: str = "dr1", progress=None, allow_scan: bool = True) -> tuple[list[dict], str]:
    hits = find_coverage_sia(ra, dec, collection=f"spherex_{release}")
    if hits is not None:
        return hits, "IRSA SIA2"
    if not allow_scan:
        return [], "unavailable"
    return find_coverage_s3(ra, dec, release, progress=progress), "AWS mirror footprint index (IRSA SIA unavailable)"


# =============================================================================
# Cutouts
# =============================================================================
def _rice_rows(rf: RemoteFITS, y0: int, y1: int) -> np.ndarray:
    """Decode rows [y0, y1) of the tile-compressed (1 row per tile) FLAGS HDU."""
    from astropy.io.fits.hdu.compressed._codecs import Rice1

    h = rf.hdu("FLAGS")
    hdr = h.header
    nx = int(hdr["ZNAXIS1"])
    assert int(hdr["ZTILE2"]) == 1 and hdr["ZCMPTYPE"].strip() == "RICE_1"
    rowlen = int(hdr["NAXIS1"])  # 8 bytes: (count, offset) int32 descriptor
    raw = fetch_range(rf.url, h.data_offset + y0 * rowlen, h.data_offset + y1 * rowlen)
    desc = np.frombuffer(raw, dtype=">i4").reshape(-1, 2)
    heap0 = h.data_offset + int(hdr.get("THEAP", rowlen * int(hdr["NAXIS2"])))
    lo, hi = int(desc[:, 1].min()), int((desc[:, 1] + desc[:, 0]).max())
    heap = fetch_range(rf.url, heap0 + lo, heap0 + hi)
    codec = Rice1(blocksize=int(hdr.get("ZVAL1", 32)), bytepix=int(hdr.get("ZVAL2", 4)), tilesize=nx)
    out = np.empty((y1 - y0, nx), dtype=np.int32)
    for i, (cnt, off) in enumerate(desc):
        out[i] = codec.decode(np.frombuffer(heap, dtype=np.uint8, count=cnt, offset=off - lo))
    return out


def _wcs_wave(rf: RemoteFITS) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Spectral WCS control grid.  Grid sizes differ by detector (e.g. 9x9 for
    D1, 9x14 for D3), so dimensions are read from TFORMn/TDIMn, never assumed."""
    import re
    hdr = rf.hdu("WCS-WAVE").header
    raw = rf.read_all("WCS-WAVE")
    nx = int(re.match(r"(\d+)", hdr["TFORM1"]).group(1))
    ny = int(re.match(r"(\d+)", hdr["TFORM2"]).group(1))
    nv = int(re.match(r"(\d+)", hdr["TFORM3"]).group(1))
    assert nv == 2 * nx * ny, (nx, ny, nv)
    o = 0
    x = np.frombuffer(raw[o:o + 4 * nx], dtype=">i4").astype(float); o += 4 * nx
    y = np.frombuffer(raw[o:o + 4 * ny], dtype=">i4").astype(float); o += 4 * ny
    # FITS TDIM3 = (2, nx, ny): fastest axis is (wave, band), then x, then y
    vals = np.frombuffer(raw[o:o + 4 * nv], dtype=">f4").reshape(ny, nx, 2).astype(float)
    return x, y, vals  # vals[iy, ix, (wave, band)] -- orientation checked in tests/test_core.py


def wave_from_lut(x: np.ndarray, y: np.ndarray, vals: np.ndarray, px: np.ndarray, py: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Bilinear interpolation of the WCS-WAVE table at 0-based pixel coords."""
    from scipy.interpolate import RegularGridInterpolator

    fx, fy = px + 1.0, py + 1.0  # table is in 1-based FITS pixel coordinates
    pts = np.stack([np.clip(fy, y[0], y[-1]).ravel(), np.clip(fx, x[0], x[-1]).ravel()], axis=-1)
    w = RegularGridInterpolator((y, x), vals[..., 0])(pts).reshape(px.shape)
    b = RegularGridInterpolator((y, x), vals[..., 1])(pts).reshape(px.shape)
    return w, b


def _epsf_row(rf: RemoteFITS, x: float, y: float) -> dict:
    h = rf.hdu("EPSF")
    nrow = int(h.header["NAXIS2"])
    rowlen = int(h.header["NAXIS1"])
    # tiling: 21x21 (det 1,2,4,5,6) or 11x41 (det 3); ordering x-fast then y.
    if nrow == 441:
        nbx = nby = 21
    else:
        nbx, nby = 11, nrow // 11
    best = None
    bx = min(nbx - 1, max(0, int(x / (NPIX / nbx))))
    by = min(nby - 1, max(0, int(y / (NPIX / nby))))
    for cand in {by * nbx + bx}:
        raw = fetch_range(rf.url, h.data_offset + cand * rowlen, h.data_offset + (cand + 1) * rowlen)
        rec = {
            "binx": int(np.frombuffer(raw[0:4], ">i4")[0]), "biny": int(np.frombuffer(raw[4:8], ">i4")[0]),
            "xc": float(np.frombuffer(raw[8:12], ">f4")[0]), "yc": float(np.frombuffer(raw[12:16], ">f4")[0]),
            "nstar": int(np.frombuffer(raw[24:28], ">i4")[0]),
            "cwave": float(np.frombuffer(raw[32:36], ">f4")[0]), "cband": float(np.frombuffer(raw[36:40], ">f4")[0]),
            "epsf": np.frombuffer(raw[40:40 + 1089 * 8], ">f8").reshape(33, 33).astype(float),
        }
        best = rec
    best["ovs"] = int(h.header.get("OVSMPX", 5))
    return best


def _meta_from_header(key: str, h: fits.Header, release: str) -> ObservationMeta:
    sc = (h.get("X_SC"), h.get("Y_SC"), h.get("Z_SC"))
    return ObservationMeta(
        obs_id=str(h.get("OBSID", "")).strip(), detector=int(h["DETECTOR"]), release=release,
        url=S3_BASE + key, mjd_beg=float(h["MJD-BEG"]), mjd_avg=float(h["MJD-AVG"]), mjd_end=float(h["MJD-END"]),
        exptime_s=float(h["XPOSURE"]), sc_pos_km=tuple(float(s) for s in sc) if None not in sc else None,
        psf_fwhm_arcsec=float(h["PSF_FWHM"]) if "PSF_FWHM" in h else None,
        omega_median_arcsec2=float(h["OMEGA_MEDIAN"]) if "OMEGA_MEDIAN" in h else None,
    )


def fetch_cutout(key: str, ra: float, dec: float, size_px: int = 64, release: str | None = None,
                 use_cache: bool = True) -> Cutout | None:
    """Fetch a square native-pixel cutout centred on (ra, dec).  None if off-detector."""
    release = release or key.split("/")[0]
    ck = hashlib.sha1(f"{key}|{ra:.6f}|{dec:.6f}|{size_px}".encode()).hexdigest()[:24]
    cpath = CACHE / "cutouts" / f"{ck}.npz"
    if use_cache and cpath.exists():
        return load_cutout(cpath)
    rf = RemoteFITS(S3_BASE + key, cache_dir=CACHE / "headers")
    hdus = rf.index()
    h = rf.hdu("IMAGE").header
    w = WCS(h)
    xc, yc = w.all_world2pix([[ra, dec]], 0)[0]
    if not (0 <= xc < NPIX and 0 <= yc < NPIX):
        return None
    half = size_px // 2
    x0, y0 = int(round(xc)) - half, int(round(yc)) - half
    x1, y1 = x0 + size_px, y0 + size_px
    ry0, ry1 = max(0, y0), min(NPIX, y1)
    rx0, rx1 = max(0, x0), min(NPIX, x1)

    def pad(rows, fill):
        out = np.full((size_px, size_px), fill, dtype=rows.dtype)
        out[ry0 - y0: ry1 - y0, rx0 - x0: rx1 - x0] = rows[:, rx0:rx1]
        return out

    with cf.ThreadPoolExecutor(4) as ex:
        f_img = ex.submit(rf.read_rows, "IMAGE", ry0, ry1)
        f_var = ex.submit(rf.read_rows, "VARIANCE", ry0, ry1)
        f_zod = ex.submit(rf.read_rows, "ZODI", ry0, ry1)
        f_flg = ex.submit(_rice_rows, rf, ry0, ry1)
        img = pad(f_img.result().astype(np.float32), np.nan)
        var = pad(f_var.result().astype(np.float32), np.nan)
        zod = pad(f_zod.result().astype(np.float32), np.nan)
        flg = pad(f_flg.result(), np.int32(1 << 9))  # off-detector -> MISSING_DATA
    lx, ly, lv = _wcs_wave(rf)
    gy, gx = np.mgrid[y0:y1, x0:x1].astype(float)
    wave, band = wave_from_lut(lx, ly, lv, gx, gy)
    ep = _epsf_row(rf, xc, yc)
    hc = h.copy()
    hc["CRPIX1"] = h["CRPIX1"] - x0
    hc["CRPIX2"] = h["CRPIX2"] - y0
    # SIP polynomials are defined relative to CRPIX, so shifting CRPIX and the
    # pixel grid together leaves the distortion model unchanged.
    wc = WCS(hc)
    meta = _meta_from_header(key, h, release)
    omega = meta.omega_median_arcsec2 or PIXSCALE_ARCSEC ** 2
    cut = Cutout(meta=meta, image_mjysr=img, variance_mjysr2=var, flags=flg, zodi_mjysr=zod,
                 wave_um=wave.astype(np.float32), band_um=band.astype(np.float32), wcs=wc, x0=x0, y0=y0,
                 psf_oversampled=ep["epsf"], psf_oversamp=ep["ovs"], psf_zone_center=(ep["xc"], ep["yc"]),
                 omega_arcsec2=omega,
                 notes=[f"ePSF zone ({ep['binx']},{ep['biny']}) nstar={ep['nstar']} cwave={ep['cwave']:.4f}um",
                        "pixel solid angle: OMEGA_MEDIAN header value (per-pixel SAPM not applied)"])
    if use_cache:
        save_cutout(cut, cpath, hc)
    return cut


def save_cutout(c: Cutout, path: Path, header: fits.Header | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    hdr = header if header is not None else c.wcs.to_header(relax=True)
    meta = asdict(c.meta)
    np.savez_compressed(
        path, image=c.image_mjysr, variance=c.variance_mjysr2, flags=c.flags, zodi=c.zodi_mjysr,
        wave=c.wave_um, band=c.band_um, psf=c.psf_oversampled,
        info=json.dumps({"meta": meta, "x0": c.x0, "y0": c.y0, "ovs": c.psf_oversamp,
                         "psf_zone_center": c.psf_zone_center, "omega": float(np.mean(c.omega_arcsec2)),
                         "wave_source": c.wave_source, "provenance": c.provenance, "notes": c.notes,
                         "header": hdr.tostring()}))


def load_cutout(path: Path) -> Cutout:
    z = np.load(path, allow_pickle=False)
    info = json.loads(str(z["info"]))
    m = info["meta"]
    if m.get("sc_pos_km") is not None:
        m["sc_pos_km"] = tuple(m["sc_pos_km"])
    return Cutout(meta=ObservationMeta(**m), image_mjysr=z["image"], variance_mjysr2=z["variance"],
                  flags=z["flags"], zodi_mjysr=z["zodi"], wave_um=z["wave"], band_um=z["band"],
                  wcs=WCS(fits.Header.fromstring(info["header"])), x0=info["x0"], y0=info["y0"],
                  psf_oversampled=z["psf"], psf_oversamp=info["ovs"],
                  psf_zone_center=tuple(info["psf_zone_center"]), omega_arcsec2=info["omega"],
                  wave_source=info["wave_source"], provenance=info["provenance"], notes=info["notes"])


def fetch_cutouts(hits: list[dict], ra: float, dec: float, size_px: int = 64, workers: int = 8,
                  progress=None) -> list[Cutout]:
    out: list[Cutout] = []
    done = 0
    with cf.ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(fetch_cutout, h["key"], ra, dec, size_px): h for h in hits}
        for f in cf.as_completed(futs):
            done += 1
            try:
                c = f.result()
                if c is not None:
                    out.append(c)
            except Exception as exc:
                log.warning("cutout failed for %s: %s", futs[f]["key"], exc)
            if progress:
                progress(f"cutouts {done}/{len(hits)}", done / max(1, len(hits)))
    out.sort(key=lambda c: (c.meta.mjd_avg, c.meta.detector))
    return out
