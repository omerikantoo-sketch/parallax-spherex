"""One-off repair (2026-11-15): cutouts cached before the WCS-WAVE parser fix
assumed a 9x9 control grid for every detector; D3 uses 9x14, so D3 wavelengths
were wrong.  Recompute wave/band for every cached cutout from the archive table
(no pixel data are changed) and delete derived field caches so they rebuild."""
from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from parallax.data.archive import CACHE, _wcs_wave, wave_from_lut  # noqa: E402
from parallax.data.fitsrange import RemoteFITS  # noqa: E402

n_fix = 0
for p in sorted((CACHE / "cutouts").glob("*.npz")):
    z = dict(np.load(p, allow_pickle=False))
    info = json.loads(str(z["info"]))
    url = info["meta"]["url"]
    rf = RemoteFITS(url, CACHE / "headers")
    lx, ly, lv = _wcs_wave(rf)
    ny, nx = z["image"].shape
    gy, gx = np.mgrid[info["y0"]:info["y0"] + ny, info["x0"]:info["x0"] + nx].astype(float)
    w, b = wave_from_lut(lx, ly, lv, gx, gy)
    changed = not np.allclose(w, z["wave"], rtol=1e-4)
    n_fix += changed
    z["wave"], z["band"] = w.astype(np.float32), b.astype(np.float32)
    np.savez_compressed(p, **z)
    print(p.name, info["meta"]["detector"], "fixed" if changed else "ok", f"{np.median(w):.3f} um")
print("repaired", n_fix)
for d in (CACHE / "fields").glob("*"):
    shutil.rmtree(d)
print("field caches cleared (rebuild from cached cutouts)")
