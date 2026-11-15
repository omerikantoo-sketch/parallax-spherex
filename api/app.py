"""PARALLAX HTTP API (FastAPI) + static Observatory frontend.

Run:  python -m uvicorn api.app:app --port 8600
"""
from __future__ import annotations

import io
import json
import logging
import re
import sys
import time
from pathlib import Path

import numpy as np
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from parallax import service as S  # noqa: E402
from parallax.data import archive as A  # noqa: E402
from parallax.data.targets import TARGETS  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s %(message)s")
logging.getLogger("httpx").setLevel(logging.WARNING)
log = logging.getLogger("parallax.api")

app = FastAPI(title="PARALLAX API", version=S.VERSION,
              description="SPHEREx time-domain viewer and physical-inference engine. "
                          "All synthetic content is labelled with provenance='synthetic'.")


def _clean(o):
    """Replace NaN/inf (not valid JSON) by None, recursively."""
    if isinstance(o, float):
        return o if np.isfinite(o) else None
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_clean(v) for v in o]
    return o


def _json(obj, status=200):
    return JSONResponse(_clean(json.loads(json.dumps(obj, default=_default))), status_code=status)


def _default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if hasattr(o, "__dict__"):
        return {k: v for k, v in o.__dict__.items() if not k.startswith("_")}
    return str(o)


# ----------------------------------------------------------------------------
_status_cache: dict = {"t": 0.0, "v": None}


@app.get("/api/status")
def status():
    """Archive reachability (cached 120 s)."""
    if time.time() - _status_cache["t"] > 120:
        def probe(url, params=None):
            try:
                r = A.client().get(url, params=params, timeout=8.0)
                return r.status_code
            except Exception:
                return None
        sia = probe(A.SIA_URL, {"COLLECTION": "spherex_dr1", "POS": "CIRCLE 0 0 0.0001", "MAXREC": "1"})
        s3 = probe(A.S3_BASE, {"list-type": "2", "max-keys": "1", "prefix": "dr1/"})
        _status_cache["v"] = {"irsa_sia": "online" if sia == 200 else f"unavailable ({sia})",
                              "aws_mirror": "online" if s3 == 200 else f"unavailable ({s3})",
                              "checked_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                              "version": S.VERSION}
        _status_cache["t"] = time.time()
    return _status_cache["v"]


@app.get("/api/targets")
def targets():
    out = []
    for t in TARGETS.values():
        ra, dec = t.at(61000.0)
        out.append({"slug": t.slug, "name": t.name, "kind": t.kind, "ra": float(ra), "dec": float(dec),
                    "epoch_note": "position propagated to MJD 61000 (2025-11-23)", "source": t.source,
                    "note": t.note, "pmra": t.pmra, "pmdec": t.pmdec, "plx": t.plx})
    return out


_SEXA = re.compile(r"^\s*(\d{1,2})[h:\s]+(\d{1,2})[m:\s]+([\d.]+)s?\s*[, ]\s*([+-]?\d{1,2})[d:\s°]+(\d{1,2})[m:\s']+([\d.]+)")


@app.get("/api/resolve")
def resolve(q: str):
    """Parse 'ra dec' (decimal degrees or sexagesimal) or a challenge-target name."""
    s = q.strip()
    for t in TARGETS.values():
        if s.lower() in (t.slug, t.name.lower()):
            ra, dec = t.at(61000.0)
            return {"ra": float(ra), "dec": float(dec), "target": t.slug, "name": t.name}
    m = _SEXA.match(s)
    if m:
        from astropy.coordinates import Angle
        ra = Angle(f"{m.group(1)}h{m.group(2)}m{m.group(3)}s").deg
        dec = Angle(f"{m.group(4)}d{m.group(5)}m{m.group(6)}s").deg
        return {"ra": float(ra), "dec": float(dec), "target": None}
    nums = re.findall(r"[+-]?\d+(?:\.\d+)?", s)
    if len(nums) >= 2:
        ra, dec = float(nums[0]), float(nums[1])
        if 0 <= ra < 360 and -90 <= dec <= 90:
            return {"ra": ra, "dec": dec, "target": None}
    raise HTTPException(400, "Could not parse coordinates. Use 'RA Dec' in degrees or 08h55m10.8s -07d14m42s.")


# ----------------------------------------------------------------------------
@app.post("/api/fields")
def create_field(ra: float, dec: float, size: int = Query(64, ge=16, le=160), release: str = "dr1",
                 target: str | None = None, max_frames: int = Query(160, ge=4, le=400), scan: bool = True):
    if release not in ("dr1", "qr3"):
        raise HTTPException(400, "release must be dr1 or qr3")
    if not (0 <= ra < 360 and -90 <= dec <= 90):
        raise HTTPException(400, "coordinates out of range")
    fid = S.field_id(ra, dec, size, release)
    if (S.FIELDS / fid / "field.json").exists():
        return {"job": None, "field": fid, "cached": True}
    job = S.submit("field", S.build_field, ra, dec, size, release, max_frames=max_frames, target=target,
                   allow_scan=scan)
    return {"job": job.id, "field": fid, "cached": False}


@app.get("/api/fields")
def fields():
    return S.list_fields()


@app.get("/api/fields/{fid}")
def field_meta(fid: str):
    try:
        return _json(S.load_field(fid))
    except FileNotFoundError:
        raise HTTPException(404, "field not built")


@app.get("/api/fields/{fid}/stack.f32")
def field_stack(fid: str):
    """All display frames as little-endian float32 [n, size, size] (NaN = no valid data)."""
    try:
        imgs, _ = S.field_display(fid)
    except FileNotFoundError:
        raise HTTPException(404, "field not built")
    return Response(imgs.astype("<f4").tobytes(), media_type="application/octet-stream",
                    headers={"X-Shape": ",".join(map(str, imgs.shape)), "Cache-Control": "max-age=3600"})


@app.get("/api/fields/{fid}/flags.i32")
def field_flags(fid: str):
    _, flg = S.field_display(fid)
    return Response(flg.astype("<i4").tobytes(), media_type="application/octet-stream",
                    headers={"X-Shape": ",".join(map(str, flg.shape))})


@app.get("/api/fields/{fid}/spectrum")
def field_spectrum(fid: str, xi: float = 0.0, eta: float = 0.0, follow: bool = False):
    return _json(S.spectrum_at(fid, xi, eta, follow_target=follow))


@app.get("/api/fields/{fid}/astrometry")
def field_astrometry(fid: str):
    return _json(S.astrometry(fid))


@app.post("/api/fields/{fid}/search")
def field_search(fid: str, vmax: float = Query(2.0, gt=0, le=20)):
    job = S.submit("search", S.search_field, fid, vmax)
    return {"job": job.id}


@app.get("/api/fields/{fid}/gravity")
def field_gravity(fid: str, xi: float = 0.0, eta: float = 0.0, follow: bool = True):
    """Gravity-mode hypothesis comparison on real forced photometry at a position."""
    from parallax.gravity.pspl import GravityModel, LightCurve
    sp = S.spectrum_at(fid, xi, eta, follow_target=follow)
    pts = sp["points"]
    if len(pts) < 12:
        raise HTTPException(422, "fewer than 12 usable photometric points")
    # real photometry: add a 1% relative calibration floor in quadrature (detector-to-detector and
    # ePSF-zone systematics are not in the archive variance)
    err = [float(np.hypot(p["err"], 0.01 * abs(p["flux"]))) for p in pts]
    lc = LightCurve([p["t"] for p in pts], [p["flux"] for p in pts], err,
                    [p["wave"] for p in pts], label=f"field {fid}", provenance="real",
                    meta={"error_floor": "1% relative, added in quadrature"})
    r = GravityModel(lc, k_spec=3).compare()
    return _json(_gravity_payload(lc, r, real=True))


def _gravity_payload(lc, r, real=False):
    def fr(f):
        return {"model": f.model, "chi2": f.chi2, "n_par": f.n_par, "params": f.params,
                "model_flux": f.model_flux.tolist()}
    return {"provenance": lc.provenance, "label": lc.label, "t": lc.t.tolist(), "flux": lc.flux.tolist(),
            "err": lc.err.tolist(), "wave": lc.wave.tolist(), "meta": {k: v for k, v in lc.meta.items() if k != "true_flux"},
            "true_flux": lc.meta.get("true_flux", np.array([])).tolist() if not real else None,
            "H0": fr(r["H0"]), "H1": fr(r["H1"]), "H2": fr(r["H2"]),
            "dchi2_H1_vs_H0": r["dchi2_H1_vs_H0"], "dchi2_H1_vs_H2": r["dchi2_H1_vs_H2"],
            "dBIC_H1_vs_H0": r["dBIC_H1_vs_H0"], "dBIC_H1_vs_H2": r["dBIC_H1_vs_H2"], "verdict": r["verdict"],
            "n": r["n"], "baseline_flux": _h1_baseline(lc, r), "chi2r": r.get("chi2r"),
            "lens_physical": r.get("lens_physical")}


def _h1_baseline(lc, r):
    """Unlensed (A = 1) flux of the best lens model: source + blend SEDs."""
    from parallax.science.spectra import spectral_basis
    k = 3
    c = r["H1"].coef
    return (spectral_basis(lc.wave, k) @ c[:k] + spectral_basis(lc.wave, 2 * k) @ c[k:]).tolist()


@app.get("/api/gravity/demo")
def gravity_demo(u0: float = Query(0.3, gt=0, le=3), tE: float = Query(15.0, gt=0.2, le=300),
                 blend: float = Query(0.5, ge=0, lt=1), cadence: str = "deep", lensed: bool = True,
                 seed: int = 7, sigma: float = Query(0.02, gt=0, le=0.5)):
    """Synthetic (clearly labelled) microlensing event under a SPHEREx-like cadence."""
    from parallax.gravity.pspl import GravityModel, mass_posterior, simulate_event
    rng = np.random.default_rng(seed)
    if cadence == "deep":   # deep field: a sample every ~0.6 d over 200 d (illustrative)
        t = np.sort(rng.uniform(0, 200, 330))
    else:                   # wide survey: two ~20-day visits 180 d apart, bursts of 4
        bursts = np.concatenate([rng.uniform(0, 20, 25), rng.uniform(180, 200, 25)])
        t = np.sort(np.concatenate([b + np.arange(4) * 0.0015 for b in bursts]))
    wave = np.exp(rng.uniform(np.log(0.75), np.log(5.0), t.size))
    t0 = float(rng.uniform(40, 160)) if cadence == "deep" else 10.0
    lc = simulate_event(t, wave, t0=t0, u0=u0, tE=tE, blend_frac_2um=blend, sigma_frac=sigma, seed=seed,
                        lensed=lensed)
    r = GravityModel(lc, k_spec=3).compare()
    out = _gravity_payload(lc, r)
    out["cadence"] = cadence
    out["mass"] = mass_posterior(r["H1"].params["tE"], max(0.05 * r["H1"].params["tE"], 0.1))
    return _json(out)


@app.get("/api/motion/demo")
def motion_demo(mag: float = Query(21.0, ge=15, le=25), seed: int = 11, sed: str = "reflected"):
    if sed not in ("reflected", "thermal", "browndwarf"):
        raise HTTPException(400, "sed must be reflected, thermal or browndwarf")
    return _json(S.motion_demo(mag, seed, sed))


@app.get("/api/reports")
def reports():
    out = {}
    for p in sorted((ROOT / "reports").glob("*.json")):
        if p.stem.startswith("_"):      # internal artefacts (e.g. reused null calibration)
            continue
        try:
            d = json.loads(p.read_text())
            d.pop("records", None)
            out[p.stem] = d
        except Exception as exc:
            out[p.stem] = {"error": str(exc)}
    return _json(out)


@app.get("/api/jobs/{jid}")
def job(jid: str):
    j = S.get_job(jid)
    if j is None:
        raise HTTPException(404, "unknown job")
    return _json({"id": j.id, "kind": j.kind, "status": j.status, "message": j.message, "progress": j.progress,
                  "error": j.error, "result": j.result if j.status == "done" else None})


# ----------------------------------------------------------------------------
WEB = ROOT / "web"


@app.get("/")
def index():
    return FileResponse(WEB / "index.html")


app.mount("/", StaticFiles(directory=WEB, html=True), name="web")
