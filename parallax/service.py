"""Application services shared by the API and scripts: background jobs,
field assembly (real SPHEREx cutouts -> display stack + science frames),
spectra / light curves, astrometry and moving-object search on a field.
"""
from __future__ import annotations

import concurrent.futures as cf
import hashlib
import json
import logging
import threading
import time
import traceback
import uuid
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from astropy.time import Time

from .data import archive as A
from .data.model import Cutout, FLAG_BITS, REJECT_MASK
from .data.targets import TARGETS
from .inference.frames import Frame, frame_from_cutout, tangent_offsets, tangent_to_sky

log = logging.getLogger("parallax.service")
ROOT = Path(__file__).resolve().parents[1]
FIELDS = ROOT / "cache" / "fields"
VERSION = "0.3.0"


# =============================================================================
# jobs
# =============================================================================
@dataclass
class Job:
    id: str
    kind: str
    status: str = "queued"   # queued | running | done | error
    message: str = ""
    progress: float = 0.0
    result: dict | None = None
    error: str | None = None
    created: float = field(default_factory=time.time)


_jobs: dict[str, Job] = {}
_pool = cf.ThreadPoolExecutor(max_workers=3)
_jobs_lock = threading.Lock()


def submit(kind: str, fn, *args, **kw) -> Job:
    job = Job(id=uuid.uuid4().hex[:12], kind=kind)
    with _jobs_lock:
        _jobs[job.id] = job

    def progress(msg: str, frac: float):
        job.message, job.progress = msg, float(max(0.0, min(1.0, frac)))

    def run():
        job.status = "running"
        try:
            job.result = fn(*args, progress=progress, **kw)
            job.status, job.progress = "done", 1.0
        except Exception as exc:  # surfaced to the UI verbatim
            job.status, job.error = "error", f"{type(exc).__name__}: {exc}"
            log.error("job %s failed: %s", job.id, traceback.format_exc())

    _pool.submit(run)
    return job


def get_job(jid: str) -> Job | None:
    return _jobs.get(jid)


# =============================================================================
# fields
# =============================================================================
def field_id(ra: float, dec: float, size: int, release: str) -> str:
    return hashlib.sha1(f"{ra:.5f}|{dec:.5f}|{size}|{release}".encode()).hexdigest()[:12]


def _display_grid(size: int, pix: float = A.PIXSCALE_ARCSEC):
    c = (size - 1) / 2.0
    i, j = np.mgrid[0:size, 0:size].astype(float)
    xi = -(j - c) * pix       # east to the left
    eta = (c - i) * pix       # north up
    return xi, eta


def _resample(c: Cutout, ra0: float, dec0: float, size: int) -> tuple[np.ndarray, np.ndarray]:
    """Bilinear resampling of (IMAGE - ZODI) onto a north-up display grid.
    For display only -- science uses native pixels."""
    from scipy.ndimage import map_coordinates
    xi, eta = _display_grid(size)
    ra, dec = tangent_to_sky(xi.ravel(), eta.ravel(), ra0, dec0)
    px, py = c.wcs.all_world2pix(np.stack([ra, dec], -1), 0, quiet=True).T
    img = (c.image_mjysr - c.zodi_mjysr).astype(float)
    good = c.good().astype(float)
    imgf = np.where(good > 0, img, 0.0)
    v = map_coordinates(imgf, [py, px], order=1, mode="constant", cval=np.nan)
    g = map_coordinates(good, [py, px], order=1, mode="constant", cval=0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(g > 0.5, v / np.maximum(g, 1e-6), np.nan)
    flg = map_coordinates(c.flags.astype(float), [py, px], order=0, mode="constant", cval=float(1 << 9))
    return out.reshape(size, size).astype(np.float32), flg.reshape(size, size).astype(np.int64)


def _target_track(ra0, dec0, mjds, target_slug):
    if not target_slug or target_slug not in TARGETS:
        return None
    t = TARGETS[target_slug]
    ra, dec = t.at(np.asarray(mjds))
    xi, eta = tangent_offsets(ra, dec, ra0, dec0)
    return {"xi": np.atleast_1d(xi).tolist(), "eta": np.atleast_1d(eta).tolist(), "name": t.name,
            "source": t.source, "note": t.note}


def build_field(ra: float, dec: float, size: int = 64, release: str = "dr1", detectors=(1, 2, 3, 4, 5, 6),
                max_frames: int = 160, target: str | None = None, progress=None, allow_scan: bool = True) -> dict:
    fid = field_id(ra, dec, size, release)
    fdir = FIELDS / fid
    if (fdir / "field.json").exists():
        meta = json.loads((fdir / "field.json").read_text())
        if target and not meta.get("target"):
            meta["target"] = target
        return meta
    progress = progress or (lambda m, f: None)
    progress("locating SPHEREx exposures", 0.02)
    hits, how = A.find_coverage(ra, dec, release, progress=lambda m, f: progress(m, 0.02 + 0.4 * f),
                                allow_scan=allow_scan)
    hits = [h for h in hits if h["detector"] in detectors]
    if not hits:
        meta = {"id": fid, "ra": ra, "dec": dec, "size": size, "release": release, "discovery": how,
                "frames": [], "empty": True, "message": "No SPHEREx exposures found covering this position."}
        return meta
    if len(hits) > max_frames:  # keep an even sampling in time
        idx = np.unique(np.linspace(0, len(hits) - 1, max_frames).round().astype(int))
        hits = [hits[i] for i in idx]
    cut_size = size + 8
    cuts = A.fetch_cutouts(hits, ra, dec, size_px=cut_size, progress=lambda m, f: progress(m, 0.45 + 0.5 * f))
    if not cuts:
        raise RuntimeError("Exposures were found but no cutout could be retrieved from the archive.")
    progress("resampling for display", 0.96)
    fdir.mkdir(parents=True, exist_ok=True)
    imgs, flgs, recs = [], [], []
    for i, c in enumerate(cuts):
        im, fl = _resample(c, ra, dec, size)
        imgs.append(im)
        flgs.append(fl)
        ny, nx = c.shape
        s = c.summary()
        rec = {**s, "i": i, "utc": Time(c.meta.mjd_avg, format="mjd").isot[:19],
               "flag_frac": float(np.mean((c.flags & REJECT_MASK) != 0)),
               "sigma_ujy": float(np.nanmedian(np.sqrt(c.variance_mjysr2)) * float(np.mean(c.mjysr_to_ujy))),
               "zodi_mjysr": float(np.nanmedian(c.zodi_mjysr)), "notes": c.notes}
        recs.append(rec)
        A.save_cutout(c, fdir / f"cut_{i:03d}.npz")
    np.savez_compressed(fdir / "display.npz", images=np.stack(imgs), flags=np.stack(flgs))
    mjds = [r["mjd"] for r in recs]
    meta = {"id": fid, "ra": ra, "dec": dec, "size": size, "release": release, "discovery": how,
            "pixscale_arcsec": A.PIXSCALE_ARCSEC, "frames": recs, "target": target,
            "track": _target_track(ra, dec, mjds, target), "created_utc": Time.now().isot[:19],
            "units": {"image": "MJy/sr (IMAGE - ZODI), bilinear-resampled north-up for display"},
            "provenance": "real", "version": VERSION}
    (fdir / "field.json").write_text(json.dumps(meta, default=float))
    return meta


def list_fields() -> list[dict]:
    out = []
    for p in sorted(FIELDS.glob("*/field.json")):
        try:
            m = json.loads(p.read_text())
            out.append({k: m.get(k) for k in ("id", "ra", "dec", "size", "release", "target", "created_utc")}
                       | {"n_frames": len(m.get("frames", []))})
        except Exception:
            pass
    return out


def load_field(fid: str) -> dict:
    return json.loads((FIELDS / fid / "field.json").read_text())


def field_display(fid: str) -> tuple[np.ndarray, np.ndarray]:
    z = np.load(FIELDS / fid / "display.npz")
    return z["images"], z["flags"]


def field_frames(fid: str, use_epsf: bool = True) -> tuple[dict, list[Frame], list[Cutout]]:
    meta = load_field(fid)
    cuts = [A.load_cutout(FIELDS / fid / f"cut_{r['i']:03d}.npz") for r in meta["frames"]]
    frames = [frame_from_cutout(c, meta["ra"], meta["dec"], use_epsf=use_epsf) for c in cuts]
    return meta, frames, cuts


# =============================================================================
# analyses on a field
# =============================================================================
def spectrum_at(fid: str, xi: float, eta: float, follow_target: bool = False) -> dict:
    from .inference.scene import forced_photometry
    meta, frames, cuts = field_frames(fid)
    if follow_target and meta.get("track"):
        tx = np.array(meta["track"]["xi"]) + xi
        ty = np.array(meta["track"]["eta"]) + eta
    else:
        tx, ty = xi, eta
    ph = forced_photometry(frames, tx, ty)
    ra, dec = tangent_to_sky(np.atleast_1d(np.mean(tx)), np.atleast_1d(np.mean(ty)), meta["ra"], meta["dec"])
    pts = [{"t": float(t), "utc": Time(t, format="mjd").isot[:19], "wave": float(w), "band": float(b),
            "flux": float(f), "err": float(e), "det": int(d), "label": l, "chi2": float(c)}
           for t, w, b, f, e, d, l, c in zip(ph["t"], ph["wave"], ph["band"], ph["flux"], ph["err"], ph["det"],
                                            ph["label"], ph["chi2"])]
    return {"field": fid, "xi": float(np.mean(tx)), "eta": float(np.mean(ty)), "ra": float(ra[0]),
            "dec": float(dec[0]), "points": pts, "units": {"flux": "uJy", "wave": "um"},
            "method": "PSF-weighted forced photometry (DR1 ePSF), ZODI and local median background subtracted",
            "provenance": "real"}


def astrometry(fid: str, target: str | None = None) -> dict:
    """Centroid a bright target in every frame and fit proper motion."""
    from .inference.scene import fit_astrometry, psf_centroid
    meta, frames, cuts = field_frames(fid)
    slug = target or meta.get("target")
    tr = meta.get("track")
    pts = []
    for k, f in enumerate(frames):
        if tr:
            x0, y0 = f.to_pix(tr["xi"][k], tr["eta"][k])
        else:
            x0, y0 = f.to_pix(0.0, 0.0)
        r = psf_centroid(f, float(x0), float(y0))
        if r is None or not (r["flux"] > 0) or r["sx"] > 0.5 or r["sy"] > 0.5:
            continue
        # pixel -> tangent offset via inverse affine
        Ainv = np.linalg.inv(f.A)
        xi, eta = Ainv @ (np.array([r["x"], r["y"]]) - f.b)
        J = Ainv @ np.diag([r["sx"], r["sy"]])
        sxi, seta = float(np.hypot(*J[0])), float(np.hypot(*J[1]))
        pts.append({"t": f.t_mjd, "xi": float(xi), "eta": float(eta), "sxi": max(sxi, 0.05), "seta": max(seta, 0.05),
                    "det": f.detector, "wave": f.wave_um, "flux": float(r["flux"]), "chi2r": r["chi2r"], "k": k})
    if len(pts) < 4:
        return {"field": fid, "n": len(pts), "error": "too few usable centroids"}
    t = np.array([p["t"] for p in pts])
    xi = np.array([p["xi"] for p in pts]); eta = np.array([p["eta"] for p in pts])
    sxi = np.array([p["sxi"] for p in pts]); seta = np.array([p["seta"] for p in pts])
    sc = np.array([c.meta.sc_pos_km or (0, 0, 0) for c in cuts], float)
    used = [p["k"] for p in pts]   # frame index of every centroid (simultaneous detectors share t)
    t_ref = float(np.round(np.mean(t), 1))
    plx = TARGETS[slug].plx if slug in TARGETS else None
    fit = fit_astrometry(t, xi, eta, sxi, seta, meta["ra"], meta["dec"], t_ref, fit_parallax=False,
                         plx_fixed_mas=plx, sc_geo_km=sc[used])
    expect = None
    if slug in TARGETS:
        tg = TARGETS[slug]
        expect = {"pmra_masyr": tg.pmra, "pmdec_masyr": tg.pmdec, "plx_mas": tg.plx, "source": tg.source}
        # predicted offsets at t_ref relative to the field centre
        ra_p, dec_p = tg.at(t_ref)
        exi, eeta = tangent_offsets(ra_p, dec_p, meta["ra"], meta["dec"])
        expect["xi0_arcsec"], expect["eta0_arcsec"] = float(exi), float(eeta)
    return {"field": fid, "target": slug, "points": pts, "fit": fit, "expected": expect,
            "method": "per-frame ePSF least-squares centroid; WLS fit of offset + proper motion with parallax "
                      "fixed to the catalogue value (observer = SPHEREx: Earth barycentre + header spacecraft position)",
            "provenance": "real"}


def search_field(fid: str, vmax_arcsec_day: float = 2.0, subtract_static: bool = True, progress=None) -> dict:
    """Moving-object search on a real field (static scene subtracted)."""
    from .inference.baselines import setup, statistics_for
    from .inference.matched import Trajectories, build_maps, gather, make_grid
    from .inference.scene import static_catalog, subtract_static as sub
    progress = progress or (lambda m, f: None)
    meta, frames, cuts = field_frames(fid)
    # restrict to one survey visit (the densest 25-day window) for linear motion
    t = np.array([f.t_mjd for f in frames])
    best_lo, best_n = t.min(), 0
    for lo in t:
        n = int(np.sum((t >= lo) & (t < lo + 25)))
        if n > best_n:
            best_lo, best_n = lo, n
    sel = [i for i in range(len(frames)) if best_lo <= t[i] < best_lo + 25]
    frames = [frames[i] for i in sel]
    progress("building matched-filter maps", 0.05)
    maps = build_maps(frames, q=2)
    cat = (np.array([]), np.array([]), np.array([]))
    if subtract_static:
        progress("modelling static sources", 0.15)
        cat = static_catalog(maps, snr=6.0)
        frames, amps = sub(frames, cat[0], cat[1])
        maps = build_maps(frames, q=2)
    if subtract_static and cat[0].size:
        # (catalogue z saturates for bright stars, so the mask radius is set from
        # the fitted PSF flux instead: 2 px at <=1 mJy, +2 px per decade)
        flux_med = np.nanmedian(np.abs(amps), axis=0)
        radius = 2.0 + 2.0 * np.log10(np.maximum(flux_med / 1000.0, 1.0))
        # mask residual cores of static sources (subtraction is imperfect for
        # bright / saturated stars): radius grows with source significance
        for f in frames:
            yy, xx = np.mgrid[0:f.shape[0], 0:f.shape[1]]
            for a_, b_, r in zip(cat[0], cat[1], radius):
                x, y = f.to_pix(a_, b_)
                f.ivar = np.where(np.hypot(xx - x, yy - y) <= r, 0.0, f.ivar)
        maps = build_maps(frames, q=2)
    size = meta["size"]
    grid = make_grid(frames, half_fov_arcsec=0.4 * size * A.PIXSCALE_ARCSEC, vmax_arcsec_day=vmax_arcsec_day)
    det, S, B2, _, best = setup(maps, None)
    px, py = grid.positions()
    vx, vy = grid.velocities()
    span = float(np.ptp([f.t_mjd for f in frames]))
    vmin = A.PIXSCALE_ARCSEC / max(span, 1e-3)        # must move >= 1 px across the visit

    def scan(mp, label, f0, f1):
        bz = np.full(px.size, -np.inf)
        bv = np.zeros(px.size, int)
        for iv in range(vx.size):
            if np.hypot(vx[iv], vy[iv]) < vmin:
                continue
            if iv % 20 == 0:
                progress(f"{label} {iv}/{vx.size}", f0 + (f1 - f0) * iv / vx.size)
            traj = Trajectories(px, py, np.full(px.size, vx[iv]), np.full(px.size, vy[iv]), grid.t_ref)
            psi, phi = gather(mp, traj)
            z = statistics_for(psi, phi, det, S, B2, None, best)["bank"]
            u = z > bz
            bz[u], bv[u] = z[u], iv
        return bz, bv

    best_z, best_v = scan(maps, "trajectory search", 0.2, 0.6)
    # empirical null on this field: same frames with their epochs randomly
    # permuted (destroys coherent motion, keeps residual/noise structure)
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(frames))
    scr = []
    for f, j in zip(frames, perm):
        nf = Frame(**{**f.__dict__})
        nf.t_mjd = frames[j].t_mjd
        scr.append(nf)
    null_z, _ = scan(build_maps(scr, q=2), "time-scrambled null", 0.6, 0.98)
    null_max = float(np.max(null_z[np.isfinite(null_z)])) if np.isfinite(null_z).any() else float("nan")
    order = np.argsort(best_z)[::-1]
    cands = []
    for i in order:
        if len(cands) >= 8:
            break
        c = (px[i], py[i], vx[best_v[i]], vy[best_v[i]], best_z[i])
        if any(np.hypot(c[0] - o["xi0"], c[1] - o["eta0"]) < 12 for o in cands):
            continue
        cands.append({"xi0": float(c[0]), "eta0": float(c[1]), "vxi": float(c[2]), "veta": float(c[3]),
                      "z_bank": float(c[4]), "above_scrambled_null": bool(c[4] > null_max)})
    # evidence for the top candidate: per-frame contributions
    ev = None
    if cands:
        c = cands[0]
        traj = Trajectories(np.array([c["xi0"]]), np.array([c["eta0"]]), np.array([c["vxi"]]), np.array([c["veta"]]),
                            grid.t_ref)
        psi, phi = gather(maps, traj)
        ev = {"t": [f.t_mjd for f in frames], "det": det.tolist(), "wave": [f.wave_um for f in frames],
              "psi": psi[:, 0].tolist(), "phi": phi[:, 0].tolist()}
    return {"field": fid, "n_frames": len(frames), "visit_start_mjd": float(best_lo), "n_trials": int(grid.n_trials),
            "vmax_arcsec_day": vmax_arcsec_day, "vmin_arcsec_day": vmin, "static_sources": int(cat[0].size),
            "scrambled_null_max": null_max, "candidates": cands,
            "evidence": ev, "t_ref": grid.t_ref,
            "significance_note": "z_bank is a single-trial Gaussian-equivalent statistic maximised over the "
                                 "template bank; it is NOT corrected for the %d trajectories searched. Compare with "
                                 "the empirically calibrated thresholds in reports/motion_benchmark.json. Field-level check: the maximum "
                                 "of the same search on time-scrambled frames is %.1f; only candidates above it "
                                 "are flagged (one scramble, so this is a rough guide, not a calibrated FAP)."
                                 % (grid.n_trials, null_max),
            "status": "Exploratory (real data, uncalibrated for this field)"}


# =============================================================================
# synthetic motion demonstration (provenance: synthetic)
# =============================================================================
def motion_demo(mag: float = 21.0, seed: int = 11, sed: str = "reflected", vxi: float = 1.6, veta: float = -0.9) -> dict:
    """A faint synthetic mover under the real SPHEREx burst cadence, with the
    per-frame evidence for its true trajectory under two flux models."""
    from .inference.baselines import setup
    from .inference.matched import Trajectories, build_maps, gather, stat_bank
    from .science.simulate import MovingSource, SimConfig, simulate
    from .science.spectra import Spectrum, band_average

    spec = {"reflected": Spectrum("reflected sunlight (slope +0.1/um)", "reflect", (0.1, 250.0, 0.0)),
            "thermal": Spectrum("reflected + thermal (T=300 K)", "reflect", (0.1, 300.0, 2.0)),
            "browndwarf": Spectrum("cold brown dwarf proxy (450 K)", "browndwarf", (450.0,))}[sed]
    cfg = SimConfig(size=24, seed=seed, spherex_bursts=True, n_bursts=20, span_days=20.0, star_density_per_deg2=0.0)
    f0, tr0 = simulate(cfg, [], include_stars=False)
    mv = MovingSource(10.0, -6.0, vxi, veta, mag, spec, t_ref=tr0["t_ref"])
    frames, _ = simulate(cfg, [mv], include_stars=False)
    order = np.argsort([f.t_mjd for f in frames])
    frames = [frames[i] for i in order]
    maps = build_maps(frames, q=2)
    det, S, B2, _, best = setup(maps, spec)
    traj = Trajectories(np.array([mv.xi0]), np.array([mv.eta0]), np.array([vxi]), np.array([veta]), mv.t_ref)
    psi, phi = gather(maps, traj)
    psi, phi = psi[:, 0], phi[:, 0]
    # cumulative statistics in time order
    cum_stack = np.cumsum(psi) / np.sqrt(np.maximum(np.cumsum(phi), 1e-30))
    num = np.cumsum(S * psi[None], axis=1)
    den = np.cumsum(S * S * phi[None], axis=1)
    cum_bank = np.max(num / np.sqrt(np.maximum(den, 1e-30)), axis=0)
    single = psi / np.sqrt(np.maximum(phi, 1e-30))
    k_best = int(np.argmax(num[:, -1] / np.sqrt(den[:, -1])))
    # stamps and matched stacks (moving vs stationary hypothesis)
    half = 5
    stamps, mov, sta = [], np.zeros((2 * half + 1,) * 2), np.zeros((2 * half + 1,) * 2)
    wm, ws = np.zeros_like(mov), np.zeros_like(sta)
    s_true = np.array([float(band_average(spec, f.wave_um, f.band_um)) for f in frames])
    for e, f in enumerate(frames):
        x, y = f.to_pix(*mv.pos(f.t_mjd))
        xs, ys = f.to_pix(*mv.pos(mv.t_ref))
        def cut(cx, cy):
            ix, iy = int(round(float(cx))), int(round(float(cy)))
            d = np.zeros((2 * half + 1,) * 2); w = np.zeros_like(d)
            ny, nx = f.shape
            for dy in range(-half, half + 1):
                for dx in range(-half, half + 1):
                    yy, xx = iy + dy, ix + dx
                    if 0 <= yy < ny and 0 <= xx < nx:
                        d[dy + half, dx + half] = f.data[yy, xx]; w[dy + half, dx + half] = f.ivar[yy, xx]
            return d, w
        d, w = cut(x, y)
        stamps.append({"t": f.t_mjd, "det": f.detector, "wave": f.wave_um, "data": (d / max(np.median(1 / np.sqrt(w[w > 0])) if np.any(w > 0) else 1, 1e-9)).round(3).tolist(),
                       "snr_frame": float(single[e]), "flux_true": float(mv.flux(f.wave_um, f.band_um))})
        mov += s_true[e] * w * d; wm += s_true[e] ** 2 * w
        d2, w2 = cut(xs, ys)
        sta += s_true[e] * w2 * d2; ws += s_true[e] ** 2 * w2
    mov = np.where(wm > 0, mov / np.sqrt(wm), 0)
    sta = np.where(ws > 0, sta / np.sqrt(ws), 0)
    return {"provenance": "synthetic", "mag_ab_1um": mag, "sed": spec.name, "velocity_arcsec_day": [vxi, veta],
            "n_frames": len(frames), "t": [f.t_mjd for f in frames], "det": det.tolist(),
            "wave": [f.wave_um for f in frames], "cum_stack": cum_stack.tolist(), "cum_bank": cum_bank.tolist(),
            "single": single.tolist(), "best_template": k_best, "stamps": stamps,
            "stack_moving": mov.round(3).tolist(), "stack_static": sta.round(3).tolist(),
            "final": {"stack": float(cum_stack[-1]), "bank": float(cum_bank[-1]), "max_single": float(single.max())},
            "note": "Stamps are per-frame data in units of the per-pixel noise. The matched stacks weight each "
                    "frame by s_e * ivar along the true (moving) or stationary hypothesis."}
