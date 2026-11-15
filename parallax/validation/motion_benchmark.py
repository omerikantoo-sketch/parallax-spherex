"""Injection-recovery benchmark for faint slow-moving sources (synthetic,
SPHEREx-like visit structure).

Protocol (fixed seeds; development / calibration / test splits by seed range):

1. Null calibration (seeds 1000+): noise-only fields with the real SPHEREx
   burst cadence; run the *full* trajectory search and record the global
   maximum of every statistic.  The per-method detection threshold is the
   (1 - FAP) quantile of that maximum, so every method operates at the same
   per-field false-alarm probability (look-elsewhere effect included
   empirically).
2. Injection (seeds 5000+, never used for calibration): one moving source per
   field at a known AB magnitude, velocity and SED class.  A method detects
   the source if its statistic on grid trajectories adjacent to the truth
   reaches its threshold.
3. Completeness vs magnitude per method, Wilson 68% intervals, and the 50%
   completeness magnitude by linear interpolation.

Usage:  python -m parallax.validation.motion_benchmark --null 40 --per-mag 20
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ..inference.baselines import link_count_at_truth, link_false_tracklets, local_stats, multi_search
from ..inference.matched import build_maps, make_grid
from ..science.simulate import MovingSource, SimConfig, simulate
from ..science.spectra import Spectrum

ROOT = Path(__file__).resolve().parents[2]

SED_CLASSES = {
    "reflected (asteroid/TNO-like)": lambda r: Spectrum("refl", "reflect", (float(r.uniform(-0.05, 0.3)), 250.0, 0.0)),
    "reflected + thermal (warm NEO-like)": lambda r: Spectrum("refl+th", "reflect", (0.1, float(r.uniform(250, 330)), float(r.uniform(0.5, 3)))),
    "cold brown dwarf proxy": lambda r: Spectrum("bd", "browndwarf", (float(r.uniform(350, 700)),)),
}


def wilson(k, n, z=1.0):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def m50(mags, comp):
    mags, comp = np.asarray(mags), np.asarray(comp)
    for i in range(len(mags) - 1):
        if comp[i] >= 0.5 > comp[i + 1]:
            f = (comp[i] - 0.5) / (comp[i] - comp[i + 1])
            return float(mags[i] + f * (mags[i + 1] - mags[i]))
    return float("nan")


def base_cfg(seed: int, size: int) -> SimConfig:
    return SimConfig(size=size, seed=seed, spherex_bursts=True, n_bursts=20, span_days=20.0,
                     star_density_per_deg2=0.0, mask_frac=0.02)


def run(n_null: int, per_mag: int, mags: list[float], size: int, vmax: float, fap: float, out: Path,
        reuse_null: Path | None = None) -> dict:
    t_start = time.time()
    # ---------------------------------------------------------------- nulls
    null_max: dict[str, list[float]] = {}
    null_tracklets = []
    t_search = []
    if reuse_null is not None and reuse_null.exists():
        prev = json.loads(reuse_null.read_text())
        null_max = prev["null_max"]
        n_null = len(next(iter(null_max.values())))
        null_tracklets = [prev.get("null_false_tracklets_mean", 0.0)] * n_null
        t_search = [prev.get("search_seconds_mean", np.nan)]
        print(f"reusing {n_null} null fields from {reuse_null}", flush=True)
    for i in range(0 if null_max else n_null):
        frames, _ = simulate(base_cfg(1000 + i, size), [], include_stars=False)
        maps = build_maps(frames, q=2)
        grid = make_grid(frames, half_fov_arcsec=0.5 * size * 6.15 * 0.8, vmax_arcsec_day=vmax)
        t0 = time.time()
        mx = multi_search(maps, grid, spectrum_true=None)
        t_search.append(time.time() - t0)
        for k, v in mx.items():
            null_max.setdefault(k, []).append(v)
        null_tracklets.append(link_false_tracklets(maps, vmax_arcsec_day=vmax))
        print(f"null {i + 1}/{n_null}  {time.time() - t0:.1f}s  " +
              " ".join(f"{k}={v:.2f}" for k, v in mx.items()), flush=True)
    thresholds = {k: float(np.quantile(v, 1 - fap)) for k, v in null_max.items()}
    # the oracle template statistic is calibrated with the solar-type bank member
    # null distribution (same Gaussian null); use the 'stack' null as proxy
    thresholds.setdefault("template", thresholds["stack"])
    print("thresholds", thresholds, flush=True)
    # ------------------------------------------------------------ injections
    rng = np.random.default_rng(5000)
    rec = []
    for im, mag in enumerate(mags):
        for j in range(per_mag):
            seed = 5000 + im * 1000 + j
            r = np.random.default_rng(seed)
            cls = list(SED_CLASSES)[j % len(SED_CLASSES)]
            spec = SED_CLASSES[cls](r)
            cfg = base_cfg(seed, size)
            lim = 0.25 * size * 6.15
            vx, vy = r.uniform(-vmax * 0.9, vmax * 0.9, 2)
            frames0, truth0 = simulate(cfg, [], include_stars=False)
            mv = MovingSource(float(r.uniform(-lim, lim)), float(r.uniform(-lim, lim)), float(vx), float(vy),
                              float(mag), spec, t_ref=truth0["t_ref"])
            frames, truth = simulate(cfg, [mv], include_stars=False)
            maps = build_maps(frames, q=2)
            grid = make_grid(frames, half_fov_arcsec=0.5 * size * 6.15 * 0.8, vmax_arcsec_day=vmax)
            # truth snapped to the grid (what an exhaustive grid search evaluates)
            ps = grid.xi[1] - grid.xi[0]
            vs = grid.vxi[1] - grid.vxi[0]
            t_ref = grid.t_ref
            xi_t, eta_t = mv.pos(t_ref)
            sx, sy = np.round(xi_t / ps) * ps, np.round(eta_t / ps) * ps
            svx, svy = np.round(vx / vs) * vs, np.round(vy / vs) * vs
            st = local_stats(maps, sx, sy, svx, svy, t_ref, spectrum_true=spec, pos_step=ps, v_step=vs)
            nlink = link_count_at_truth(maps, xi_t, eta_t, vx, vy, t_ref)
            det = {k: bool(v >= thresholds[k]) for k, v in st.items()}
            det["link"] = bool(nlink >= 3)
            # analytic expectation: SNR of the constant-flux statistic relative to
            # the matched (true-SED) statistic, R = sum(s phi)/sqrt(sum phi sum s^2 phi)
            from ..inference.matched import Trajectories, gather, spectral_factors
            trj = Trajectories(np.array([xi_t]), np.array([eta_t]), np.array([vx]), np.array([vy]), t_ref)
            _, phi_t = gather(maps, trj)
            sf = spectral_factors(frames, spec)
            ph = phi_t[:, 0]
            R = float(np.sum(sf * ph) / np.sqrt(np.sum(ph) * np.sum(sf * sf * ph))) if np.sum(ph) > 0 else np.nan
            rec.append({"mag": mag, "class": cls, "stats": st, "detected": det, "n_link_exp": nlink,
                        "v": [float(vx), float(vy)], "n_frames": len(frames), "mismatch_R": R})
        comp = {k: np.mean([x["detected"][k] for x in rec if x["mag"] == mag]) for k in rec[-1]["detected"]}
        print(f"mag {mag}: " + " ".join(f"{k}={v:.2f}" for k, v in comp.items()), flush=True)
    # ------------------------------------------------------------- summary
    methods = list(rec[0]["detected"])
    summary = {}
    for k in methods:
        comp, lo, hi = [], [], []
        for mag in mags:
            sel = [x["detected"][k] for x in rec if x["mag"] == mag]
            comp.append(float(np.mean(sel)))
            a, b = wilson(sum(sel), len(sel))
            lo.append(a)
            hi.append(b)
        by_class = {}
        for cls in SED_CLASSES:
            cc = [float(np.mean([x["detected"][k] for x in rec if x["mag"] == mag and x["class"] == cls]))
                  for mag in mags]
            by_class[cls] = {"completeness": cc, "m50": m50(mags, cc)}
        summary[k] = {"completeness": comp, "lo68": lo, "hi68": hi, "m50": m50(mags, comp), "by_class": by_class,
                      "threshold": thresholds.get(k)}
    theory = {}
    for cls in SED_CLASSES:
        R = np.array([x["mismatch_R"] for x in rec if x["class"] == cls and np.isfinite(x["mismatch_R"])])
        theory[cls] = {"mean_R_constant_vs_matched": float(np.mean(R)) if R.size else None,
                       "expected_dm50_matched_minus_constant": float(2.5 * np.log10(1 / np.mean(R))) if R.size else None}
    res = {
        "theory_spectral_mismatch": theory,
        "experiment": "E1 synthetic injection-recovery, SPHEREx burst cadence, noise-limited fields",
        "provenance": "synthetic",
        "config": {"n_null": n_null, "per_mag": per_mag, "mags": mags, "size_px": size, "vmax_arcsec_day": vmax,
                   "fap_per_field": fap, "n_bursts": 20, "span_days": 20.0,
                   "noise_sigma_ujy": base_cfg(0, size).sigma_ujy, "psf": "pixel-integrated Gaussian FWHM 1.22 px",
                   "link_rule": ">=3 exposures with per-frame matched-filter SNR>=5 within 1.5 px of truth"},
        "null_max": null_max, "thresholds": thresholds,
        "null_false_tracklets_mean": float(np.mean(null_tracklets)),
        "null_fields_with_false_tracklet": float(np.mean(np.array(null_tracklets) > 0)),
        "search_seconds_mean": float(np.mean(t_search)), "summary": summary, "records": rec,
        "runtime_s": time.time() - t_start,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=float))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", type=int, default=40)
    ap.add_argument("--per-mag", type=int, default=21)
    ap.add_argument("--mags", type=str, default="19.5,20,20.5,21,21.5,22,22.5")
    ap.add_argument("--size", type=int, default=24)
    ap.add_argument("--vmax", type=float, default=3.0)
    ap.add_argument("--fap", type=float, default=0.05)
    ap.add_argument("--out", type=str, default=str(ROOT / "reports" / "motion_benchmark.json"))
    ap.add_argument("--reuse-null", type=str, default=None)
    a = ap.parse_args()
    mags = [float(m) for m in a.mags.split(",")] if "," in a.mags else list(np.arange(*map(float, a.mags.split(":"))))
    run(a.null, a.per_mag, mags, a.size, a.vmax, a.fap, Path(a.out),
        Path(a.reuse_null) if a.reuse_null else None)


if __name__ == "__main__":
    main()
