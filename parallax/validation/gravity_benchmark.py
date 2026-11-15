"""E3: Gravity-mode injection-recovery (synthetic, labelled).

Question: does a joint spectro-temporal point-lens fit with *chromatic*
blending detect and characterise lensing better than (a) an achromatic
("grey") multi-band fit and (b) a classic single-band fit, at the same
null-calibrated false-alarm probability?  And how often does an ordinary
symmetric variable masquerade as a lens?

Methods (all fitted to exactly the same light curves):

* joint  : F = A(t) S(lambda) + Bl(lambda), S and Bl in a smooth spectral basis
           (GravityModel, k=3).  Statistic: delta-chi2(H1 vs H0).
* grey   : F = F0(lambda) [1 + a (A(t) - 1)],  F0 = best static SED (2K basis);
           one blend fraction for all wavelengths (what one gets by
           normalising every band and fitting a single light curve).
* single : only samples with 1.6 < lambda < 2.6 um (one broad band),
           F = f_b(lambda) + f_s(lambda) A(t), each linear in lambda across the
           band  -- classic single-filter PSPL.

Thresholds: (1 - FAP) quantile of each method's statistic on null (constant)
light curves with the same cadence and noise.  Seeds: nulls 100000+, lenses
200000+, variables 300000+.

Usage: python -m parallax.validation.gravity_benchmark --cadence deep
"""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from ..gravity.pspl import GravityModel, LightCurve, _wls, magnification, simulate_event, u_of_t

ROOT = Path(__file__).resolve().parents[2]
U0_VALUES = [1.0, 1.5, 2.0, 2.5, 3.0]


def cadence(kind: str, rng):
    if kind == "deep":
        t = np.sort(rng.uniform(0, 200, 330))
    else:   # wide: two ~20-d visits 180 d apart, 25 bursts of 4 exposures each
        bursts = np.concatenate([rng.uniform(0, 20, 25), rng.uniform(180, 200, 25)])
        t = np.sort(np.concatenate([b + np.arange(4) * 0.0015 for b in bursts]))
    wave = np.exp(rng.uniform(np.log(0.75), np.log(5.0), t.size))
    return t, wave


def _scan(lc: LightCurve, design, coarse=(40, 20)):
    t = lc.t
    w = 1.0 / lc.err ** 2
    best = (np.inf, None)
    for t0 in np.linspace(t.min() - 10, t.max() + 10, coarse[0]):
        for tE in np.geomspace(1.0, 150, coarse[1]):
            for u0 in (0.05, 0.2, 0.5, 1.0, 1.6, 2.4):
                X = design(magnification(u_of_t(t, t0, u0, tE)))
                x2 = _wls(X, lc.flux, w)[1]
                if x2 < best[0]:
                    best = (x2, (t0, np.log(u0), np.log(tE)))
    f = lambda p: _wls(design(magnification(u_of_t(t, p[0], np.exp(p[1]), np.exp(p[2])))), lc.flux, w)[1]
    r = minimize(f, np.array(best[1]), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3, "maxiter": 1500})
    return float(min(r.fun, best[0])), {"t0": float(r.x[0]), "u0": float(np.exp(r.x[1])), "tE": float(np.exp(r.x[2]))}


def stats(lc: LightCurve) -> dict:
    out = {}
    gm = GravityModel(lc, k_spec=3)
    w = gm.w
    h0 = gm.fit_constant()
    # joint chromatic (coarse scan with the same generic scanner for fairness)
    x2j, pj = _scan(lc, lambda A: np.hstack([gm.Bs * A[:, None], gm.B2]))
    out["joint"] = {"dchi2": h0.chi2 - x2j, **pj}
    # grey: one blend fraction for every wavelength
    F0 = h0.model_flux
    x2g, pg = _scan(lc, lambda A: np.hstack([gm.B2, ((A - 1.0) * F0)[:, None]]))
    out["grey"] = {"dchi2": h0.chi2 - x2g, **pg}
    # single band
    sel = (lc.wave > 1.6) & (lc.wave < 2.6)
    if sel.sum() >= 8:
        sub = LightCurve(lc.t[sel], lc.flux[sel], lc.err[sel], lc.wave[sel])
        ws = 1.0 / sub.err ** 2
        # in-band SED slope: linear in (lambda - 2um) for both the blend and the source
        dl = sub.wave - 2.0
        x20 = _wls(np.stack([np.ones_like(dl), dl], axis=1), sub.flux, ws)[1]
        x2s, ps = _scan(sub, lambda A: np.stack([np.ones_like(A), dl, A, A * dl], axis=1))
        out["single"] = {"dchi2": x20 - x2s, "n": int(sel.sum()), **ps}
    else:
        out["single"] = {"dchi2": 0.0, "n": int(sel.sum())}
    # hypothesis comparison (lens vs symmetric variability) for the joint model
    h2 = gm.fit_variability()
    n = lc.t.size
    out["dBIC_lens_vs_var"] = float((h2.chi2 + h2.n_par * np.log(n)) - (x2j + (3 * gm.k + 3) * np.log(n)))
    out["dBIC_lens_vs_const"] = float((h0.chi2 + h0.n_par * np.log(n)) - (x2j + (3 * gm.k + 3) * np.log(n)))
    # the app's full decision (adequacy + physical-validity guards + BIC)
    out["app_verdict"] = gm.compare()["verdict"]
    return out


def _variable(t, wave, rng, sigma):
    """Symmetric achromatic brightening (e.g. an outburst / spot-modulation
    proxy) with the same SED machinery: F = S(lambda) [1 + a exp(-(t-c)^2/2w^2)]."""
    from ..science.spectra import ab_to_ujy, planck_nu
    S = ab_to_ujy(16.0) * planck_nu(wave, 4500.0) / planck_nu(2.0, 4500.0)
    c, wd, a = rng.uniform(t.min(), t.max()), np.exp(rng.uniform(np.log(3), np.log(40))), rng.uniform(0.1, 0.6)
    true = S * (1 + a * np.exp(-0.5 * ((t - c) / wd) ** 2))
    err = sigma * S
    return LightCurve(t, true + rng.normal(0, err), err, wave, provenance="synthetic",
                      meta={"centre": c, "width": wd, "amp": a})


def job(args):
    kind, seed, cad, sigma, u0 = args
    rng = np.random.default_rng(seed)
    t, wave = cadence(cad, rng)
    if kind == "var":
        lc = _variable(t, wave, rng, sigma)
        meta = lc.meta
    else:
        tE = float(np.exp(rng.uniform(np.log(5), np.log(60))))
        t0 = float(rng.uniform(t.min() - 0.5 * tE, t.max() + 0.5 * tE))
        blend = float(rng.uniform(0.3, 0.8))
        lc = simulate_event(t, wave, t0=t0, u0=u0 if kind == "lens" else 1.0, tE=tE, blend_frac_2um=blend,
                            src_T=6000.0, blend_T=3200.0, sigma_frac=sigma, seed=seed + 1, lensed=(kind == "lens"))
        meta = {"t0": t0, "u0": u0, "tE": tE, "blend": blend}
    return {"kind": kind, "seed": seed, "u0": u0, "meta": meta, "stats": stats(lc)}


def run(cad: str, n_null: int, per_u0: int, n_var: int, sigma: float, fap: float, workers: int, out: Path):
    t0 = time.time()
    tasks = [("null", 100000 + i, cad, sigma, None) for i in range(n_null)]
    tasks += [("lens", 200000 + 1000 * k + j, cad, sigma, u) for k, u in enumerate(U0_VALUES) for j in range(per_u0)]
    tasks += [("var", 300000 + i, cad, sigma, None) for i in range(n_var)]
    with ProcessPoolExecutor(workers) as ex:
        rec = []
        for i, r in enumerate(ex.map(job, tasks, chunksize=2)):
            rec.append(r)
            if i % 20 == 0:
                print(f"{i + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    methods = ["joint", "grey", "single"]
    nulls = [r for r in rec if r["kind"] == "null"]
    thr = {m: float(np.quantile([r["stats"][m]["dchi2"] for r in nulls], 1 - fap)) for m in methods}
    eff = {m: [] for m in methods}
    te_err = {m: [] for m in methods}
    for u in U0_VALUES:
        L = [r for r in rec if r["kind"] == "lens" and r["u0"] == u]
        for m in methods:
            det = [r["stats"][m]["dchi2"] >= thr[m] for r in L]
            eff[m].append(float(np.mean(det)))
            rel = [abs(np.log(r["stats"][m].get("tE", np.nan) / r["meta"]["tE"])) for r, d in zip(L, det) if d]
            te_err[m].append(float(np.median(rel)) if rel else None)
    lenses = [r for r in rec if r["kind"] == "lens"]
    vars_ = [r for r in rec if r["kind"] == "var"]

    def verdict(r):
        s = r["stats"]
        if s["dBIC_lens_vs_const"] < 10:
            return "none"
        return "lens" if s["dBIC_lens_vs_var"] > 10 else ("var" if s["dBIC_lens_vs_var"] < -10 else "ambiguous")

    def app(r):
        v = r["stats"]["app_verdict"]
        return ("lens" if v.startswith("Microlensing") else "var" if v.startswith("Alternative") else
                "none" if v.startswith("Insufficient") else "inadequate" if v.startswith("Model inadequate") else
                "ambiguous")

    vl = [app(r) for r in lenses]
    vv = [app(r) for r in vars_]
    vn = [app(r) for r in nulls]
    bic_only = {"lens_injections": [verdict(r) for r in lenses], "variable_injections": [verdict(r) for r in vars_],
                "null_injections": [verdict(r) for r in nulls]}
    res = {
        "experiment": f"E3 gravity injection-recovery ({cad} cadence)", "provenance": "synthetic", "cadence": cad,
        "description": (f"{n_null} null, {per_u0} lens events per u0 and {n_var} symmetric achromatic variables, "
                        f"{'330 samples over 200 d' if cad == 'deep' else 'two 20-day visits of 25 bursts x 4 exposures'}, "
                        f"{sigma * 100:.0f}% per-point noise, random wavelength per sample (0.75-5 um), source 6000 K, "
                        "blend 3200 K, blend fraction 0.3-0.8 at 2 um, tE log-uniform 5-60 d, t0 uniform over the span."),
        "fap": fap, "u0_values": U0_VALUES, "thresholds": thr,
        "methods": {m: {"label": {"joint": "Joint chromatic (PARALLAX)", "grey": "Grey multi-band",
                                  "single": "Single band 1.6-2.6 um"}[m],
                        "threshold": thr[m], "efficiency": eff[m], "median_abs_ln_tE_error": te_err[m]} for m in methods},
        "classification": {
            "lens_injections": {k: vl.count(k) for k in ("lens", "var", "ambiguous", "none", "inadequate")},
            "variable_injections": {k: vv.count(k) for k in ("lens", "var", "ambiguous", "none", "inadequate")},
            "null_injections": {k: vn.count(k) for k in ("lens", "var", "ambiguous", "none", "inadequate")},
            "rule": "app decision: adequacy (reduced chi2 <= 3), physical lens (0<f_s<=1.05, tE<=3 span), BIC 10",
        },
        "classification_bic_only": {k: {c: v.count(c) for c in ("lens", "var", "ambiguous", "none")}
                                    for k, v in bic_only.items()},
        "records": rec, "runtime_s": time.time() - t0,
    }
    res["conclusion"] = conclude(res)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1, default=float))
    print(res["conclusion"], flush=True)
    return res


def conclude(res: dict) -> str:
    """Plain-language summary generated only from the stored numbers."""
    e, u = res["methods"], res["u0_values"]
    gain = [a - b for a, b in zip(e["joint"]["efficiency"], e["grey"]["efficiency"])]
    single = [a - b for a, b in zip(e["joint"]["efficiency"], e["single"]["efficiency"])]
    cl = res["classification"]
    nv = sum(cl["variable_injections"].values())
    nn = sum(cl["null_injections"].values())
    tj, tg = e["joint"]["median_abs_ln_tE_error"][0], e["grey"]["median_abs_ln_tE_error"][0]
    txt = (f"Detection at FAP {res['fap']} per light curve, efficiency difference by u0 {u}: joint - grey "
           + ", ".join(f"{g:+.2f}" for g in gain) + "; joint - single band " + ", ".join(f"{g:+.2f}" for g in single) + ". ")
    if tj is not None and tg is not None:
        txt += (f"Timescale accuracy for detected u0={u[0]} events: median |ln(tE_fit/tE_true)| {tj:.2f} (joint) vs "
                f"{tg:.2f} (grey). ")
    txt += (f"Symmetric achromatic variables given a lens verdict: {cl['variable_injections']['lens']}/{nv}; "
            f"nulls given any lens verdict: {cl['null_injections']['lens']}/{nn}")
    bo = res.get("classification_bic_only")
    if bo:
        txt += f" (BIC alone, without the physical-validity guard: {bo['variable_injections']['lens']}/{nv} variables)"
    txt += "."
    return txt


def resummarize():
    for c in ("deep", "wide"):
        p = ROOT / "reports" / f"gravity_benchmark_{c}.json"
        if p.exists():
            d = json.loads(p.read_text())
            d["conclusion"] = conclude(d)
            p.write_text(json.dumps(d, indent=1, default=float))
            print(c, d["conclusion"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cadence", default="deep", choices=["deep", "wide"])
    ap.add_argument("--null", type=int, default=300)
    ap.add_argument("--per-u0", type=int, default=60)
    ap.add_argument("--var", type=int, default=100)
    ap.add_argument("--sigma", type=float, default=0.02)
    ap.add_argument("--fap", type=float, default=0.01)
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--out", default=None)
    ap.add_argument("--resummarize", action="store_true")
    a = ap.parse_args()
    if a.resummarize:
        return resummarize()
    out = Path(a.out) if a.out else ROOT / "reports" / f"gravity_benchmark_{a.cadence}.json"
    run(a.cadence, a.null, a.per_u0, a.var, a.sigma, a.fap, a.workers, out)


if __name__ == "__main__":
    main()
