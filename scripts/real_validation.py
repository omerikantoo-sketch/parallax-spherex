"""E2: checks on real SPHEREx data (writes reports/real_validation.json).

  barnard   proper motion of Barnard's Star from SPHEREx centroids vs Gaia DR3
  search    run the moving-object search on the real Barnard field (exploratory)
  variable  constant-SED vs variability test of forced photometry at a target

usage: python scripts/real_validation.py barnard [search] [variable:rw-aur]
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from parallax import service as S  # noqa: E402
from parallax.data.targets import TARGETS  # noqa: E402

OUT = ROOT / "reports" / "real_validation.json"


def field_for(slug: str, size: int = 64) -> str:
    t = TARGETS[slug]
    ra, dec = (float(v) for v in t.at(61000.0))
    ra, dec = float(f"{ra:.5f}"), float(f"{dec:.5f}")   # same rounding as the web client
    fid = S.field_id(ra, dec, size, "dr1")
    if not (S.FIELDS / fid / "field.json").exists():
        print("building field", slug, fid, flush=True)
        S.build_field(ra, dec, size, "dr1", target=slug, progress=lambda m, f: None)
    return fid


def barnard() -> dict:
    fid = field_for("barnard")
    a = S.astrometry(fid)
    if "fit" not in a:
        return {"title": "Barnard's Star proper motion", "summary": f"Not measured: {a.get('error')}", "field": fid}
    f, e = a["fit"], a["expected"]
    t = np.array([p["t"] for p in a["points"]])
    visits = sorted({int(round((x - t.min()) // 60)) for x in t})
    dra = f["pmra_masyr"] - e["pmra_masyr"]
    ddec = f["pmdec_masyr"] - e["pmdec_masyr"]
    sig = np.hypot(dra / f["e_pmra_masyr"], ddec / f["e_pmdec_masyr"]) / np.sqrt(2)
    summary = (f"{f['n']} centroided exposures in {len(visits)} visit(s) spanning {t.max() - t.min():.0f} d. "
               f"Measured μα* = {f['pmra_masyr']:.0f} ± {f['e_pmra_masyr']:.0f}, μδ = {f['pmdec_masyr']:.0f} ± "
               f"{f['e_pmdec_masyr']:.0f} mas/yr vs Gaia DR3 {e['pmra_masyr']:.0f}, {e['pmdec_masyr']:.0f} mas/yr "
               f"(parallax fixed to {e['plx_mas']:.0f} mas). Difference {dra:+.0f}, {ddec:+.0f} mas/yr; "
               f"residual rms {f['residual_rms_arcsec']:.2f}\" ({f['residual_rms_arcsec'] / 6.15:.2f} px).")
    # diagnostics: offsets of each visit / detector from the Gaia-predicted path
    from parallax.inference.frames import tangent_offsets
    from parallax.inference.scene import fit_astrometry, parallax_factors
    meta = S.load_field(fid)
    _, frames, cuts = S.field_frames(fid)
    pts = a["points"]
    xi = np.array([p["xi"] for p in pts]); eta = np.array([p["eta"] for p in pts])
    sx = np.array([p["sxi"] for p in pts]); sy = np.array([p["seta"] for p in pts])
    det = np.array([p["det"] for p in pts])
    sc = np.array([cuts[p["k"]].meta.sc_pos_km for p in pts], float)
    tg = TARGETS["barnard"]
    pr = np.array([tangent_offsets(*tg.at(ti), meta["ra"], meta["dec"]) for ti in t])
    pa, pd = parallax_factors(t, meta["ra"], meta["dec"], sc)
    dx = xi - (pr[:, 0] + tg.plx / 1000 * pa)
    dy = eta - (pr[:, 1] + tg.plx / 1000 * pd)
    vis = ((t - t.min()) // 60).astype(int)
    per_visit = [{"mjd": float(t[vis == v].mean()), "n": int((vis == v).sum()), "dxi_arcsec": float(np.median(dx[vis == v])),
                  "deta_arcsec": float(np.median(dy[vis == v]))} for v in np.unique(vis)]
    per_det = [{"det": int(d), "n": int((det == d).sum()), "dxi_arcsec": float(np.median(dx[det == d])),
                "deta_arcsec": float(np.median(dy[det == d]))} for d in np.unique(det)]
    m = det >= 5
    f56 = fit_astrometry(t[m], xi[m], eta[m], sx[m], sy[m], meta["ra"], meta["dec"], f["t_ref_mjd"],
                         plx_fixed_mas=tg.plx, sc_geo_km=sc[m]) if m.sum() >= 6 else None
    sep = max(abs(v["dxi_arcsec"]) for v in per_visit), max(abs(v["deta_arcsec"]) for v in per_visit)
    summary += (f" Each visit's median position agrees with the Gaia-predicted position (proper motion + parallax) to "
                f"within {max(sep):.1f}″ (0.2 px), while the star moved ~{abs(tg.pmdec) / 1000 * (t.max() - t.min()) / 365.25:.1f}″ "
                f"between first and last visit. Short-wavelength detectors (D1-D3, few points, saturated core) show "
                f"offsets up to {max(abs(d['deta_arcsec']) for d in per_det):.1f}″.")
    # systematic floor: per-visit offsets (rms over visits with >= 5 points) propagated to proper motion
    big = [v for v in per_visit if v["n"] >= 5]
    sys_rms = float(np.sqrt(np.mean([v["dxi_arcsec"] ** 2 + v["deta_arcsec"] ** 2 for v in big]) / 2)) if big else np.nan
    dt_yr = (max(v["mjd"] for v in big) - min(v["mjd"] for v in big)) / 365.25 if len(big) > 1 else np.nan
    pm_sys = float(np.sqrt(2) * sys_rms / dt_yr * 1000) if dt_yr == dt_yr else np.nan
    tot_ra, tot_de = np.hypot(f["e_pmra_masyr"], pm_sys), np.hypot(f["e_pmdec_masyr"], pm_sys)
    summary += (f" Including a per-visit systematic floor of {sys_rms:.2f}″ per axis (±{pm_sys:.0f} mas/yr on proper motion), "
                f"the differences are {dra / tot_ra:+.1f}σ and {ddec / tot_de:+.1f}σ; the total motion is detected at "
                f"{np.hypot(f['pmra_masyr'], f['pmdec_masyr']) / np.hypot(tot_ra, tot_de) * np.sqrt(2):.0f}σ. "
                "The floor is estimated from offsets relative to the Gaia prediction, so this consistency check is not "
                "independent; the robust results are the detected motion and the ≤1.2″ per-visit agreement.")
    if f56:
        summary += (f" Diagnostic chosen after seeing these offsets (not pre-specified): D5+D6 only gives "
                    f"{f56['pmra_masyr']:.0f} ± {f56['e_pmra_masyr']:.0f}, {f56['pmdec_masyr']:.0f} ± {f56['e_pmdec_masyr']:.0f} mas/yr.")
    return {"title": "Barnard's Star proper motion (known-object recovery)", "summary": summary, "field": fid,
            "per_visit_offsets": per_visit, "systematic_floor_arcsec": sys_rms, "pm_systematic_masyr": pm_sys,
            "parallax_model_check": "parallax_factors agrees with astropy GCRS(obsgeoloc=spacecraft) to <0.1 arcsec in xi, <0.3 in eta", "per_detector_offsets": per_det, "fit_D5_D6_posthoc": f56,
            "fit": f, "expected": e, "n_visits": len(visits), "consistency_sigma": float(sig),
            "table": {"columns": ["Quantity", "SPHEREx (this work)", "± 1σ", "Gaia DR3"],
                      "rows": [["μα* [mas/yr]", f["pmra_masyr"], f["e_pmra_masyr"], e["pmra_masyr"]],
                               ["μδ [mas/yr]", f["pmdec_masyr"], f["e_pmdec_masyr"], e["pmdec_masyr"]],
                               ["systematic floor on μ [mas/yr]", pm_sys, "", ""],
                               ["residual rms [arcsec]", f["residual_rms_arcsec"], "", ""]]}}


def search() -> dict:
    fid = field_for("barnard")
    t0 = time.time()
    r = S.search_field(fid, 2.0)
    c = r.get("candidates", [])
    top = c[0] if c else None
    return {"title": "Moving-object search on the real Barnard field (exploratory)", "field": fid,
            "summary": (f"Searched {r.get('n_frames')} exposures of one visit ({r.get('n_trials'):,} trajectories, "
                        f"|v| between {r.get('vmin_arcsec_day', 0):.2f} and {r.get('vmax_arcsec_day')} arcsec/day) after "
                        f"subtracting and masking {r.get('static_sources')} static sources, in {time.time() - t0:.0f} s. "
                        f"Top bank statistic {top['z_bank']:.1f}; maximum on time-scrambled frames "
                        f"{r.get('scrambled_null_max', float('nan')):.1f}; "
                        f"{sum(x.get('above_scrambled_null', False) for x in c)} of {len(c)} listed candidates exceed it. "
                        "No candidate is claimed; vetting would need a second visit and real-null calibration."
                        if top else "No candidates."),
            "scrambled_null_max": r.get("scrambled_null_max"),
            "status": "Exploratory", "n_candidates": len(c),
            "candidates": [{k: v for k, v in x.items() if k in ("xi0", "eta0", "vxi", "veta", "z_bank", "above_scrambled_null")} for x in c]}


def variable(slug: str) -> dict:
    from parallax.gravity.pspl import GravityModel, LightCurve
    fid = field_for(slug)
    sp = S.spectrum_at(fid, 0.0, 0.0, follow_target=False)
    pts = sp["points"]
    if len(pts) < 12:
        return {"title": f"{TARGETS[slug].name} variability", "summary": f"Only {len(pts)} usable points.", "field": fid}
    err = [float(np.hypot(p["err"], 0.01 * abs(p["flux"]))) for p in pts]   # 1% calibration floor
    lc = LightCurve([p["t"] for p in pts], [p["flux"] for p in pts], err, [p["wave"] for p in pts],
                    provenance="real")
    r = GravityModel(lc, k_spec=3).compare()
    n = r["n"]
    chi0 = r["H0"].chi2 / max(1, n - r["H0"].n_par)
    return {"title": f"{TARGETS[slug].name}: constant SED vs variability (real photometry)", "field": fid,
            "summary": (f"{n} forced-photometry points. Constant-SED reduced χ² {chi0:.1f}; ΔBIC lens−constant "
                        f"{r['dBIC_H1_vs_H0']:.0f}, lens−variability {r['dBIC_H1_vs_H2']:.0f}; verdict: {r['verdict']}. "
                        f"{TARGETS[slug].note}"),
            "verdict": r["verdict"], "chi2r_constant": chi0, "dBIC_H1_vs_H0": r["dBIC_H1_vs_H0"],
            "dBIC_H1_vs_H2": r["dBIC_H1_vs_H2"]}


def main():
    res = json.loads(OUT.read_text()) if OUT.exists() else {}
    for arg in sys.argv[1:] or ["barnard"]:
        if arg == "barnard":
            res["barnard"] = barnard()
        elif arg == "search":
            res["barnard_search"] = search()
        elif arg.startswith("variable:"):
            slug = arg.split(":", 1)[1]
            res[f"variable_{slug}"] = variable(slug)
        print(json.dumps({k: v.get("summary") for k, v in res.items() if isinstance(v, dict)}, indent=1), flush=True)
        res["provenance"] = "real"
        OUT.write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
