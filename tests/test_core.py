"""Unit and integration tests for PARALLAX.  Offline except where marked
(tests needing cached real cutouts are skipped when the cache is empty)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from parallax.data.model import FLAG_BITS, REJECT_MASK  # noqa: E402
from parallax.gravity.pspl import GravityModel, magnification, mass_posterior, simulate_event  # noqa: E402
from parallax.inference.frames import tangent_offsets, tangent_to_sky  # noqa: E402
from parallax.inference.matched import (Trajectories, branch_and_bound, build_maps, gather, make_grid,  # noqa: E402
                                        spectral_factors, stat_constant, stat_template)
from parallax.science.psf import EPSF, GaussianPSF  # noqa: E402
from parallax.science.simulate import MovingSource, SimConfig, simulate  # noqa: E402
from parallax.science.spectra import Spectrum, ab_to_ujy, spectral_basis, ujy_to_ab  # noqa: E402
from parallax.validation.motion_benchmark import m50, wilson  # noqa: E402

CUTOUTS = sorted((ROOT / "cache" / "cutouts").glob("*.npz"))


# ------------------------------------------------------------------ units
def test_ab_roundtrip_and_zero_point():
    assert ab_to_ujy(23.9) == pytest.approx(1.0, rel=1e-6)          # AB 23.9 = 1 uJy
    m = np.array([15.0, 20.5, 24.0])
    assert np.allclose(ujy_to_ab(ab_to_ujy(m)), m)


def test_reject_mask_contains_required_bits():
    for name in ("TRANSIENT", "OVERFLOW", "NONFUNC", "MISSING_DATA", "HOT", "PERSIST", "OUTLIER"):
        assert REJECT_MASK & (1 << FLAG_BITS[name])
    assert not REJECT_MASK & (1 << FLAG_BITS["SOURCE"])             # sources are not bad pixels


# ------------------------------------------------------------------ PSF
def test_gaussian_psf_normalised_and_centred():
    p = GaussianPSF(1.22)
    img = p.render((31, 31), 15.3, 14.6)
    assert img.sum() == pytest.approx(1.0, abs=1e-6)
    yy, xx = np.mgrid[:31, :31]
    # first moment of an undersampled (FWHM 1.22 px) pixel-integrated PSF is
    # biased by a few millipixels (pixel-phase effect), hence 0.01 px tolerance
    assert (img * xx).sum() == pytest.approx(15.3, abs=1e-2)
    assert (img * yy).sum() == pytest.approx(14.6, abs=1e-2)


@pytest.mark.skipif(not CUTOUTS, reason="no cached real cutouts")
def test_real_epsf_native_sum_is_one():
    from parallax.data.archive import load_cutout
    c = load_cutout(CUTOUTS[0])
    p = EPSF(c.psf_oversampled, c.psf_oversamp)
    img = p.render((41, 41), 20.2, 19.7)
    assert img.sum() == pytest.approx(1.0, abs=0.02)
    assert 0.8 < p.fwhm_px() < 2.5                                   # ~7.5" / 6.15"/px


@pytest.mark.skipif(not CUTOUTS, reason="no cached real cutouts")
def test_real_cutout_units_and_wavelength_gradient():
    from parallax.data.archive import load_cutout
    c = load_cutout(CUTOUTS[0])
    f, v, g = c.flux_ujy()
    k = c.mjysr_to_ujy
    # zodi is subtracted exactly once
    assert np.allclose(f[g], ((c.image_mjysr - c.zodi_mjysr) * k)[g], equal_nan=True)
    assert np.all(v[g] > 0)
    # 1 MJy/sr over ~37.9 arcsec^2 is ~890 uJy
    assert 700 < float(np.median(k)) < 1100
    # linear variable filter: wavelength varies smoothly and monotonically along one axis
    w = c.wave_um
    gy, gx = np.nanmean(np.abs(np.diff(w, axis=0))), np.nanmean(np.abs(np.diff(w, axis=1)))
    assert max(gy, gx) > 5 * min(gy, gx)
    lo, hi = {1: (0.73, 1.12), 2: (1.08, 1.65), 3: (1.6, 2.45), 4: (2.4, 3.85), 5: (3.8, 4.45), 6: (4.4, 5.05)}[c.meta.detector]
    assert lo < np.nanmedian(w) < hi


# --------------------------------------------------------------- geometry
def test_tangent_roundtrip():
    ra0, dec0 = 269.45, 4.74
    ra, dec = np.array([269.46, 269.40]), np.array([4.75, 4.70])
    xi, eta = tangent_offsets(ra, dec, ra0, dec0)
    r2, d2 = tangent_to_sky(xi, eta, ra0, dec0)
    assert np.allclose(r2, ra, atol=1e-9) and np.allclose(d2, dec, atol=1e-9)
    # east is +xi
    assert xi[0] > 0


def test_spectral_basis_shape_and_rank():
    w = np.geomspace(0.75, 5.0, 200)
    B = spectral_basis(w, 4)
    assert B.shape == (200, 4) and np.linalg.matrix_rank(B) == 4


# --------------------------------------------------------------- inference
def _small_field(mag=19.0, seed=3):
    cfg = SimConfig(size=16, seed=seed, spherex_bursts=True, n_bursts=6, span_days=6.0,
                    star_density_per_deg2=0.0, mask_frac=0.0)
    f0, t0 = simulate(cfg, [], include_stars=False)
    mv = MovingSource(5.0, -4.0, 1.2, 0.6, mag, Spectrum("bd", "browndwarf", (500.0,)), t_ref=t0["t_ref"])
    frames, truth = simulate(cfg, [mv], include_stars=False)
    return frames, mv


def test_matched_filter_peaks_at_truth():
    frames, mv = _small_field(mag=17.5)
    maps = build_maps(frames, q=2)
    grid = make_grid(frames, half_fov_arcsec=30, vmax_arcsec_day=2)
    s = spectral_factors(frames, mv.spectrum)
    xi, eta = mv.pos(grid.t_ref)
    tr = Trajectories(np.array([xi, xi + 18]), np.array([eta, eta - 18]), np.array([mv.vxi, mv.vxi]),
                      np.array([mv.veta, mv.veta]), grid.t_ref)
    psi, phi = gather(maps, tr)
    z_true, _ = stat_template(psi, phi, s)
    assert z_true[0] > 5 and z_true[0] > z_true[1] + 3


def test_template_beats_constant_for_non_flat_spectrum():
    frames, mv = _small_field(mag=18.0)
    maps = build_maps(frames, q=2)
    grid = make_grid(frames, half_fov_arcsec=30, vmax_arcsec_day=2)
    xi, eta = mv.pos(grid.t_ref)
    tr = Trajectories(np.array([xi]), np.array([eta]), np.array([mv.vxi]), np.array([mv.veta]), grid.t_ref)
    psi, phi = gather(maps, tr)
    zt, _ = stat_template(psi, phi, spectral_factors(frames, mv.spectrum))
    zc, _ = stat_constant(psi, phi)
    assert zt[0] > zc[0]


def test_branch_and_bound_is_exact():
    frames, mv = _small_field(mag=18.5)
    maps = build_maps(frames, q=2)
    grid = make_grid(frames, half_fov_arcsec=18, vmax_arcsec_day=1.5)
    s = spectral_factors(frames, mv.spectrum)
    thr = 4.0
    found, _ = branch_and_bound(maps, grid, s, thr)
    px, py = grid.positions()
    exhaustive = set()
    for iy, vy in enumerate(grid.veta):
        for ix, vx in enumerate(grid.vxi):
            tr = Trajectories(px, py, np.full(px.size, vx), np.full(px.size, vy), grid.t_ref)
            z, _ = stat_template(*gather(maps, tr), s)
            for j in np.where(z >= thr)[0]:
                exhaustive.add((int(j), iy * grid.vxi.size + ix))
    assert exhaustive, "test needs at least one above-threshold trajectory"
    assert {(a, b) for a, b, _ in found} == exhaustive


# ------------------------------------------------------------------ gravity
def test_magnification_limits():
    assert magnification(1.0) == pytest.approx(3 / np.sqrt(5))
    assert magnification(100.0) == pytest.approx(1.0, abs=1e-6)


def test_gravity_models_are_nested_and_recover_lens():
    t = np.sort(np.random.default_rng(0).uniform(0, 200, 330))
    w = np.exp(np.random.default_rng(1).uniform(np.log(0.75), np.log(5), t.size))
    lc = simulate_event(t, w, t0=90.0, u0=0.4, tE=12.0, blend_frac_2um=0.6, sigma_frac=0.01, seed=2)
    gm = GravityModel(lc, k_spec=3)
    r = gm.compare()
    assert r["H1"].chi2 <= r["H0"].chi2 + 1e-6
    assert r["H1"].params["tE"] == pytest.approx(12.0, rel=0.15)
    assert r["H1"].params["t0"] == pytest.approx(90.0, abs=0.5)
    assert r["verdict"].startswith("Microlensing")


def test_gravity_null_is_not_called_a_lens():
    t = np.sort(np.random.default_rng(5).uniform(0, 200, 330))
    w = np.exp(np.random.default_rng(6).uniform(np.log(0.75), np.log(5), t.size))
    lc = simulate_event(t, w, t0=90.0, u0=0.4, tE=12.0, sigma_frac=0.02, seed=7, lensed=False)
    r = GravityModel(lc, k_spec=3).compare()
    assert r["dchi2_H1_vs_H0"] >= -1e-6
    assert r["verdict"].startswith("Insufficient")


def test_gravity_refuses_inadequate_models():
    """Underestimated errors (as for saturated real stars) must not yield a lens verdict."""
    t = np.sort(np.random.default_rng(8).uniform(0, 200, 200))
    w = np.exp(np.random.default_rng(9).uniform(np.log(0.75), np.log(5), t.size))
    lc = simulate_event(t, w, t0=90.0, u0=1.0, tE=20.0, sigma_frac=0.05, seed=10, lensed=False)
    lc.err = lc.err / 30.0                      # errors 30x too small
    r = GravityModel(lc, k_spec=3).compare()
    assert r["verdict"].startswith("Model inadequate")


def test_mass_posterior_is_broad():
    q = mass_posterior(20.0, 1.0)["M_quantiles_Msun"]
    assert q["95%"] / q["5%"] > 10            # tE alone cannot pin the mass


# ------------------------------------------------------------------- stats
def test_wilson_and_m50():
    lo, hi = wilson(5, 10)
    assert lo < 0.5 < hi
    assert m50([20, 21, 22], [0.9, 0.5, 0.1]) == pytest.approx(21.0)
    assert np.isnan(m50([20, 21], [0.4, 0.1]))


# --------------------------------------------------------------------- API
def test_api_smoke():
    from fastapi.testclient import TestClient
    from api.app import app
    c = TestClient(app)
    t = c.get("/api/targets").json()
    assert any(x["slug"] == "barnard" for x in t)
    r = c.get("/api/resolve", params={"q": "08h55m10.8s -07d14m42s"}).json()
    assert r["ra"] == pytest.approx(133.795, abs=1e-3) and r["dec"] == pytest.approx(-7.245, abs=1e-3)
    assert c.get("/api/resolve", params={"q": "nonsense"}).status_code == 400
    g = c.get("/api/gravity/demo", params={"u0": 0.3, "tE": 15, "seed": 1}).json()
    assert g["provenance"] == "synthetic" and len(g["t"]) == len(g["flux"])
    assert c.post("/api/fields", params={"ra": 400, "dec": 0}).status_code == 400
    assert c.get("/").status_code == 200


def test_psiphi_statistic_matches_direct_pixel_sum():
    """Acceptance 7: the reported statistic equals the underlying computation.
    The psi/phi-map trajectory statistic must equal a direct pixel-level
    matched filter sum_e sum_p w P d / sqrt(sum_e sum_p w P^2) at the same
    (q-grid-rounded) positions."""
    frames, mv = _small_field(mag=18.0)
    maps = build_maps(frames, q=2)
    t_ref = 0.5 * (frames[0].t_mjd + frames[-1].t_mjd)
    xi, eta = mv.pos(t_ref)
    tr = Trajectories(np.array([xi]), np.array([eta]), np.array([mv.vxi]), np.array([mv.veta]), t_ref)
    psi, phi = gather(maps, tr)
    z_maps, _ = stat_constant(psi, phi)
    num = den = 0.0
    for m, f in zip(maps, frames):
        x, y = f.to_pix(*tr.at(f.t_mjd))
        xq, yq = np.round(np.asarray(x) * m.q) / m.q, np.round(np.asarray(y) * m.q) / m.q   # same rounding as maps
        P = f.psf.render(f.shape, float(np.ravel(xq)[0]), float(np.ravel(yq)[0]))
        num += float(np.sum(f.ivar * P * f.data))
        den += float(np.sum(f.ivar * P * P))
    z_direct = num / np.sqrt(den)
    assert z_maps[0] == pytest.approx(z_direct, rel=0.02)
