r"""Gravity Mode: point-source point-lens (PSPL) microlensing with spectral
blending, evaluated directly against multi-wavelength flux measurements.

Measurement model (one flux sample i per exposure/detector at wavelength
lambda_i and time t_i, Gaussian errors sigma_i):

    F_i = A(t_i; t0, u0, tE) * S(lambda_i) + Bl(lambda_i) + noise

with achromatic magnification (Paczynski 1986)

    A(u) = (u^2 + 2) / (u sqrt(u^2 + 4)),   u(t) = sqrt(u0^2 + ((t - t0)/tE)^2)

S and Bl are the unknown source and blend SEDs.  We expand both in a smooth
K-dimensional spectral basis, so for fixed nonlinear parameters (t0, u0, tE)
the model is *linear* in 2K amplitudes, which we profile out analytically by
weighted least squares.  The nonlinear search is then only 3-D.

Hypotheses (all with analytic amplitude profiling):

* H0  constant source:        F = S(lambda)
* H1  PSPL + blend:           F = A(t) S(lambda) + Bl(lambda)
* H2  intrinsic variability:  F = S(lambda) [1 + sum_j c_j g_j(t)]  -- smooth
      achromatic variability on a basis of Gaussian bumps spanning the data
      (the same *number* of extra linear parameters as H1 has nonlinear ones,
      chosen so that a single symmetric bump can mimic a microlensing peak).
* H2c chromatic variability:  F = S(lambda) + g(t) C(lambda)

Blending makes the *observed* fractional brightening wavelength dependent even
though lensing is achromatic -- this is what a joint spectro-temporal fit can
exploit and what per-band light-curve fitting discards.

Physical interpretation:  tE = theta_E / mu_rel, theta_E = sqrt(kappa M pi_rel)
with kappa = 8.144 mas / M_sun.  A light curve alone constrains only tE, so mass
is degenerate with distance and proper motion; ``mass_posterior`` makes that
explicit with a Galactic-model prior rather than quoting a single mass.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from ..science.spectra import spectral_basis


def magnification(u):
    u = np.maximum(np.asarray(u, float), 1e-12)
    return (u * u + 2.0) / (u * np.sqrt(u * u + 4.0))


def u_of_t(t, t0, u0, tE):
    return np.sqrt(u0 * u0 + ((np.asarray(t, float) - t0) / tE) ** 2)


@dataclass
class LightCurve:
    t: np.ndarray          # MJD
    flux: np.ndarray       # uJy
    err: np.ndarray        # uJy
    wave: np.ndarray       # um
    label: str = ""
    provenance: str = "synthetic"
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        for k in ("t", "flux", "err", "wave"):
            setattr(self, k, np.asarray(getattr(self, k), float))


def _wls(X: np.ndarray, y: np.ndarray, w: np.ndarray):
    """Weighted least squares; returns (coef, chi2, cov)."""
    Xw = X * w[:, None]
    G = X.T @ Xw
    G += np.eye(G.shape[0]) * 1e-12 * np.trace(G)
    coef = np.linalg.solve(G, Xw.T @ y)
    r = y - X @ coef
    return coef, float(np.sum(w * r * r)), np.linalg.inv(G)


@dataclass
class FitResult:
    model: str
    chi2: float
    n_par: int
    params: dict
    coef: np.ndarray
    model_flux: np.ndarray

    def bic(self, n):
        return self.chi2 + self.n_par * np.log(n)


class GravityModel:
    def __init__(self, lc: LightCurve, k_spec: int = 2):
        self.lc = lc
        self.k = k_spec
        self.Bs = spectral_basis(lc.wave, k_spec)
        # the static SED of (source + blend) gets the same spectral freedom in
        # every hypothesis (2K dims) so that H1 cannot win merely by having a
        # more flexible SED than H0/H2 (found in the first null test).
        self.B2 = spectral_basis(lc.wave, 2 * k_spec)
        self.w = 1.0 / lc.err ** 2
        self.n = lc.t.size

    # ------------------------------------------------------------ hypotheses
    def fit_constant(self) -> FitResult:
        c, x2, _ = _wls(self.B2, self.lc.flux, self.w)
        return FitResult("H0 constant", x2, 2 * self.k, {}, c, self.B2 @ c)

    def design_pspl(self, t0, u0, tE):
        A = magnification(u_of_t(self.lc.t, t0, u0, tE))
        # blend/static term uses the 2K basis so H0 is nested in H1 (A -> 1)
        return np.hstack([self.Bs * A[:, None], self.B2])

    def chi2_pspl(self, t0, u0, tE) -> tuple[float, np.ndarray]:
        X = self.design_pspl(t0, u0, tE)
        c, x2, _ = _wls(X, self.lc.flux, self.w)
        return x2, c

    def fit_pspl(self, t0_grid=None, tE_grid=None, u0_grid=None) -> FitResult:
        t = self.lc.t
        t0_grid = t0_grid if t0_grid is not None else np.linspace(t.min(), t.max(), 60)
        tE_grid = tE_grid if tE_grid is not None else np.geomspace(0.5, 200, 30)
        u0_grid = u0_grid if u0_grid is not None else np.array([0.01, 0.05, 0.1, 0.2, 0.4, 0.7, 1.0, 1.5])
        best = (np.inf, None)
        for t0 in t0_grid:
            for tE in tE_grid:
                for u0 in u0_grid:
                    x2, _ = self.chi2_pspl(t0, u0, tE)
                    if x2 < best[0]:
                        best = (x2, (t0, np.log(u0), np.log(tE)))
        # local refinement in (t0, ln u0, ln tE)
        f = lambda p: self.chi2_pspl(p[0], np.exp(p[1]), np.exp(p[2]))[0]
        res = minimize(f, np.array(best[1]), method="Nelder-Mead",
                       options={"xatol": 1e-4, "fatol": 1e-3, "maxiter": 3000})
        t0, u0, tE = res.x[0], float(np.exp(res.x[1])), float(np.exp(res.x[2]))
        x2, c = self.chi2_pspl(t0, u0, tE)
        X = self.design_pspl(t0, u0, tE)
        k = self.k
        wref = 2.0
        src = float((spectral_basis(np.array([wref]), k) @ c[:k])[0])
        bl = float((spectral_basis(np.array([wref]), 2 * k) @ c[k:])[0])
        fs = src / (src + bl) if (src + bl) != 0 else np.nan
        return FitResult("H1 point lens + blend", x2, 3 * k + 3,
                         {"t0": t0, "u0": u0, "tE": tE, "Amax": float(magnification(u0)),
                          "source_frac_2um": fs}, c, X @ c)

    def fit_variability(self, n_bumps: int | None = None, width_days: float | None = None) -> FitResult:
        """Achromatic smooth variability: S(lambda) * (1 + sum c_j g_j(t)).

        This is nonlinear (product), so we use the standard linearisation that
        is exact for an achromatic modulation: F = sum_j S(lambda) g_j(t) with
        g_0 = 1, i.e. the design matrix is the outer product of the spectral
        basis and the temporal basis.  The best single bump (centre, width) is
        found on a grid, giving H2 the same nonlinear freedom as the lens.
        """
        t = self.lc.t
        best = None
        centres = np.linspace(t.min(), t.max(), 60)
        widths = np.geomspace(0.5, 200, 30) if width_days is None else [width_days]
        for c0 in centres:
            for wd in widths:
                g = np.exp(-0.5 * ((t - c0) / wd) ** 2)
                X = np.hstack([self.B2, self.Bs * g[:, None]])
                coef, x2, _ = _wls(X, self.lc.flux, self.w)
                if best is None or x2 < best[0]:
                    best = (x2, c0, wd, coef, X)
        x2, c0, wd, coef, X = best
        return FitResult("H2 achromatic variability (Gaussian bump)", x2, 3 * self.k + 2,
                         {"centre": c0, "width": wd}, coef, X @ coef)

    def fit_periodic(self, periods=None) -> FitResult:
        t = self.lc.t
        periods = periods if periods is not None else np.geomspace(1.0, 600.0, 200)
        best = None
        for P in periods:
            ph = 2 * np.pi * t / P
            X = np.hstack([self.B2, self.Bs * np.sin(ph)[:, None], self.Bs * np.cos(ph)[:, None]])
            coef, x2, _ = _wls(X, self.lc.flux, self.w)
            if best is None or x2 < best[0]:
                best = (x2, P, coef, X)
        x2, P, coef, X = best
        return FitResult("H2p periodic variability", x2, 4 * self.k + 1, {"period": P}, coef, X @ coef)

    def compare(self) -> dict:
        h0 = self.fit_constant()
        h1 = self.fit_pspl()
        h2 = self.fit_variability()
        n = self.n
        res = {"n": n, "H0": h0, "H1": h1, "H2": h2,
               "dchi2_H1_vs_H0": h0.chi2 - h1.chi2, "dchi2_H1_vs_H2": h2.chi2 - h1.chi2,
               "dBIC_H1_vs_H0": h0.bic(n) - h1.bic(n), "dBIC_H1_vs_H2": h2.bic(n) - h1.bic(n)}
        # --- guards (added after a real saturated star produced a nonsensical
        # "lens" fit): model adequacy and physical validity come before BIC.
        span = float(np.ptp(self.lc.t)) if n > 1 else 0.0
        chi2r = {k: res[k].chi2 / max(1, n - res[k].n_par) for k in ("H0", "H1", "H2")}
        p1 = h1.params
        physical = (0.0 < (p1.get("source_frac_2um") or -1) <= 1.05) and (p1["tE"] <= 3.0 * max(span, 1e-3))             and (self.lc.t.min() - p1["tE"] <= p1["t0"] <= self.lc.t.max() + p1["tE"])
        res["chi2r"] = chi2r
        res["lens_physical"] = bool(physical)
        best_r = min(chi2r["H1"], chi2r["H2"])
        if best_r > 3.0:
            verdict = (f"Model inadequate (best reduced chi2 = {best_r:.0f}); systematics, saturation or an "
                       "unmodelled spectrum dominate. No inference made")
        elif res["dBIC_H1_vs_H0"] < 10:
            verdict = "Insufficient evidence for variability"
        elif res["dBIC_H1_vs_H2"] > 10 and physical:
            verdict = "Microlensing-compatible signal (point lens favoured over smooth variability)"
        elif res["dBIC_H1_vs_H2"] > 10:
            verdict = "Variable; best lens fit is unphysical (t_E, t0 or source fraction out of range)"
        elif res["dBIC_H1_vs_H2"] < -10:
            verdict = "Alternative model favoured (intrinsic variability)"
        else:
            verdict = "Variable; point lens and intrinsic variability not distinguished"
        res["verdict"] = verdict
        return res


# ------------------------------------------------------------------ simulation
def simulate_event(t, wave, t0, u0, tE, src_mag_2um=16.0, blend_frac_2um=0.5, src_T=4500.0, blend_T=3200.0,
                   sigma_frac=0.01, seed=0, lensed: bool = True) -> LightCurve:
    """Synthetic multi-wavelength light curve (provenance synthetic)."""
    from ..science.spectra import Spectrum, ab_to_ujy, planck_nu

    rng = np.random.default_rng(seed)
    t = np.asarray(t, float)
    wave = np.asarray(wave, float)
    s_sed = planck_nu(wave, src_T) / planck_nu(2.0, src_T)
    b_sed = planck_nu(wave, blend_T) / planck_nu(2.0, blend_T)
    f2 = ab_to_ujy(src_mag_2um)
    fs = f2 * (1 - blend_frac_2um) * s_sed
    fb = f2 * blend_frac_2um * b_sed
    A = magnification(u_of_t(t, t0, u0, tE)) if lensed else np.ones_like(t)
    true = A * fs + fb
    err = sigma_frac * (fs + fb) + 0.0
    flux = true + rng.normal(0, err)
    return LightCurve(t, flux, err, wave, label="synthetic PSPL" if lensed else "synthetic constant",
                      provenance="synthetic",
                      meta={"t0": t0, "u0": u0, "tE": tE, "blend_frac_2um": blend_frac_2um, "lensed": lensed,
                            "true_flux": true})


# ------------------------------------------------------------- mass degeneracy
KAPPA = 8.144  # mas / M_sun


def mass_posterior(tE_days: float, tE_err: float, n: int = 200_000, seed: int = 0, ds_kpc: float = 8.0) -> dict:
    """Monte-Carlo lens-mass distribution given only tE (illustrates the
    mass-distance-velocity degeneracy).  Priors (explicit, simple):
    lens distance uniform-in-volume along a line of sight to a source at
    ``ds_kpc``; relative proper motion ~ Rayleigh with sigma=4 mas/yr; mass is
    whatever makes tE match.  Returned quantiles are prior-dominated by design.
    """
    rng = np.random.default_rng(seed)
    x = rng.uniform(0.0, 1.0, n) ** (1 / 3)              # D_L / D_S, volume weighting
    x = np.clip(x, 0.02, 0.98)
    dl = x * ds_kpc
    pirel = (1.0 / dl - 1.0 / ds_kpc)                     # mas (kpc^-1 -> mas)
    mu = rng.rayleigh(4.0, n)                             # mas/yr
    tE = rng.normal(tE_days, tE_err, n)
    tE = tE[tE > 0][: n]
    m = min(tE.size, mu.size)
    thetaE = mu[:m] * tE[:m] / 365.25                     # mas
    M = thetaE ** 2 / (KAPPA * pirel[:m])
    q = np.quantile(M, [0.05, 0.16, 0.5, 0.84, 0.95])
    return {"M_quantiles_Msun": dict(zip(["5%", "16%", "50%", "84%", "95%"], map(float, q))),
            "assumptions": "D_S=%.0f kpc, volume-weighted D_L, Rayleigh mu_rel (sigma=4 mas/yr); "
                           "no parallax or finite-source constraint" % ds_kpc}
