"""Source spectral models evaluated through the SPHEREx linear-variable-filter
response.

SPHEREx does not observe a source in 102 bands simultaneously.  At any one
exposure a source lands on one pixel of each of the six detectors, and that
pixel's linear-variable filter passes a narrow band (R ~ 35-130) centred on
the pixel's ``CWAVE``.  A moving object therefore yields one flux sample per
(exposure, detector) at a wavelength set by *where it fell on the array*.

All flux densities here are F_nu in micro-Jansky (uJy).  AB magnitude:
m_AB = 23.9 - 2.5 log10(F_nu / uJy).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

H = 6.62607015e-34
C = 2.99792458e8
KB = 1.380649e-23


def ab_to_ujy(m):
    return 10 ** ((23.9 - np.asarray(m, dtype=float)) / 2.5)


def ujy_to_ab(f):
    f = np.asarray(f, dtype=float)
    return 23.9 - 2.5 * np.log10(np.where(f > 0, f, np.nan))


def planck_nu(wave_um, T):
    """B_nu(T) up to a constant (shape only)."""
    lam = np.asarray(wave_um, dtype=float) * 1e-6
    nu = C / lam
    x = np.clip(H * nu / (KB * T), 1e-8, 700)
    return nu ** 3 / np.expm1(x)


def band_average(fn, wave_um, band_um, n: int = 7):
    """Average F_nu over a top-hat band [w - b/2, w + b/2] (Simpson-like)."""
    w = np.asarray(wave_um, dtype=float)
    b = np.asarray(band_um, dtype=float)
    t = np.linspace(-0.5, 0.5, n)
    wts = np.ones(n)
    wts[1:-1:2], wts[2:-1:2] = 4, 2
    wts /= wts.sum()
    out = 0.0
    for ti, wi in zip(t, wts):
        out = out + wi * fn(w + ti * b)
    return out


# --------------------------------------------------------------------------
# Spectral families (each returns relative F_nu, normalised to 1 at 1.0 um
# unless noted, so the linear amplitude carries the brightness).
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class Spectrum:
    name: str
    kind: str
    params: tuple

    def __call__(self, wave_um):
        w = np.asarray(wave_um, dtype=float)
        if self.kind == "flat":
            return np.ones_like(w)
        if self.kind == "powerlaw":  # F_nu ∝ λ^alpha
            (alpha,) = self.params
            return w ** alpha
        if self.kind == "blackbody":
            (T,) = self.params
            return planck_nu(w, T) / planck_nu(1.0, T)
        if self.kind == "reflect":  # solar-type reflectance with slope + optional thermal tail
            slope, t_eq, thermal_frac = self.params
            sun = planck_nu(w, 5778.0) / planck_nu(1.0, 5778.0)
            refl = sun * (1.0 + slope * (w - 1.0))
            if thermal_frac > 0:
                th = planck_nu(w, t_eq) / planck_nu(4.5, t_eq)
                refl = refl + thermal_frac * th * (planck_nu(4.5, 5778.0) / planck_nu(1.0, 5778.0))
            return np.clip(refl, 0, None)
        if self.kind == "browndwarf":  # crude T/Y dwarf proxy: CH4/H2O suppressed except 4-5 um peak
            (teff,) = self.params
            bb = planck_nu(w, teff) / planck_nu(4.6, teff)
            win = 0.15 + 0.85 * np.exp(-0.5 * ((w - 4.4) / 0.35) ** 2) + 0.25 * np.exp(-0.5 * ((w - 1.27) / 0.05) ** 2)
            return bb * win
        raise ValueError(self.kind)


def template_bank() -> list[Spectrum]:
    """Physically motivated bank spanning reflected light, stars and cold sources."""
    bank = [Spectrum("flat (F_nu const)", "flat", ())]
    for s in (-0.1, 0.1, 0.3):
        bank.append(Spectrum(f"reflectance slope {s:+.1f}/um", "reflect", (s, 250.0, 0.0)))
    for te, fr in ((250.0, 0.5), (300.0, 2.0)):
        bank.append(Spectrum(f"reflectance + thermal T={te:.0f}K", "reflect", (0.1, te, fr)))
    for T in (1500.0, 3000.0, 6000.0):
        bank.append(Spectrum(f"blackbody {T:.0f} K", "blackbody", (T,)))
    bank.append(Spectrum("cold brown dwarf proxy 400 K", "browndwarf", (400.0,)))
    return bank


def spectral_basis(wave_um: np.ndarray, k: int) -> np.ndarray:
    """Smooth low-dimensional basis for unknown spectra: log-λ Legendre polynomials
    multiplied by a solar envelope (so k=1 is a solar-reflectance spectrum).

    Returns an (n, k) matrix.  The amplitude vector is fitted linearly; the
    detection statistic then follows a chi^2_k null distribution.
    """
    w = np.asarray(wave_um, dtype=float)
    x = (np.log(w) - np.log(0.75)) / (np.log(5.0) - np.log(0.75)) * 2 - 1
    env = planck_nu(w, 5778.0) / planck_nu(1.0, 5778.0)
    cols = [np.polynomial.legendre.legval(x, np.eye(k)[i]) * env for i in range(k)]
    return np.stack(cols, axis=-1)
