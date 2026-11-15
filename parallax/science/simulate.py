"""Synthetic SPHEREx-like observations (provenance = "synthetic").

The simulator reproduces the properties that matter for moving-object
inference:

* six detectors exposed simultaneously, each at a *different* wavelength,
  set by where the source falls on that detector's linear-variable filter;
* irregular exposure times and arbitrary field rotation per exposure;
* an undersampled pixel-integrated PSF (FWHM ~ 1.2 px);
* per-pixel Gaussian noise with detector-dependent depth, masked pixels;
* a static star field with stellar spectra, plus optional moving sources.

Real-background injections (``inject_into_frames``) reuse real SPHEREx frames
and add only the synthetic source, so noise, artefacts and crowding are real.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..inference.frames import Frame
from .psf import GaussianPSF, add_source
from .spectra import Spectrum, ab_to_ujy, band_average

DET_RANGE = {1: (0.75, 1.09), 2: (1.10, 1.62), 3: (1.63, 2.41), 4: (2.42, 3.82), 5: (3.83, 4.41), 6: (4.42, 5.00)}
DET_R = {1: 39, 2: 41, 3: 41, 4: 35, 5: 112, 6: 128}
# Per-pixel 1-sigma noise in uJy for a single exposure.  Defaults are of the
# order measured on DR1 cutouts (D1 ~ 16 uJy at ecliptic latitude -25 deg); the
# benchmark scripts overwrite them with values measured from real cutouts.
DEFAULT_SIGMA_UJY = {1: 16.0, 2: 14.0, 3: 12.0, 4: 11.0, 5: 22.0, 6: 26.0}
PIXSCALE = 6.15


@dataclass
class MovingSource:
    xi0: float           # arcsec offset at t_ref
    eta0: float
    vxi: float           # arcsec / day
    veta: float
    mag_ab_1um: float    # AB magnitude of F_nu at 1.0 um
    spectrum: Spectrum
    t_ref: float
    axi: float = 0.0     # arcsec / day^2
    aeta: float = 0.0

    def pos(self, t):
        dt = np.asarray(t, float) - self.t_ref
        return (self.xi0 + self.vxi * dt + 0.5 * self.axi * dt ** 2,
                self.eta0 + self.veta * dt + 0.5 * self.aeta * dt ** 2)

    def flux(self, wave_um, band_um):
        return ab_to_ujy(self.mag_ab_1um) * band_average(self.spectrum, wave_um, band_um)


@dataclass
class SimConfig:
    n_exposures: int = 12
    span_days: float = 1.5
    detectors: tuple = (1, 2, 3, 4, 5, 6)
    size: int = 48
    sigma_ujy: dict = field(default_factory=lambda: dict(DEFAULT_SIGMA_UJY))
    mask_frac: float = 0.02
    star_density_per_deg2: float = 6000.0
    star_mag_range: tuple = (12.0, 21.0)
    psf_fwhm_px: float = 1.22
    seed: int = 0
    irregular: bool = True
    # Real SPHEREx visit structure (measured from DR1 headers at Barnard's Star):
    # bursts of 2-4 exposures 0.00149 d apart, each step moving the source
    # ~115 px along the LVF (one spectral channel); only one detector pair
    # (1,4), (2,5) or (3,6) sees the source per burst; ~25 bursts in ~20 days.
    spherex_bursts: bool = False
    n_bursts: int = 25
    burst_len: tuple = (2, 4)


PAIRS = [(1, 4), (2, 5), (3, 6)]
LVF_STEP = 115.0 / 2040.0
BURST_DT = 0.0014923


def make_schedule_spherex(cfg: "SimConfig", rng: np.random.Generator) -> list[dict]:
    t_b = np.sort(rng.uniform(0, cfg.span_days, cfg.n_bursts))
    sched = []
    iexp = 0
    for tb in t_b:
        pair = PAIRS[rng.integers(0, 3)]
        n = int(rng.integers(cfg.burst_len[0], cfg.burst_len[1] + 1))
        frac0 = rng.uniform(0.02, 0.98 - LVF_STEP * (n - 1))
        pa = rng.uniform(0, 2 * np.pi)
        for k in range(n):
            frac = frac0 + LVF_STEP * k
            for d in pair:
                if d not in cfg.detectors:
                    continue
                lo, hi = DET_RANGE[d]
                w = lo * (hi / lo) ** frac
                sched.append({"t": 60800.0 + tb + k * BURST_DT, "det": d, "wave": w, "band": w / DET_R[d],
                              "pa": pa, "exp": iexp})
            iexp += 1
    return sched


def make_schedule(cfg: SimConfig, rng: np.random.Generator) -> list[dict]:
    """Exposure schedule: times, per-detector wavelengths, rotations."""
    if cfg.spherex_bursts:
        return make_schedule_spherex(cfg, rng)
    if cfg.irregular:
        t = np.sort(rng.uniform(0, cfg.span_days, cfg.n_exposures))
    else:
        t = np.linspace(0, cfg.span_days, cfg.n_exposures)
    sched = []
    for i, ti in enumerate(t):
        pa = rng.uniform(0, 2 * np.pi)
        frac = rng.uniform(0.02, 0.98)  # position along the LVF (same for all detectors: dichroic)
        for d in cfg.detectors:
            lo, hi = DET_RANGE[d]
            w = lo * (hi / lo) ** frac
            sched.append({"t": 60800.0 + ti, "det": d, "wave": w, "band": w / DET_R[d], "pa": pa, "exp": i})
    return sched


def _affine(pa: float, center: float) -> tuple[np.ndarray, np.ndarray]:
    c, s = np.cos(pa), np.sin(pa)
    # east-left parity: pixel x increases toward -xi (like a sky image), rotated by pa
    A = np.array([[-c, s], [s, c]]) / PIXSCALE
    return A, np.array([center, center])


def simulate(cfg: SimConfig, movers: list[MovingSource] | None = None, include_stars: bool = True,
             psf=None) -> tuple[list[Frame], dict]:
    rng = np.random.default_rng(cfg.seed)
    psf = psf or GaussianPSF(cfg.psf_fwhm_px)
    sched = make_schedule(cfg, rng)
    half_fov = cfg.size * PIXSCALE / 2
    # static star field over a slightly larger area (stars just outside can bleed in)
    area_deg2 = (2 * half_fov * 1.2 / 3600.0) ** 2
    nstar = rng.poisson(cfg.star_density_per_deg2 * area_deg2) if include_stars else 0
    # N(<m) ∝ 10^(0.3 m): draw magnitudes from that law within range
    lo, hi = cfg.star_mag_range
    u = rng.uniform(size=nstar)
    mags = np.log10(10 ** (0.3 * lo) + u * (10 ** (0.3 * hi) - 10 ** (0.3 * lo))) / 0.3
    sxi = rng.uniform(-1.2 * half_fov, 1.2 * half_fov, nstar)
    seta = rng.uniform(-1.2 * half_fov, 1.2 * half_fov, nstar)
    stemps = rng.uniform(3000, 7000, nstar)
    stars = [Spectrum(f"bb{T:.0f}", "blackbody", (T,)) for T in stemps]
    frames: list[Frame] = []
    movers = movers or []
    t_ref = float(np.mean([s["t"] for s in sched]))
    for s in sched:
        A, b = _affine(s["pa"], (cfg.size - 1) / 2)
        img = np.zeros((cfg.size, cfg.size))
        for k in range(nstar):
            x, y = A @ np.array([sxi[k], seta[k]]) + b
            if -4 < x < cfg.size + 4 and -4 < y < cfg.size + 4:
                f = ab_to_ujy(mags[k]) * float(band_average(stars[k], s["wave"], s["band"]))
                add_source(img, psf, x, y, f, half=5)
        for m in movers:
            mx, my = m.pos(s["t"])
            x, y = A @ np.array([mx, my]) + b
            add_source(img, psf, x, y, float(m.flux(s["wave"], s["band"])), half=5)
        sig = cfg.sigma_ujy[s["det"]]
        noise = rng.normal(0, sig, img.shape)
        ivar = np.full(img.shape, 1.0 / sig ** 2)
        bad = rng.uniform(size=img.shape) < cfg.mask_frac
        ivar[bad] = 0.0
        data = np.where(bad, 0.0, img + noise)
        frames.append(Frame(t_mjd=s["t"], data=data, ivar=ivar, wave_um=s["wave"], band_um=s["band"],
                            detector=s["det"], psf=psf, A=A, b=b, label=f"sim e{s['exp']:02d} D{s['det']}",
                            provenance="synthetic", meta={"exp": s["exp"], "pa": s["pa"]}))
    truth = {"t_ref": t_ref, "stars": {"xi": sxi, "eta": seta, "mag": mags, "T": stemps}, "movers": movers,
             "noiseless": None}
    return frames, truth


def inject_into_frames(frames: list[Frame], movers: list[MovingSource], seed: int = 0) -> list[Frame]:
    """Return copies of (real) frames with synthetic moving sources added.

    Provenance becomes "real+injected"; masked pixels stay masked.
    """
    out = []
    for f in frames:
        img = np.zeros_like(f.data)
        for m in movers:
            mx, my = m.pos(f.t_mjd)
            x, y = f.to_pix(mx, my)
            add_source(img, f.psf, float(x), float(y), float(m.flux(f.wave_um, f.band_um)), half=5)
        nf = Frame(**{**f.__dict__})
        nf.data = np.where(f.ivar > 0, f.data + img, 0.0)
        nf.provenance = "real+injected" if f.provenance.startswith("real") else f.provenance
        out.append(nf)
    return out
