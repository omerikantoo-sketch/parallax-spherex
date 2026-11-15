"""Inference-ready frames.

A ``Frame`` is one exposure on one detector, expressed in physical units and
linked to a common reference tangent plane.  The *data stay on their native
pixel grid*: hypotheses are projected into each frame (forward modelling)
rather than resampling the data, which would correlate the noise.

Reference plane: gnomonic offsets (xi, eta) in arcsec about (ra0, dec0), with
xi toward +RA (east) and eta toward +Dec (north).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from astropy.stats import sigma_clipped_stats

from ..data.model import Cutout
from ..science.psf import EPSF, GaussianPSF


def tangent_offsets(ra, dec, ra0, dec0):
    """ICRS (deg) -> gnomonic offsets (arcsec) about (ra0, dec0)."""
    ra, dec, ra0, dec0 = (np.radians(np.asarray(v, dtype=float)) for v in (ra, dec, ra0, dec0))
    cosc = np.sin(dec0) * np.sin(dec) + np.cos(dec0) * np.cos(dec) * np.cos(ra - ra0)
    xi = np.cos(dec) * np.sin(ra - ra0) / cosc
    eta = (np.cos(dec0) * np.sin(dec) - np.sin(dec0) * np.cos(dec) * np.cos(ra - ra0)) / cosc
    return np.degrees(xi) * 3600.0, np.degrees(eta) * 3600.0


def tangent_to_sky(xi, eta, ra0, dec0):
    xi, eta = np.radians(np.asarray(xi, float) / 3600.0), np.radians(np.asarray(eta, float) / 3600.0)
    ra0, dec0 = np.radians(ra0), np.radians(dec0)
    rho = np.hypot(xi, eta)
    c = np.arctan(rho)
    with np.errstate(invalid="ignore", divide="ignore"):
        dec = np.arcsin(np.cos(c) * np.sin(dec0) + np.where(rho > 0, eta * np.sin(c) * np.cos(dec0) / rho, 0))
        ra = ra0 + np.arctan2(xi * np.sin(c), rho * np.cos(dec0) * np.cos(c) - eta * np.sin(dec0) * np.sin(c))
    return np.degrees(ra) % 360.0, np.degrees(dec)


@dataclass
class Frame:
    t_mjd: float
    data: np.ndarray          # uJy / pixel, background subtracted
    ivar: np.ndarray          # 1/uJy^2, 0 where masked
    wave_um: float            # wavelength at the frame (cutout) centre
    band_um: float
    detector: int
    psf: object               # EPSF or GaussianPSF
    A: np.ndarray             # 2x2: d(pixel)/d(arcsec offset)
    b: np.ndarray             # pixel position of (xi, eta) = (0, 0)
    label: str = ""
    provenance: str = "real"
    wave_map: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    def to_pix(self, xi, eta):
        xi = np.asarray(xi, float)
        eta = np.asarray(eta, float)
        x = self.A[0, 0] * xi + self.A[0, 1] * eta + self.b[0]
        y = self.A[1, 0] * xi + self.A[1, 1] * eta + self.b[1]
        return x, y

    @property
    def shape(self):
        return self.data.shape

    @property
    def pixscale(self) -> float:
        """arcsec per pixel (geometric mean)."""
        return float(1.0 / np.sqrt(abs(np.linalg.det(self.A))))


def frame_from_cutout(c: Cutout, ra0: float, dec0: float, use_epsf: bool = True,
                      subtract_background: bool = True) -> Frame:
    """Convert a calibrated cutout into a Frame (uJy, inverse variance, affine map).

    * ZODI is subtracted exactly once (inside ``Cutout.flux_ujy``).
    * A residual constant background is removed with a sigma-clipped median of
      unflagged pixels (the SOURCE flag is ignored for this estimate only by
      excluding SOURCE-flagged pixels from the statistic, not from the data).
    * The sky->pixel map is the local affine Jacobian of the full TAN-SIP WCS,
      evaluated by finite differences at the reference point; its error over a
      64-px cutout is checked in tests (< 0.02 px).
    """
    flux, var, good = c.flux_ujy(subtract_zodi=True)
    data = np.where(good, flux, 0.0).astype(float)
    ivar = np.where(good, 1.0 / np.where(var > 0, var, np.inf), 0.0)
    if subtract_background:
        src = (c.flags & (1 << 21)) != 0
        sel = good & ~src
        if sel.sum() < 50:
            sel = good
        _, med, _ = sigma_clipped_stats(flux[sel], sigma=3.0, maxiters=5)
        data = np.where(good, data - med, 0.0)
    # affine map from tangent offsets to cutout pixels
    eps = 3.0  # arcsec
    from .frames import tangent_to_sky as t2s  # local import keeps module import order simple
    pts_xi = np.array([0.0, eps, 0.0])
    pts_eta = np.array([0.0, 0.0, eps])
    ra, dec = t2s(pts_xi, pts_eta, ra0, dec0)
    pix = c.wcs.all_world2pix(np.stack([ra, dec], -1), 0)
    b = pix[0]
    A = np.stack([(pix[1] - pix[0]) / eps, (pix[2] - pix[0]) / eps], axis=-1)
    ny, nx = c.shape
    psf = EPSF(c.psf_oversampled, c.psf_oversamp) if use_epsf else GaussianPSF()
    return Frame(t_mjd=c.meta.mjd_avg, data=data, ivar=ivar, wave_um=float(c.wave_um[ny // 2, nx // 2]),
                 band_um=float(c.band_um[ny // 2, nx // 2]), detector=c.meta.detector, psf=psf, A=A, b=b,
                 label=f"{c.meta.obs_id} D{c.meta.detector}", provenance=c.provenance, wave_map=c.wave_um,
                 meta={"obs_id": c.meta.obs_id, "url": c.meta.url, "exptime_s": c.meta.exptime_s})
