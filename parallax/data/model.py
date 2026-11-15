"""Typed observation data model.

Units are carried explicitly in field names (``_mjysr``, ``_ujy``, ``_um``) and
in the ``units`` dict so that physically meaningless combinations are hard to
write by accident.  All arrays use numpy (row, column) = (y, x) order; pixel
coordinates are 0-based.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from astropy.wcs import WCS

ARCSEC2_TO_SR = (np.pi / 180.0 / 3600.0) ** 2
MJY_TO_UJY = 1e12

Provenance = Literal["real", "synthetic", "real+injected"]

# SPHEREx FLAGS bit planes (from the MSKNnnnn/MSKMnnnn header cards of the
# DR1 Level-2 FLAGS extension; verified against a real file header).
FLAG_BITS = {
    "TRANSIENT": 0, "OVERFLOW": 1, "SUR_ERROR": 2, "NONFUNC": 6, "DICHROIC": 7,
    "MISSING_DATA": 9, "HOT": 10, "COLD": 11, "FULLSAMPLE": 12, "PHANMISS": 14,
    "NONLINEAR": 15, "PERSIST": 17, "OUTLIER": 19, "CROSSTALK": 20, "SOURCE": 21,
    "GHOST": 22, "GHOST_FPA": 23, "GHOST_EXT": 24, "STREAK": 25, "BLOOM": 26,
    "SNOWBALL": 27, "HALO": 28,
}
# Pixels with these bits are excluded from quantitative inference.  SOURCE
# (bit 21) merely marks a known source and FULLSAMPLE (12) is informational, so
# they are *not* rejected -- rejecting SOURCE would mask every star.
REJECT_BITS = ["TRANSIENT", "OVERFLOW", "SUR_ERROR", "NONFUNC", "MISSING_DATA", "HOT",
               "COLD", "NONLINEAR", "PERSIST", "OUTLIER"]
REJECT_MASK = int(sum(1 << FLAG_BITS[b] for b in REJECT_BITS))


@dataclass
class ObservationMeta:
    """Metadata for one SPHEREx detector exposure (one Level-2 MEF)."""

    obs_id: str
    detector: int
    release: str
    url: str
    mjd_beg: float
    mjd_avg: float
    mjd_end: float
    exptime_s: float
    time_scale: str = "utc"
    sc_pos_km: tuple[float, float, float] | None = None  # geocentric, mid-exposure
    psf_fwhm_arcsec: float | None = None
    omega_median_arcsec2: float | None = None
    source: str = "NASA/IPAC IRSA SPHEREx (AWS Open Data mirror)"

    @property
    def band_name(self) -> str:
        return {1: "D1 0.75-1.09um", 2: "D2 1.10-1.62um", 3: "D3 1.63-2.41um",
                4: "D4 2.42-3.82um", 5: "D5 3.83-4.41um", 6: "D6 4.42-5.00um"}[self.detector]


@dataclass
class Cutout:
    """A calibrated cutout of one exposure with everything needed for inference."""

    meta: ObservationMeta
    image_mjysr: np.ndarray          # IMAGE, surface brightness (zodi NOT subtracted)
    variance_mjysr2: np.ndarray      # VARIANCE
    flags: np.ndarray                # FLAGS bitmap (int32)
    zodi_mjysr: np.ndarray           # ZODI model (to subtract exactly once)
    wave_um: np.ndarray              # per-pixel central wavelength
    band_um: np.ndarray              # per-pixel bandwidth
    wcs: WCS                         # celestial WCS of the cutout (SIP included)
    x0: int                          # cutout origin in the full detector frame
    y0: int
    psf_oversampled: np.ndarray      # ePSF, oversampled
    psf_oversamp: int
    psf_zone_center: tuple[float, float]
    omega_arcsec2: np.ndarray | float  # pixel solid angle
    wave_source: str = "WCS-WAVE (bilinear, ~1 nm)"
    provenance: Provenance = "real"
    notes: list[str] = field(default_factory=list)

    # ---------------------------------------------------------- unit handling
    @property
    def mjysr_to_ujy(self) -> np.ndarray | float:
        """Multiply surface brightness [MJy/sr] by this to get flux per pixel [uJy]."""
        return np.asarray(self.omega_arcsec2) * ARCSEC2_TO_SR * MJY_TO_UJY

    def good(self) -> np.ndarray:
        g = (self.flags & REJECT_MASK) == 0
        g &= np.isfinite(self.image_mjysr) & np.isfinite(self.variance_mjysr2)
        g &= self.variance_mjysr2 > 0
        return g

    def flux_ujy(self, subtract_zodi: bool = True) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (flux, variance, good) per pixel in uJy and uJy^2.

        The ZODI extension is subtracted here and nowhere else (it is *not*
        pre-subtracted in IMAGE), and variance is transformed with the square
        of the same conversion factor.
        """
        k = self.mjysr_to_ujy
        img = self.image_mjysr - (self.zodi_mjysr if subtract_zodi else 0.0)
        return img * k, self.variance_mjysr2 * np.square(k), self.good()

    @property
    def shape(self) -> tuple[int, int]:
        return self.image_mjysr.shape  # type: ignore[return-value]

    def pixel_of(self, ra: float, dec: float) -> tuple[float, float]:
        x, y = self.wcs.all_world2pix([[ra, dec]], 0)[0]
        return float(x), float(y)

    def summary(self) -> dict:
        ny, nx = self.shape
        cy, cx = ny // 2, nx // 2
        return {
            "obs_id": self.meta.obs_id, "detector": self.meta.detector,
            "release": self.meta.release, "mjd": self.meta.mjd_avg,
            "exptime_s": self.meta.exptime_s, "wave_um": float(self.wave_um[cy, cx]),
            "band_um": float(self.band_um[cy, cx]), "provenance": self.provenance,
            "source": self.meta.source, "url": self.meta.url, "shape": [ny, nx],
            "x0": self.x0, "y0": self.y0, "wave_source": self.wave_source,
        }
