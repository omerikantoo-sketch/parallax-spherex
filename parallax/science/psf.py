"""Point-spread-function models.

SPHEREx is strongly undersampled (FWHM ~ 7.5" on 6.15" pixels), so a PSF must
be evaluated as a *pixel-integrated* response at a sub-pixel source position.
Two models are provided:

* ``EPSF`` -- the mission's effective PSF (Anderson & King 2000 convention),
  oversampled by ``ovs`` (5 for DR1).  An ePSF already includes the pixel
  response, so the value recorded in a native pixel whose centre is offset by
  (dx, dy) from the source is ``ovs**2 * epsf(c + ovs*dx, c + ovs*dy)`` when
  the oversampled ePSF sums to 1 (verified on DR1 files: native-phase sums
  are 0.040 = 1/25 for every phase).
* ``GaussianPSF`` -- an analytically pixel-integrated circular Gaussian (erf
  form), used for synthetic data and as an exact reference implementation.

Both expose ``render(shape, x, y)`` returning a unit-flux image (sum ~ 1 when
the PSF is fully inside the stamp), with x, y 0-based pixel coordinates of the
source in the stamp frame.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.special import erf


@dataclass
class GaussianPSF:
    fwhm_px: float = 1.22  # 7.5" / 6.15"

    @property
    def sigma(self) -> float:
        return self.fwhm_px / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    def _axis(self, n: int, c: float) -> np.ndarray:
        s = self.sigma * np.sqrt(2.0)
        edges = np.arange(n + 1) - 0.5 - c
        return 0.5 * np.diff(erf(edges / s))

    def render(self, shape: tuple[int, int], x: float, y: float) -> np.ndarray:
        ny, nx = shape
        return np.outer(self._axis(ny, y), self._axis(nx, x))

    def render_local(self, x: float, y: float, half: int = 3):
        """Stamp of size (2*half+1)^2 around the nearest pixel -> (img, ix0, iy0)."""
        ix, iy = int(np.floor(x + 0.5)), int(np.floor(y + 0.5))
        n = 2 * half + 1
        return self.render((n, n), x - (ix - half), y - (iy - half)), ix - half, iy - half


@dataclass
class EPSF:
    data: np.ndarray  # oversampled ePSF, normalised to sum 1
    ovs: int = 5

    def __post_init__(self):
        d = np.asarray(self.data, dtype=float)
        d = np.where(np.isfinite(d), d, 0.0)
        s = d.sum()
        self.data = d / s if s > 0 else d
        self.c = (self.data.shape[0] - 1) / 2.0

    @property
    def half_extent_px(self) -> float:
        return self.c / self.ovs

    def render(self, shape: tuple[int, int], x: float, y: float) -> np.ndarray:
        ny, nx = shape
        gy, gx = np.mgrid[0:ny, 0:nx].astype(float)
        u = self.c + self.ovs * (gx - x)
        v = self.c + self.ovs * (gy - y)
        val = map_coordinates(self.data, [v.ravel(), u.ravel()], order=3, mode="constant", cval=0.0)
        return (self.ovs ** 2) * val.reshape(ny, nx)

    def render_local(self, x: float, y: float, half: int = 3):
        ix, iy = int(np.floor(x + 0.5)), int(np.floor(y + 0.5))
        n = 2 * half + 1
        return self.render((n, n), x - (ix - half), y - (iy - half)), ix - half, iy - half

    def fwhm_px(self) -> float:
        """Approximate FWHM from the oversampled core (for display/diagnostics)."""
        prof = self.data[int(round(self.c)), :]
        half = prof.max() / 2
        above = np.where(prof >= half)[0]
        return float((above[-1] - above[0] + 1) / self.ovs)


def add_source(img: np.ndarray, psf, x: float, y: float, flux: float, half: int = 4) -> None:
    """Add ``flux`` (any linear unit) of a point source at (x, y) into ``img`` in place."""
    stamp, ix0, iy0 = psf.render_local(x, y, half)
    ny, nx = img.shape
    sy0, sx0 = max(0, -iy0), max(0, -ix0)
    y0, x0 = max(0, iy0), max(0, ix0)
    y1, x1 = min(ny, iy0 + stamp.shape[0]), min(nx, ix0 + stamp.shape[1])
    if y1 > y0 and x1 > x0:
        img[y0:y1, x0:x1] += flux * stamp[sy0:sy0 + (y1 - y0), sx0:sx0 + (x1 - x0)]
