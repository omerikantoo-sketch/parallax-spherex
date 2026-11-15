"""Reference targets named in the Space Apps 'Planet X and SPHEREx' challenge.

Astrometry is from the literature (values typed in, *not* fetched: SIMBAD was
unreachable during development).  Each entry records its source and an
uncertainty note.  Positions are propagated to an observation epoch with
astropy's rigorous space-motion model (``SkyCoord.apply_space_motion``),
including parallax when a distance is given.  The SPHEREx data themselves are
used to *test* these predictions (see ``science.astrometry``).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time


@dataclass(frozen=True)
class Target:
    slug: str
    name: str
    ra: float          # deg, ICRS, at ``epoch``
    dec: float
    epoch: float       # Julian year
    pmra: float = 0.0  # mas/yr (includes cos dec)
    pmdec: float = 0.0
    plx: float | None = None  # mas
    rv: float | None = None   # km/s
    kind: str = "star"
    note: str = ""
    source: str = ""

    def coord(self) -> SkyCoord:
        kw = {}
        if self.plx:
            kw["distance"] = Distance(parallax=self.plx * u.mas)
            kw["radial_velocity"] = (self.rv or 0.0) * u.km / u.s
        return SkyCoord(ra=self.ra * u.deg, dec=self.dec * u.deg, pm_ra_cosdec=self.pmra * u.mas / u.yr,
                        pm_dec=self.pmdec * u.mas / u.yr, obstime=Time(self.epoch, format="jyear"),
                        frame="icrs", **kw)

    def at(self, mjd: float | np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Barycentric (no parallactic shift) ICRS position at the given MJD(s)."""
        c = self.coord().apply_space_motion(new_obstime=Time(mjd, format="mjd"))
        return np.asarray(c.ra.deg), np.asarray(c.dec.deg)


TARGETS: dict[str, Target] = {t.slug: t for t in [
    Target("barnard", "Barnard's Star", 269.44850252543, 4.73942005, 2016.0, -801.551, 10362.394, 546.976, -110.1,
           kind="high proper-motion M dwarf", source="Gaia DR3 4472832130942575872",
           note="~10.4 arcsec/yr = 1.7 SPHEREx pixels/yr; very bright (K~4.5) so saturation/bloom flags expected."),
    Target("eps-indi-b", "Epsilon Indi B", 331.043833, -56.782694, 2000.0, 3981.6, -2466.7, 274.8, -40.0,
           kind="binary brown dwarf (T1+T6)", source="literature J2000 position, Hipparcos-system proper motion (approx.)",
           note="Position approximate (~1 arcsec); 402 arcsec from Epsilon Indi A (K~2), whose PSF wings/ghosts may intrude."),
    Target("wise-0855", "WISE 0855-0714", 133.795129, -7.245147, 2014.0, -8118.0, 673.0, 439.0, None,
           kind="Y dwarf (coldest known brown dwarf)", source="challenge coordinates; pm/plx Luhman & Esplin (2016), approx.",
           note="Reference epoch of challenge coordinates assumed 2014.0 (+/-1 yr -> +/-8 arcsec). "
                "Essentially invisible below ~3 um; only the longest SPHEREx channels (D5/D6) can detect it."),
    Target("rw-aur", "RW Aurigae", 76.956417, 30.401417, 2000.0, 3.0, -25.0, None, None,
           kind="variable T Tauri binary", source="literature J2000 (approx.)",
           note="Irregular variable with deep dimming episodes."),
    Target("r-dra", "R Draconis", 248.167375, 66.754917, 2000.0, -8.0, 9.0, None, None,
           kind="Mira variable (P~246 d)", source="literature J2000 (approx.)",
           note="Very bright Mira; expect saturation in short-wavelength channels. Near the north ecliptic pole."),
]}
