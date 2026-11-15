"""Static-scene modelling, forced photometry and astrometry on real frames.

* ``forced_photometry``  -- PSF-weighted flux of a point source at a fixed
  (or proper-motion-propagated) sky position in every frame.  This is the
  on-the-fly spectrum / light curve.
* ``static_catalog`` + ``subtract_static`` -- Model A: static point sources
  found in the v=0 stack, each with a *free amplitude per frame* (absorbing
  the SED sampled at each frame's wavelength and any variability), fitted
  jointly by linear least squares; returns residual frames for moving-object
  search.  Joint fitting (rather than sequential subtraction) keeps blended
  neighbours from leaking flux into each other.
* ``fit_astrometry`` -- per-frame PSF centroid fits of a bright star and a
  weighted least-squares fit of position offset + proper motion + parallax
  (observer = SPHEREx, using the spacecraft geocentric position in the
  header + Earth's barycentric position from astropy).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .frames import Frame, tangent_offsets
from .matched import PsiPhi, Trajectories, gather


def forced_photometry(frames: list[Frame], xi: float | np.ndarray, eta: float | np.ndarray, half: int = 3) -> dict:
    """PSF-weighted flux at tangent-plane offsets (scalar, or one per frame)."""
    xi = np.broadcast_to(np.asarray(xi, float), (len(frames),))
    eta = np.broadcast_to(np.asarray(eta, float), (len(frames),))
    out = {k: [] for k in ("t", "wave", "band", "flux", "err", "det", "label", "x", "y", "nbad", "chi2")}
    for f, a, b in zip(frames, xi, eta):
        x, y = f.to_pix(a, b)
        x, y = float(x), float(y)
        stamp, ix0, iy0 = f.psf.render_local(x, y, half)
        ny, nx = f.shape
        ys, xs = slice(max(0, iy0), min(ny, iy0 + stamp.shape[0])), slice(max(0, ix0), min(nx, ix0 + stamp.shape[1]))
        if ys.stop <= ys.start or xs.stop <= xs.start:
            continue
        P = stamp[ys.start - iy0: ys.stop - iy0, xs.start - ix0: xs.stop - ix0]
        d = f.data[ys, xs]
        w = f.ivar[ys, xs]
        core = P > 0.02
        nbad = int(np.sum((w <= 0) & core))
        den = float(np.sum(w * P * P))
        if den <= 0 or nbad > 2:
            continue
        flux = float(np.sum(w * P * d) / den)
        chi2 = float(np.sum(w * (d - flux * P) ** 2))
        for k, v in zip(out, (f.t_mjd, f.wave_um, f.band_um, flux, 1 / np.sqrt(den), f.detector, f.label, x, y, nbad,
                              chi2)):
            out[k].append(v)
    return {k: np.array(v) if k != "label" else v for k, v in out.items()}


def static_catalog(maps: list[PsiPhi], snr: float = 5.0, max_n: int = 200, sep_arcsec: float = 9.0):
    """Detect static sources on the v=0 stack (constant flux statistic on the
    reference plane).  Returns arrays (xi, eta, snr)."""
    f0 = maps[0].frame
    pix = f0.pixscale
    ny, nx = f0.shape
    half = 0.5 * min(ny, nx) * pix
    step = 0.5 * pix
    g = np.arange(-half, half + step, step)
    X, Y = np.meshgrid(g, g)
    traj = Trajectories(X.ravel(), Y.ravel(), np.zeros(X.size), np.zeros(X.size), f0.t_mjd)
    psi, phi = gather(maps, traj)
    # per-frame free amplitudes would be chi2_E; for *detection* use the
    # per-detector chi2 combined as Gaussian-equivalent z (robust to SED).
    from .matched import chi2_to_z, stat_per_detector
    det = np.array([m.frame.detector for m in maps])
    D, k = stat_per_detector(psi, phi, det)
    z = chi2_to_z(D, k).reshape(X.shape)
    from scipy.ndimage import maximum_filter
    pk = (z == maximum_filter(z, size=5)) & (z >= snr)
    iy, ix = np.nonzero(pk)
    order = np.argsort(z[iy, ix])[::-1]
    xs, ys, zs = [], [], []
    for i in order:
        a, b = X[iy[i], ix[i]], Y[iy[i], ix[i]]
        if any(np.hypot(a - p, b - q) < sep_arcsec for p, q in zip(xs, ys)):
            continue
        xs.append(a); ys.append(b); zs.append(z[iy[i], ix[i]])
        if len(xs) >= max_n:
            break
    return np.array(xs), np.array(ys), np.array(zs)


def subtract_static(frames: list[Frame], xi: np.ndarray, eta: np.ndarray, half: int = 4) -> tuple[list[Frame], np.ndarray]:
    """Joint per-frame linear fit of static PSF amplitudes + constant; returns
    residual frames and the amplitude matrix (frames x sources)."""
    out = []
    amps = np.zeros((len(frames), len(xi)))
    for e, f in enumerate(frames):
        ny, nx = f.shape
        cols = []
        for a, b in zip(xi, eta):
            x, y = f.to_pix(a, b)
            img = np.zeros(f.shape)
            from ..science.psf import add_source
            add_source(img, f.psf, float(x), float(y), 1.0, half=half)
            cols.append(img.ravel())
        cols.append(np.ones(ny * nx))
        Xd = np.stack(cols, axis=1)
        w = f.ivar.ravel()
        sel = w > 0
        if sel.sum() < Xd.shape[1] + 5:
            out.append(f)
            continue
        Xw = Xd[sel] * np.sqrt(w[sel])[:, None]
        yw = f.data.ravel()[sel] * np.sqrt(w[sel])
        coef, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
        model = (Xd @ coef).reshape(ny, nx)
        nf = Frame(**{**f.__dict__})
        nf.data = np.where(f.ivar > 0, f.data - model, 0.0)
        nf.meta = {**f.meta, "static_model": "joint PSF fit, %d sources" % len(xi)}
        out.append(nf)
        amps[e] = coef[:-1]
    return out, amps


# ----------------------------------------------------------------------------
# astrometry
# ----------------------------------------------------------------------------
def psf_centroid(f: Frame, x0: float, y0: float, half: int = 3):
    """Least-squares PSF fit (flux, x, y, background) near (x0, y0)."""
    ny, nx = f.shape
    ix, iy = int(round(x0)), int(round(y0))
    ys, xs = slice(max(0, iy - half), min(ny, iy + half + 1)), slice(max(0, ix - half), min(nx, ix + half + 1))
    d = f.data[ys, xs]
    w = f.ivar[ys, xs]
    if np.sum(w > 0) < 12:
        return None

    def model(p):
        flux, x, y, bg = p
        P = f.psf.render(d.shape, x - xs.start, y - ys.start)
        return flux * P + bg

    def resid(p):
        return ((d - model(p)) * np.sqrt(w)).ravel()

    f0 = max(float(np.sum(d[w > 0])), 1.0)
    try:
        r = least_squares(resid, [f0, x0, y0, 0.0], x_scale=[f0, 0.1, 0.1, max(1.0, f0 * 0.01)])
    except Exception:
        return None
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J) * max(1.0, 2 * r.cost / max(1, J.shape[0] - 4))
    except np.linalg.LinAlgError:
        return None
    return {"flux": r.x[0], "x": r.x[1], "y": r.x[2], "sx": float(np.sqrt(cov[1, 1])), "sy": float(np.sqrt(cov[2, 2])),
            "chi2r": float(2 * r.cost / max(1, J.shape[0] - 4))}


def parallax_factors(mjd: np.ndarray, ra: float, dec: float, sc_geo_km: np.ndarray | None = None):
    """(P_alpha*, P_delta) such that the parallactic offset is plx * P (Green 1985),
    for an observer at Earth's barycentric position plus the spacecraft offset."""
    from astropy.coordinates import get_body_barycentric
    from astropy.time import Time
    t = Time(np.asarray(mjd), format="mjd", scale="utc")
    e = get_body_barycentric("earth", t)
    X = e.x.to_value("au")
    Y = e.y.to_value("au")
    Z = e.z.to_value("au")
    if sc_geo_km is not None:
        au = 1.495978707e8
        X = X + sc_geo_km[:, 0] / au
        Y = Y + sc_geo_km[:, 1] / au
        Z = Z + sc_geo_km[:, 2] / au
    a, d = np.radians(ra), np.radians(dec)
    pa = X * np.sin(a) - Y * np.cos(a)
    pd = X * np.cos(a) * np.sin(d) + Y * np.sin(a) * np.sin(d) - Z * np.cos(d)
    return pa, pd


def fit_astrometry(t: np.ndarray, xi: np.ndarray, eta: np.ndarray, sxi: np.ndarray, seta: np.ndarray,
                   ra0: float, dec0: float, t_ref: float, fit_parallax: bool = False, plx_fixed_mas: float | None = None,
                   sc_geo_km: np.ndarray | None = None) -> dict:
    """Weighted LSQ for offsets (arcsec) + proper motion (mas/yr) [+ parallax]."""
    dt = (t - t_ref) / 365.25
    pa, pd = parallax_factors(t, ra0, dec0, sc_geo_km)
    rows, ys, ws = [], [], []
    for i in range(t.size):
        r1 = [1, 0, dt[i], 0]
        r2 = [0, 1, 0, dt[i]]
        if fit_parallax:
            r1.append(pa[i]); r2.append(pd[i])
        rows += [r1, r2]
        yx, yy = xi[i], eta[i]
        if plx_fixed_mas is not None and not fit_parallax:
            yx -= plx_fixed_mas / 1000 * pa[i]
            yy -= plx_fixed_mas / 1000 * pd[i]
        ys += [yx, yy]
        ws += [1 / sxi[i] ** 2, 1 / seta[i] ** 2]
    X, y, w = np.array(rows, float), np.array(ys), np.array(ws)
    G = X.T @ (X * w[:, None])
    coef = np.linalg.solve(G, X.T @ (w * y))
    cov = np.linalg.inv(G)
    res = y - X @ coef
    chi2 = float(np.sum(w * res ** 2))
    dof = max(1, y.size - X.shape[1])
    scale = max(1.0, chi2 / dof)
    err = np.sqrt(np.diag(cov) * scale)
    out = {"xi0_arcsec": coef[0], "eta0_arcsec": coef[1], "pmra_masyr": coef[2] * 1000, "pmdec_masyr": coef[3] * 1000,
           "e_pmra_masyr": err[2] * 1000, "e_pmdec_masyr": err[3] * 1000, "chi2": chi2, "dof": dof,
           "error_scale": float(np.sqrt(scale)), "n": int(t.size), "t_ref_mjd": t_ref,
           "residual_rms_arcsec": float(np.sqrt(np.mean(res ** 2)))}
    if fit_parallax:
        out["plx_mas"] = coef[4] * 1000
        out["e_plx_mas"] = err[4] * 1000
    return out
