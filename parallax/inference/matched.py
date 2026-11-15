r"""Matched-filter sufficient statistics and detection statistics for moving
point sources observed in many frames at different wavelengths.

Model (per frame e, pixel p), Gaussian noise with inverse variance w_{e,p}:

    d_{e,p} = f_e P_e(p - x_e(theta)) + n_{e,p}

For a fixed trajectory theta, the log-likelihood is quadratic in the per-frame
fluxes f_e and depends on the data only through

    psi_e = sum_p w P d          phi_e = sum_p w P^2

(Kaiser 2004; Whidden et al. 2019 "KBMOD").  Different flux models then give
different statistics, all with closed forms:

* constant flux  f_e = A                     ->  S = sum psi / sqrt(sum phi)          (N(0,1) null)
* known spectrum f_e = A s_e                 ->  S = sum s psi / sqrt(sum s^2 phi)    (N(0,1) null)
* linear spectral subspace f = B a  (K dims)  ->  D = psi^T B (B^T Phi B)^{-1} B^T psi (chi^2_K null)
* independent per-detector constants          ->  D = sum_d (sum psi)^2 / sum phi    (chi^2_{Ndet})
* free per-frame fluxes                       ->  D = sum psi^2 / phi                (chi^2_E)

The spectral factor s_e is the source SED evaluated at frame e's wavelength,
which for SPHEREx differs between detectors *and* between exposures.

The per-frame maps psi_e(x), phi_e(x) are pre-computed by correlating the
weighted data with the PSF evaluated at q x q sub-pixel phases, so evaluating a
trajectory is a gather of E numbers ("shift-and-stack" on psi/phi images).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import correlate, maximum_filter, minimum_filter
from scipy.stats import chi2, norm

from .frames import Frame


# =============================================================================
# psi / phi maps
# =============================================================================
@dataclass
class PsiPhi:
    frame: Frame
    q: int = 2
    half: int = 4
    psi: np.ndarray = field(init=False)
    phi: np.ndarray = field(init=False)

    def __post_init__(self):
        f = self.frame
        ny, nx = f.shape
        q = self.q
        n = 2 * self.half + 1
        wd = f.data * f.ivar
        self.psi = np.zeros((ny * q, nx * q), dtype=np.float32)
        self.phi = np.zeros((ny * q, nx * q), dtype=np.float32)
        for iy in range(q):
            for ix in range(q):
                st = f.psf.render((n, n), self.half + ix / q, self.half + iy / q)
                self.psi[iy::q, ix::q] = correlate(wd, st, mode="constant", cval=0.0)
                self.phi[iy::q, ix::q] = correlate(f.ivar, st * st, mode="constant", cval=0.0)
        self._pool: dict[int, tuple[np.ndarray, np.ndarray]] = {}

    def index(self, x, y):
        """Nearest q-grid indices for pixel coords; returns (iy, ix, valid)."""
        q = self.q
        ix = np.rint(np.asarray(x) * q).astype(np.int64)
        iy = np.rint(np.asarray(y) * q).astype(np.int64)
        ny, nx = self.psi.shape
        ok = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)
        return np.clip(iy, 0, ny - 1), np.clip(ix, 0, nx - 1), ok

    def sample(self, x, y):
        iy, ix, ok = self.index(x, y)
        return np.where(ok, self.psi[iy, ix], 0.0), np.where(ok, self.phi[iy, ix], 0.0)

    def pooled(self, r: int) -> tuple[np.ndarray, np.ndarray]:
        """(max-pooled psi, min-pooled phi) over (2r+1)^2 windows.  Outside the
        frame psi=phi=0 is a legal value (no data), so padding uses 0."""
        if r not in self._pool:
            if r == 0:
                self._pool[r] = (self.psi, self.phi)
            else:
                k = 2 * r + 1
                self._pool[r] = (maximum_filter(self.psi, size=k, mode="constant", cval=0.0),
                                 minimum_filter(self.phi, size=k, mode="constant", cval=0.0))
        return self._pool[r]

    def sample_pooled(self, x, y, r: int):
        mp, mf = self.pooled(r)
        iy, ix, ok = self.index(x, y)
        # if the window centre is off-frame but within r, the window may still
        # overlap data; be conservative: for off-frame centres use clipped index
        # (the pooled value at the clipped index covers the overlap) only when
        # within r q-pixels of the frame, else zero.
        ny, nx = mp.shape
        rx = np.rint(np.asarray(x) * self.q).astype(np.int64)
        ry = np.rint(np.asarray(y) * self.q).astype(np.int64)
        near = (rx >= -r) & (rx < nx + r) & (ry >= -r) & (ry < ny + r)
        pmax = np.where(near, mp[iy, ix], 0.0)
        # clipped index may not see the whole window -> use global max for edge cases
        edge = near & ~ok
        if np.any(edge):
            pmax = np.where(edge, max(float(self.psi.max()), 0.0), pmax)
        pmin = np.where(ok, mf[iy, ix], 0.0)
        return pmax, pmin


def build_maps(frames: list[Frame], q: int = 2) -> list[PsiPhi]:
    return [PsiPhi(f, q=q) for f in frames]


# =============================================================================
# trajectories
# =============================================================================
@dataclass
class Trajectories:
    xi0: np.ndarray    # arcsec at t_ref
    eta0: np.ndarray
    vxi: np.ndarray    # arcsec/day
    veta: np.ndarray
    t_ref: float

    def __len__(self):
        return int(np.size(self.xi0))

    def at(self, t):
        dt = t - self.t_ref
        return self.xi0 + self.vxi * dt, self.eta0 + self.veta * dt


def gather(maps: list[PsiPhi], traj: Trajectories) -> tuple[np.ndarray, np.ndarray]:
    """psi, phi arrays of shape (E, N)."""
    E, N = len(maps), len(traj)
    psi = np.empty((E, N), dtype=np.float64)
    phi = np.empty((E, N), dtype=np.float64)
    for e, m in enumerate(maps):
        xi, eta = traj.at(m.frame.t_mjd)
        x, y = m.frame.to_pix(xi, eta)
        psi[e], phi[e] = m.sample(x, y)
    return psi, phi


# =============================================================================
# statistics
# =============================================================================
def stat_constant(psi, phi):
    den = np.sqrt(np.maximum(phi.sum(0), 1e-30))
    return psi.sum(0) / den, psi.sum(0) / np.maximum(phi.sum(0), 1e-30)


def stat_template(psi, phi, s):
    s = np.asarray(s, float)[:, None]
    num = (s * psi).sum(0)
    den = (s * s * phi).sum(0)
    return num / np.sqrt(np.maximum(den, 1e-30)), num / np.maximum(den, 1e-30)


def stat_bank(psi, phi, S):
    """Max over a template bank S (T, E).  Returns (stat, best template index, amplitude)."""
    num = S @ psi                      # (T, N)
    den = (S * S) @ phi
    z = num / np.sqrt(np.maximum(den, 1e-30))
    k = np.argmax(z, axis=0)
    idx = np.arange(psi.shape[1])
    return z[k, idx], k, (num / np.maximum(den, 1e-30))[k, idx]


def stat_subspace(psi, phi, B, nonneg_first: bool = True):
    """Delta chi^2 for flux = B a (B: (E, K)).  chi^2_K under the null."""
    E, K = B.shape
    b = B.T @ psi                                       # (K, N)
    G = np.einsum("ek,el,en->nkl", B, B, phi)           # (N, K, K)
    G += np.eye(K)[None] * 1e-12 * np.maximum(np.trace(G, axis1=1, axis2=2), 1e-30)[:, None, None]
    a = np.linalg.solve(G, b.T[..., None])[..., 0]      # (N, K)
    D = np.einsum("nk,kn->n", a, b)
    if nonneg_first:  # a physical source has positive flux at the envelope
        D = np.where(a[:, 0] > 0, D, 0.0)
    return D, a


def stat_per_detector(psi, phi, det):
    det = np.asarray(det)
    D = 0.0
    for d in np.unique(det):
        m = det == d
        D = D + psi[m].sum(0) ** 2 / np.maximum(phi[m].sum(0), 1e-30)
    return D, len(np.unique(det))


def stat_free(psi, phi):
    return (psi ** 2 / np.maximum(phi, 1e-30)).sum(0)


def chi2_to_z(D, k):
    """Gaussian-equivalent one-sided significance of a chi^2_k value.

    Exact (via the survival function) where representable; for extreme values
    (p < 1e-250, z > ~33) the Wilson-Hilferty cube-root transform is used so the
    statistic keeps increasing with D instead of saturating -- a saturated z
    made bright static stars indistinguishable in the static catalogue.
    """
    D = np.maximum(np.asarray(D, float), 0.0)
    k = np.asarray(k, float)
    p = chi2.sf(D, k)
    z = norm.isf(np.clip(p, 1e-300, 1.0))
    wh = ((D / k) ** (1.0 / 3.0) - (1.0 - 2.0 / (9.0 * k))) / np.sqrt(2.0 / (9.0 * k))
    return np.where(p < 1e-250, np.maximum(wh, z), z)


# =============================================================================
# exhaustive grid search
# =============================================================================
@dataclass
class SearchGrid:
    xi: np.ndarray       # start positions (1-D arrays, outer product forms the grid)
    eta: np.ndarray
    vxi: np.ndarray
    veta: np.ndarray
    t_ref: float

    def positions(self):
        X, Y = np.meshgrid(self.xi, self.eta)
        return X.ravel(), Y.ravel()

    def velocities(self):
        VX, VY = np.meshgrid(self.vxi, self.veta)
        return VX.ravel(), VY.ravel()

    @property
    def n_trials(self):
        return self.xi.size * self.eta.size * self.vxi.size * self.veta.size


def make_grid(frames: list[Frame], half_fov_arcsec: float, vmax_arcsec_day: float,
              pos_step: float | None = None, v_step: float | None = None, vmin_abs: float = 0.0) -> SearchGrid:
    """Grid sized so that adjacent trajectories differ by <= ~half a pixel at
    every epoch (pos_step = 0.5 px, v_step = 0.5 px / max|t - t_ref|)."""
    t = np.array([f.t_mjd for f in frames])
    t_ref = float(0.5 * (t.min() + t.max()))
    pix = float(np.median([f.pixscale for f in frames]))
    pos_step = pos_step or 0.5 * pix
    span = max(float(np.max(np.abs(t - t_ref))), 1e-3)
    v_step = v_step or 0.5 * pix / span
    nv = int(np.ceil(vmax_arcsec_day / v_step))
    v = np.arange(-nv, nv + 1) * v_step
    npos = int(np.ceil(half_fov_arcsec / pos_step))
    p = np.arange(-npos, npos + 1) * pos_step
    return SearchGrid(p, p, v, v, t_ref)


@dataclass
class Candidate:
    xi0: float
    eta0: float
    vxi: float
    veta: float
    stat: float
    flux: float
    method: str
    extra: dict = field(default_factory=dict)


def spectral_factors(frames: list[Frame], spectrum) -> np.ndarray:
    from ..science.spectra import band_average
    return np.array([float(band_average(spectrum, f.wave_um, f.band_um)) for f in frames])


def exhaustive_search(maps: list[PsiPhi], grid: SearchGrid, method: str = "constant",
                      S: np.ndarray | None = None, B: np.ndarray | None = None,
                      chunk: int = 200_000, keep: int = 50) -> tuple[list[Candidate], dict]:
    """Evaluate every trajectory in ``grid``.  Returns top candidates (after a
    simple velocity-position non-maximum suppression) and the full max map."""
    px, py = grid.positions()
    vx, vy = grid.velocities()
    det = np.array([m.frame.detector for m in maps])
    best = np.full(px.size, -np.inf)
    best_v = np.zeros(px.size, dtype=np.int64)
    best_amp = np.zeros(px.size)
    n_eval = 0
    for iv in range(vx.size):
        traj = Trajectories(px, py, np.full(px.size, vx[iv]), np.full(px.size, vy[iv]), grid.t_ref)
        psi, phi = gather(maps, traj)
        n_eval += px.size
        amp = np.zeros(px.size)
        if method == "constant":
            z, amp = stat_constant(psi, phi)
        elif method == "template":
            z, amp = stat_template(psi, phi, S[0])
        elif method == "bank":
            z, _, amp = stat_bank(psi, phi, S)
        elif method == "subspace":
            D, a = stat_subspace(psi, phi, B)
            z, amp = chi2_to_z(D, B.shape[1]), a[:, 0]
        elif method == "perdet":
            D, k = stat_per_detector(psi, phi, det)
            z = chi2_to_z(D, k)
        elif method == "free":
            z = chi2_to_z(stat_free(psi, phi), len(maps))
        else:
            raise ValueError(method)
        upd = z > best
        best[upd] = z[upd]
        best_v[upd] = iv
        best_amp[upd] = amp[upd] if np.ndim(amp) else 0.0
    order = np.argsort(best)[::-1]
    cands: list[Candidate] = []
    pix = float(np.median([m.frame.pixscale for m in maps]))
    for i in order:
        if len(cands) >= keep:
            break
        c = Candidate(float(px[i]), float(py[i]), float(vx[best_v[i]]), float(vy[best_v[i]]), float(best[i]),
                      float(best_amp[i]), method)
        if any(np.hypot(c.xi0 - o.xi0, c.eta0 - o.eta0) < 1.5 * pix for o in cands):
            continue
        cands.append(c)
    return cands, {"n_eval": n_eval, "max_map": best.reshape(grid.eta.size, grid.xi.size)}


# =============================================================================
# branch-and-bound over velocity (template statistics)
# =============================================================================
def branch_and_bound(maps: list[PsiPhi], grid: SearchGrid, s: np.ndarray, threshold: float,
                     coarse_levels: int = 3) -> tuple[list[tuple[int, int, float]], dict]:
    r"""Find *all* (position, velocity) grid trajectories with template SNR >=
    ``threshold`` using velocity-cell upper bounds.

    For a velocity cell C (centre v_c, half-widths h) and frame e the predicted
    positions of every trajectory in C lie within a box of half-size
    |t_e - t_ref| * h around the centre prediction.  With s_e >= 0,

        SNR(theta) <= sum_e s_e max_{box_e} psi_e / sqrt( sum_e s_e^2 min_{box_e} phi_e )

    whenever the numerator bound is positive (if it is <= 0 no trajectory in C
    can reach a positive threshold).  Positions are rounded to the q-grid, so
    the box is enlarged by one q-pixel.  The bound is therefore rigorous for
    the discrete trajectory set evaluated by ``exhaustive_search``.
    """
    assert np.all(s >= 0) and threshold > 0
    px, py = grid.positions()
    vx, vy = grid.vxi, grid.veta
    dvx = float(vx[1] - vx[0]) if vx.size > 1 else 1.0
    dvy = float(vy[1] - vy[0]) if vy.size > 1 else 1.0
    # root cells: split the velocity index ranges into 2^L blocks
    nvx, nvy = vx.size, vy.size
    cells = [(0, nvx, 0, nvy, np.arange(px.size))]
    n_bound = 0
    n_exact = 0
    found: list[tuple[int, int, float]] = []
    s2 = s * s
    while cells:
        new = []
        for (ix0, ix1, iy0, iy1, pos_idx) in cells:
            if pos_idx.size == 0:
                continue
            if (ix1 - ix0) * (iy1 - iy0) <= 1:
                # leaf: exact evaluation at this grid velocity
                traj = Trajectories(px[pos_idx], py[pos_idx], np.full(pos_idx.size, vx[ix0]),
                                    np.full(pos_idx.size, vy[iy0]), grid.t_ref)
                psi, phi = gather(maps, traj)
                z, _ = stat_template(psi, phi, s)
                n_exact += pos_idx.size
                for j in np.where(z >= threshold)[0]:
                    found.append((int(pos_idx[j]), int(iy0 * nvx + ix0), float(z[j])))
                continue
            # cell centre and half-widths in velocity (covering grid points ix0..ix1-1)
            vcx = 0.5 * (vx[ix0] + vx[ix1 - 1])
            vcy = 0.5 * (vy[iy0] + vy[iy1 - 1])
            hx = 0.5 * (vx[ix1 - 1] - vx[ix0])
            hy = 0.5 * (vy[iy1 - 1] - vy[iy0])
            num = np.zeros(pos_idx.size)
            den = np.zeros(pos_idx.size)
            for e, m in enumerate(maps):
                dt = m.frame.t_mjd - grid.t_ref
                x, y = m.frame.to_pix(px[pos_idx] + vcx * dt, py[pos_idx] + vcy * dt)
                A = np.abs(m.frame.A)
                rx = abs(dt) * (A[0, 0] * hx + A[0, 1] * hy)
                ry = abs(dt) * (A[1, 0] * hx + A[1, 1] * hy)
                r = int(np.ceil(m.q * max(rx, ry))) + 1
                pmax, pmin = m.sample_pooled(x, y, r)
                num += s[e] * pmax
                den += s2[e] * pmin
            n_bound += pos_idx.size
            ub = np.where(num > 0, num / np.sqrt(np.maximum(den, 1e-30)), -np.inf)
            ub = np.where((num > 0) & (den <= 0), np.inf, ub)
            keep_idx = pos_idx[ub >= threshold]
            if keep_idx.size == 0:
                continue
            # split the larger velocity dimension
            if (ix1 - ix0) >= (iy1 - iy0):
                mid = (ix0 + ix1) // 2
                new += [(ix0, mid, iy0, iy1, keep_idx), (mid, ix1, iy0, iy1, keep_idx)]
            else:
                mid = (iy0 + iy1) // 2
                new += [(ix0, ix1, iy0, mid, keep_idx), (ix0, ix1, mid, iy1, keep_idx)]
        cells = new
    full = px.size * vx.size * vy.size
    return found, {"n_bound_evals": n_bound, "n_exact_evals": n_exact, "n_exhaustive": full,
                   "work_fraction": (n_bound + n_exact) / full}
