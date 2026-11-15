"""Baselines and a multi-statistic search used for controlled comparisons.

All methods consume the *same* psi/phi gathers on the *same* trajectory grid,
so differences between them come only from how flux is modelled across frames.

Methods
-------
link        conventional: per-frame detection (matched-filter SNR >= 5) and
            linking of >= ``min_det`` detections along a linear trajectory.
stack       constant-flux shift-and-stack over every frame (KBMOD statistic).
best_det    shift-and-stack using only the single most sensitive detector.
perdet      independent constant flux per detector, chi^2_Ndet combination.
template    known source SED (oracle upper reference, not a fair method).
bank        max over a physically motivated SED template bank.
subspace2   unknown SED in a 2-D smooth basis (chi^2_2).
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import maximum_filter

from ..science.spectra import spectral_basis, template_bank
from .matched import (PsiPhi, SearchGrid, Trajectories, chi2_to_z, gather, spectral_factors, stat_bank,
                      stat_constant, stat_per_detector, stat_subspace, stat_template)

METHODS = ["stack", "best_det", "perdet", "bank", "subspace2", "template"]


def statistics_for(psi, phi, det, S_bank, B2, s_true=None, best_detector=None) -> dict[str, np.ndarray]:
    out = {}
    out["stack"], _ = stat_constant(psi, phi)
    if best_detector is not None:
        m = det == best_detector
        out["best_det"], _ = stat_constant(psi[m], phi[m])
    D, k = stat_per_detector(psi, phi, det)
    out["perdet"] = chi2_to_z(D, k)
    out["bank"], _, _ = stat_bank(psi, phi, S_bank)
    D2, _ = stat_subspace(psi, phi, B2)
    out["subspace2"] = chi2_to_z(D2, B2.shape[1])
    if s_true is not None:
        out["template"], _ = stat_template(psi, phi, s_true)
    return out


def setup(maps: list[PsiPhi], spectrum_true=None):
    frames = [m.frame for m in maps]
    det = np.array([f.detector for f in frames])
    bank = template_bank()
    S = np.stack([spectral_factors(frames, sp) for sp in bank])
    B2 = spectral_basis(np.array([f.wave_um for f in frames]), 2)
    s_true = spectral_factors(frames, spectrum_true) if spectrum_true is not None else None
    # most sensitive detector for a solar-type spectrum: max sum(s^2 / sigma^2)
    sol = spectral_factors(frames, bank[1])
    sens = {d: float(np.sum([sol[i] ** 2 * np.median(frames[i].ivar[frames[i].ivar > 0])
                             for i in range(len(frames)) if det[i] == d])) for d in np.unique(det)}
    best = max(sens, key=sens.get)
    return det, S, B2, s_true, best


def multi_search(maps: list[PsiPhi], grid: SearchGrid, spectrum_true=None) -> dict[str, float]:
    """Global maximum of each statistic over the full trajectory grid."""
    det, S, B2, s_true, best = setup(maps, spectrum_true)
    px, py = grid.positions()
    vx, vy = grid.velocities()
    best_val: dict[str, float] = {}
    for iv in range(vx.size):
        traj = Trajectories(px, py, np.full(px.size, vx[iv]), np.full(px.size, vy[iv]), grid.t_ref)
        psi, phi = gather(maps, traj)
        st = statistics_for(psi, phi, det, S, B2, s_true, best)
        for k, v in st.items():
            m = float(np.max(v))
            if m > best_val.get(k, -np.inf):
                best_val[k] = m
    return best_val


def local_stats(maps: list[PsiPhi], truth_xi, truth_eta, truth_vxi, truth_veta, t_ref, spectrum_true=None,
                pos_tol=3.0, v_steps=(-1, 0, 1), v_step=1.0, pos_step=3.0) -> dict[str, float]:
    """Maximum of each statistic over grid trajectories adjacent to the truth
    (what a grid search would report for the injected object)."""
    det, S, B2, s_true, best = setup(maps, spectrum_true)
    offs = np.array([-1, 0, 1]) * pos_step
    out: dict[str, float] = {}
    for dvx in v_steps:
        for dvy in v_steps:
            X, Y = np.meshgrid(truth_xi + offs, truth_eta + offs)
            n = X.size
            traj = Trajectories(X.ravel(), Y.ravel(), np.full(n, truth_vxi + dvx * v_step),
                                np.full(n, truth_veta + dvy * v_step), t_ref)
            psi, phi = gather(maps, traj)
            st = statistics_for(psi, phi, det, S, B2, s_true, best)
            for k, v in st.items():
                out[k] = max(out.get(k, -np.inf), float(np.max(v)))
    return out


# ----------------------------------------------------------------------------
# conventional detect-then-link baseline
# ----------------------------------------------------------------------------
def frame_detections(m: PsiPhi, snr_thresh: float = 5.0):
    snr = m.psi / np.sqrt(np.maximum(m.phi, 1e-30))
    snr[m.phi <= 0] = 0
    peaks = (snr == maximum_filter(snr, size=2 * m.q + 1)) & (snr >= snr_thresh)
    iy, ix = np.nonzero(peaks)
    return ix / m.q, iy / m.q, snr[iy, ix]


def link_count_at_truth(maps: list[PsiPhi], truth_xi, truth_eta, truth_vxi, truth_veta, t_ref,
                        snr_thresh: float = 5.0, tol_px: float = 1.5) -> int:
    """Number of distinct *exposures* with a >=snr_thresh detection within
    tol_px of the true predicted position (necessary condition for linking)."""
    exps = set()
    for m in maps:
        dt = m.frame.t_mjd - t_ref
        x0, y0 = m.frame.to_pix(truth_xi + truth_vxi * dt, truth_eta + truth_veta * dt)
        x, y, _ = frame_detections(m, snr_thresh)
        if x.size and np.min(np.hypot(x - x0, y - y0)) <= tol_px:
            exps.add(round(m.frame.t_mjd, 5))
    return len(exps)


def link_false_tracklets(maps: list[PsiPhi], vmax_arcsec_day: float, snr_thresh: float = 5.0,
                         tol_arcsec: float = 9.0, min_det: int = 3) -> int:
    """Count linear tracklets with >= min_det exposures built purely from
    detections in (noise-only) frames; used to set the link baseline's false
    alarm rate.  Detections are expressed on the common tangent plane."""
    pts = []  # (t, xi, eta)
    for m in maps:
        x, y, _ = frame_detections(m, snr_thresh)
        if x.size == 0:
            continue
        Ainv = np.linalg.inv(m.frame.A)
        xy = np.stack([x - m.frame.b[0], y - m.frame.b[1]])
        xe = Ainv @ xy
        for i in range(x.size):
            pts.append((m.frame.t_mjd, xe[0, i], xe[1, i]))
    if len(pts) < min_det:
        return 0
    P = np.array(pts)
    n_track = 0
    used = set()
    for i in range(len(P)):
        for j in range(i + 1, len(P)):
            dt = P[j, 0] - P[i, 0]
            if abs(dt) < 1e-3:
                continue
            v = (P[j, 1:] - P[i, 1:]) / dt
            if np.hypot(*v) > vmax_arcsec_day:
                continue
            pred = P[i, 1:] + np.outer(P[:, 0] - P[i, 0], v)
            ok = np.hypot(*(P[:, 1:] - pred).T) < tol_arcsec
            nexp = len(set(np.round(P[ok, 0], 5)))
            if nexp >= min_det:
                key = tuple(sorted(np.where(ok)[0][:3]))
                if key not in used:
                    used.add(key)
                    n_track += 1
    return n_track
