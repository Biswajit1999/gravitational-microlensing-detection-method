"""Synthetic cadence-phase and blended point-lens Fisher diagnostics.

No Gaussian perturbation is represented as a full binary-lens solution.
"""
import numpy as np


def pspl(t, t0, u0, t_e):
    if u0 <= 0 or t_e <= 0:
        raise ValueError("u0 and t_e must be positive")
    u = np.hypot(u0, (np.asarray(t, float) - t0) / t_e)
    return (u*u + 2) / (u * np.sqrt(u*u + 4))


def blended_flux(t, t0, u0, t_e, source_fraction):
    if not 0 < source_fraction <= 1:
        raise ValueError("source_fraction must be in (0,1]")
    return 1 + source_fraction * (pspl(t, t0, u0, t_e) - 1)


def fisher_blending(t, sigma, t0=0., u0=.2, t_e=20., source_fraction=.7):
    """Gaussian Fisher matrix for (u0, tE, source fraction), fixed t0.

    The reported covariance is conditional on known independent errors and
    a correct blended point-source point-lens model.
    """
    t = np.asarray(t, float)
    sigma = np.broadcast_to(np.asarray(sigma, float), t.shape)
    if t.ndim != 1 or t.size < 4 or not np.all(np.isfinite(sigma)) or np.any(sigma <= 0):
        raise ValueError("invalid sampling or noise")
    params = np.array([u0, t_e, source_fraction], float)
    steps = np.array([1e-5, 1e-3, 1e-5])
    cols = []
    for i in range(3):
        p1 = params.copy()
        p2 = params.copy()
        p1[i] += steps[i]
        p2[i] -= steps[i]
        cols.append((blended_flux(t, t0, *p1) - blended_flux(t, t0, *p2)) / (2 * steps[i]) / sigma)
    j = np.stack(cols, axis=1)
    fisher = j.T @ j
    singular = np.linalg.svd(fisher, compute_uv=False)
    covariance = np.linalg.pinv(fisher, rtol=1e-12)
    return {"singular_values": singular,
            "condition_number": float(singular[0] / singular[-1]),
            "u0_sigma": float(np.sqrt(covariance[0, 0])),
            "fisher": fisher}


def cadence_phase_completeness(cadence_days, width_days, threshold=3, phases=10000):
    """Phase-grid fraction with >= threshold samples inside a top-hat window.

    This geometric sampling diagnostic excludes noise, weather, and false alarms.
    """
    if cadence_days <= 0 or width_days <= 0 or phases < 2 or threshold < 1:
        raise ValueError("invalid cadence experiment")
    phase = (np.arange(phases) + .5) * cadence_days / phases
    count = (np.floor((width_days/2 - phase) / cadence_days)
             - np.ceil((-width_days/2 - phase) / cadence_days) + 1)
    return float(np.mean(count >= threshold))
