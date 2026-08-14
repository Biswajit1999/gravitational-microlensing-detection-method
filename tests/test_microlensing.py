"""Executable checks on the PSPL magnification formula and the blind
perturbation-detection logic (the part that fixes the earlier circular
validation, where the injected truth leaked into the search)."""

import numpy as np
import microlensing_demo as ml


def test_pspl_magnification_matches_analytic_formula():
    u0 = 0.12
    a = ml.pspl_magnification(np.array([0.0]), 0.0, u0, 20.0)[0]
    expected = (u0**2 + 2) / (u0 * np.sqrt(u0**2 + 4))
    assert np.isclose(a, expected, rtol=1e-10)


def test_pspl_magnification_decays_to_unity_far_from_peak():
    far = ml.pspl_magnification(np.array([1000.0]), 0.0, 0.12, 20.0)[0]
    assert np.isclose(far, 1.0, atol=1e-3)


def test_pspl_magnification_peaks_at_t0():
    time = np.linspace(-5, 5, 2001)
    a = ml.pspl_magnification(time, 0.0, 0.12, 20.0)
    assert time[np.argmax(a)] == 0.0 or abs(time[np.argmax(a)]) < 0.01


def test_planetary_perturbation_peak_matches_scaling_relation():
    q = 1.5e-3
    peak = ml.planetary_perturbation(np.array([5.0]), 5.0, 20.0, q)[0]
    expected = ml.AMP_COEFF * np.sqrt(q)
    assert np.isclose(peak, expected, rtol=1e-10)


def test_find_perturbation_window_detects_injected_bump_blind():
    # The detector is handed only a noisy residual with an injected bump --
    # no true time or mass ratio is passed in, matching how the real
    # pipeline in main() calls this function.
    rng = np.random.default_rng(7)
    time = np.arange(0, 40, 0.05)
    residual = rng.normal(0, 0.01, size=time.size)
    bump_mask = (time >= 20) & (time <= 21)
    residual[bump_mask] += 0.2
    resid_err = np.full_like(time, 0.01)

    window = ml.find_perturbation_window(time, residual, resid_err)
    assert window is not None
    win_start, win_end = window
    # The detected window should fall inside (or very close to) the true
    # injected bump, recovered purely from the data.
    assert 19.5 <= win_start <= 20.5
    assert 20.5 <= win_end <= 21.5


def test_find_perturbation_window_returns_none_on_pure_noise():
    rng = np.random.default_rng(3)
    time = np.arange(0, 40, 0.05)
    residual = rng.normal(0, 0.01, size=time.size)
    resid_err = np.full_like(time, 0.01)
    window = ml.find_perturbation_window(time, residual, resid_err, sigma_thresh=5.0, min_run=5)
    assert window is None
