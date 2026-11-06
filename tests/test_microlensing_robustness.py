import numpy as np
import pytest
from microlensing_robustness import pspl, blended_flux, fisher_blending, cadence_phase_completeness


def test_pspl_symmetric():
    assert np.allclose(pspl(np.array([-2., 2.]), 0, .2, 20),
                       pspl(np.array([2., -2.]), 0, .2, 20))


def test_blending_limit():
    t = np.linspace(-10, 10, 30)
    assert np.allclose(blended_flux(t, 0, .2, 20, 1), pspl(t, 0, .2, 20))


def test_fisher_positive_and_noise_scaling():
    t = np.linspace(-30, 30, 101)
    a = fisher_blending(t, .02)
    b = fisher_blending(t, .01)
    assert np.all(a["singular_values"] > 0)
    assert b["u0_sigma"] == pytest.approx(a["u0_sigma"]/2, rel=1e-8)


def test_phase_completeness():
    assert cadence_phase_completeness(1, 4, 3) == 1
    assert cadence_phase_completeness(1, .5, 3) == 0
    assert cadence_phase_completeness(1, 2.5, 3) == pytest.approx(.5, abs=.001)


def test_reject_bad_parameters():
    with pytest.raises(ValueError):
        pspl(np.array([0]), 0, 0, 1)
    with pytest.raises(ValueError):
        cadence_phase_completeness(0, 1)
