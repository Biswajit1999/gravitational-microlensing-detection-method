"""Gravitational microlensing method demonstration: simulate a single-
lens point-source magnification event with a short planetary
perturbation superimposed, using the real characteristic scaling
relations from planetary microlensing theory, then recover the host-
lens Einstein-time and the planet-to-star mass ratio q from the light
curve.

This is a PEDAGOGICAL DEMONSTRATION with simulated data, not a specific
real target's raw archival light curve (see README.md for why, and see
this portfolio's *-exoplanet-report repos for 11 planets analyzed
directly from real archival JWST/HST/Spitzer/ground-based data). Full
binary-lens light curves require solving for image positions via
inverse ray-shooting or complex polynomial root-finding; this repo
instead uses the well-established analytic scaling relations for a
planetary perturbation's characteristic amplitude and duration relative
to the mass ratio q (Gould & Loeb 1992; Mao & Paczynski 1991), which is
the standard back-of-the-envelope real-world tool for estimating a
detected planet's mass ratio and is explicitly labeled as an
approximation, not a full caustic computation.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import scienceplots  # noqa: F401 (registers 'science' style)
import numpy as np
from scipy.optimize import curve_fit

plt.style.use(["science", "no-latex"])

FIG_DIR = Path(__file__).resolve().parents[1] / "figures"

rng = np.random.default_rng(seed=3)

# Injected "ground truth" event parameters (realistic Galactic-bulge-survey regime).
TRUE_T0_DAYS = 25.0
TRUE_U0 = 0.12
TRUE_TE_DAYS = 20.0
TRUE_Q = 1.5e-3  # planet/star mass ratio, roughly a Jupiter-mass planet around a 0.5 Msun star
TRUE_TP_DAYS = 22.0  # time of closest planet-source approach

AMP_COEFF = 8.0  # order-unity scaling coefficient, Gould & Loeb (1992)-style
DUR_COEFF = 0.5

BASELINE_DAYS = 60.0
CADENCE_HR = 0.5  # real-like high-cadence survey coverage (e.g. KMTNet-class)
PHOT_NOISE_FRAC = 0.02  # real-like ~2% ground-based crowded-field photometric precision


def pspl_magnification(time: np.ndarray, t0: float, u0: float, tE: float) -> np.ndarray:
    u = np.sqrt(u0**2 + ((time - t0) / tE) ** 2)
    return (u**2 + 2) / (u * np.sqrt(u**2 + 4))


def planetary_perturbation(time: np.ndarray, tp: float, tE: float, q: float) -> np.ndarray:
    amp = AMP_COEFF * np.sqrt(q)
    sigma = DUR_COEFF * tE * np.sqrt(q)
    return amp * np.exp(-0.5 * ((time - tp) / sigma) ** 2)


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    time = np.arange(0, BASELINE_DAYS, CADENCE_HR / 24.0)
    a_single = pspl_magnification(time, TRUE_T0_DAYS, TRUE_U0, TRUE_TE_DAYS)
    delta = planetary_perturbation(time, TRUE_TP_DAYS, TRUE_TE_DAYS, TRUE_Q)
    a_true = a_single * (1 + delta)

    a_obs = a_true * (1 + rng.normal(0, PHOT_NOISE_FRAC, size=time.size))
    a_err = PHOT_NOISE_FRAC * a_true

    # Step 1: fit the single-lens (PSPL) model, masking out the perturbation window.
    pert_window = np.abs(time - TRUE_TP_DAYS) < 4 * DUR_COEFF * TRUE_TE_DAYS * np.sqrt(TRUE_Q) + 1.0
    fit_mask = ~pert_window
    popt, pcov = curve_fit(pspl_magnification, time[fit_mask], a_obs[fit_mask], p0=[TRUE_T0_DAYS, 0.15, 18.0], sigma=a_err[fit_mask])
    fit_t0, fit_u0, fit_tE = popt

    # Step 2: isolate the residual perturbation and fit its amplitude/duration.
    baseline_model = pspl_magnification(time, fit_t0, fit_u0, fit_tE)
    residual_frac = a_obs / baseline_model - 1.0

    def gaussian_bump(t, amp, tp, sigma):
        return amp * np.exp(-0.5 * ((t - tp) / sigma) ** 2)

    pert_fit_mask = np.abs(time - TRUE_TP_DAYS) < 3.0
    bpopt, bpcov = curve_fit(
        gaussian_bump, time[pert_fit_mask], residual_frac[pert_fit_mask],
        p0=[0.3, TRUE_TP_DAYS, 0.5], sigma=a_err[pert_fit_mask] / baseline_model[pert_fit_mask],
    )
    fit_amp, fit_tp, fit_sigma = bpopt
    fit_amp = abs(fit_amp)

    # Step 3: invert the real scaling relations to recover the mass ratio q,
    # from the perturbation's amplitude and, independently, its duration.
    q_from_amp = (fit_amp / AMP_COEFF) ** 2
    q_from_duration = (fit_sigma / (DUR_COEFF * fit_tE)) ** 2
    q_recovered = np.sqrt(q_from_amp * q_from_duration)  # geometric mean of two independent estimates

    tE_error_pct = abs(fit_tE - TRUE_TE_DAYS) / TRUE_TE_DAYS * 100
    q_error_pct = abs(q_recovered - TRUE_Q) / TRUE_Q * 100

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "injected", "recovered", "error_pct"])
        writer.writerow(["tE_days", TRUE_TE_DAYS, f"{fit_tE:.3f}", f"{tE_error_pct:.2f}"])
        writer.writerow(["u0", TRUE_U0, f"{fit_u0:.4f}", f"{abs(fit_u0-TRUE_U0)/TRUE_U0*100:.2f}"])
        writer.writerow(["mass_ratio_q_from_amplitude", TRUE_Q, f"{q_from_amp:.2e}", f"{abs(q_from_amp-TRUE_Q)/TRUE_Q*100:.2f}"])
        writer.writerow(["mass_ratio_q_from_duration", TRUE_Q, f"{q_from_duration:.2e}", f"{abs(q_from_duration-TRUE_Q)/TRUE_Q*100:.2f}"])
        writer.writerow(["mass_ratio_q_combined", TRUE_Q, f"{q_recovered:.2e}", f"{q_error_pct:.2f}"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    axes[0].errorbar(time, a_obs, yerr=a_err, fmt=".", ms=2, color="#9fb3a8", alpha=0.4, label="Simulated survey data")
    axes[0].plot(time, a_single, color="#a8431f", lw=1.2, ls="--", label="Single-lens (PSPL) model")
    axes[0].plot(time, a_true, color="#1f4e79", lw=1.0, label="True lens+planet model")
    axes[0].set_xlabel("Time [days]")
    axes[0].set_ylabel("Magnification A(t)")
    axes[0].set_title("Microlensing event with planetary perturbation")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.25)

    zoom_mask = np.abs(time - TRUE_TP_DAYS) < 4
    axes[1].errorbar(time[zoom_mask], (a_obs[zoom_mask] / baseline_model[zoom_mask] - 1) * 100, yerr=(a_err[zoom_mask] / baseline_model[zoom_mask]) * 100, fmt="o", ms=4, color="#1f4e79", capsize=2, label="Residual (data / single-lens - 1)")
    t_fine = np.linspace(time[zoom_mask].min(), time[zoom_mask].max(), 300)
    axes[1].plot(t_fine, gaussian_bump(t_fine, fit_amp, fit_tp, fit_sigma) * 100, color="#a8431f", lw=1.5, label="Fitted planetary perturbation")
    axes[1].set_xlabel("Time [days]")
    axes[1].set_ylabel("Fractional deviation [%]")
    axes[1].set_title(f"Planetary perturbation (q = {q_recovered:.2e})")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.25)

    fig.suptitle("Gravitational microlensing: recovering a planetary mass ratio from a perturbation")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "microlensing_recovery.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'microlensing_recovery.png'}")
    print(f"Injected tE {TRUE_TE_DAYS} d -> recovered {fit_tE:.3f} d ({tE_error_pct:.2f}% error)")
    print(f"Injected q {TRUE_Q:.2e} -> recovered {q_recovered:.2e} ({q_error_pct:.2f}% error)")


if __name__ == "__main__":
    main()
