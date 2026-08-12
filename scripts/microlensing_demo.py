"""Gravitational microlensing method demonstration: simulate a single-
lens point-source magnification event with a short planetary
perturbation superimposed, then find and characterize that perturbation
without ever referencing the injected parameters that generated it --
only the light curve itself is used for detection and fitting.

This is a PEDAGOGICAL DEMONSTRATION with simulated data, not a specific
real target's raw archival light curve (see README.md for why, and see
this portfolio's *-exoplanet-report repos for 11 planets analyzed
directly from real archival JWST/HST/Spitzer/ground-based data). Full
binary-lens light curves require solving for image positions via
inverse ray-shooting or complex polynomial root-finding; this repo
instead uses the well-established analytic scaling relations for a
planetary perturbation's characteristic amplitude and duration relative
to the mass ratio q (Gould & Loeb 1992; Mao & Paczynski 1991), which is
a standard back-of-the-envelope tool for estimating a detected planet's
mass ratio and is explicitly labeled as an approximation, not a full
caustic computation.

Detection pipeline (the injected TRUE_* values below are used only to
generate the data, and are never passed into any of the fitting or
window-finding steps):
  1. Fit the single-lens model to the full light curve with data-driven
     starting guesses (peak time, an approximate impact parameter from
     the peak magnification, an approximate timescale from the light
     curve width).
  2. Flag points that deviate from that fit by more than a noise
     threshold and locate the longest contiguous run of flagged points
     -- this is the detected perturbation window.
  3. Re-fit the single-lens model with that window masked out, to get a
     baseline unaffected by the planetary signal.
  4. Fit a Gaussian bump to the residual inside (and just around) the
     detected window, seeded from the window's own location and width.
  5. Invert the perturbation's fitted amplitude and duration through
     the scaling relations to recover the mass ratio q.
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

# Injected "ground truth" event parameters (realistic Galactic-bulge-survey
# regime). Used ONLY to build the simulated light curve below -- none of
# the detection or fitting code references these.
TRUE_T0_DAYS = 25.0
TRUE_U0 = 0.12
TRUE_TE_DAYS = 20.0
TRUE_Q = 1.5e-3  # planet/star mass ratio, roughly a Jupiter-mass planet around a 0.5 Msun star
TRUE_TP_DAYS = 22.0  # time of closest planet-source approach

AMP_COEFF = 8.0  # order-unity scaling coefficient, Gould & Loeb (1992)-style
DUR_COEFF = 0.5

BASELINE_DAYS = 60.0
CADENCE_HR = 0.5  # high-cadence survey coverage (e.g. KMTNet-class)
PHOT_NOISE_FRAC = 0.02  # ~2% ground-based crowded-field photometric precision


def pspl_magnification(time: np.ndarray, t0: float, u0: float, tE: float) -> np.ndarray:
    u = np.sqrt(u0**2 + ((time - t0) / tE) ** 2)
    return (u**2 + 2) / (u * np.sqrt(u**2 + 4))


def planetary_perturbation(time: np.ndarray, tp: float, tE: float, q: float) -> np.ndarray:
    amp = AMP_COEFF * np.sqrt(q)
    sigma = DUR_COEFF * tE * np.sqrt(q)
    return amp * np.exp(-0.5 * ((time - tp) / sigma) ** 2)


def data_driven_initial_guess(time: np.ndarray, flux: np.ndarray) -> tuple[float, float, float]:
    """Rough starting parameters read off the light curve itself, the way
    a real search would seed a fit: no injected values used."""
    t0_guess = float(time[np.argmax(flux)])
    peak = float(flux.max())
    u0_guess = max(1.0 / peak, 0.02)  # small-u approximation A(u) ~ 1/u
    above_half = time[flux > 1 + (peak - 1) / 2]
    tE_guess = max(float(above_half.max() - above_half.min()) / 2, 3.0) if above_half.size > 2 else 15.0
    return t0_guess, u0_guess, tE_guess


def find_perturbation_window(time: np.ndarray, residual: np.ndarray, resid_err: np.ndarray, sigma_thresh: float = 3.0, min_run: int = 3):
    """Flag points that deviate from the baseline fit by more than
    sigma_thresh times the local noise, and return the time span of the
    longest contiguous run of flagged points. No injected truth is used."""
    flagged = np.abs(residual) > sigma_thresh * resid_err
    best_start, best_len = None, 0
    run_start = None
    for i, flag in enumerate(np.append(flagged, False)):
        if flag and run_start is None:
            run_start = i
        elif not flag and run_start is not None:
            if i - run_start > best_len:
                best_start, best_len = run_start, i - run_start
            run_start = None
    if best_start is None or best_len < min_run:
        return None
    window = time[best_start:best_start + best_len]
    return float(window.min()), float(window.max())


def main() -> None:
    FIG_DIR.mkdir(exist_ok=True)

    time = np.arange(0, BASELINE_DAYS, CADENCE_HR / 24.0)
    a_single = pspl_magnification(time, TRUE_T0_DAYS, TRUE_U0, TRUE_TE_DAYS)
    delta = planetary_perturbation(time, TRUE_TP_DAYS, TRUE_TE_DAYS, TRUE_Q)
    a_true = a_single * (1 + delta)

    a_obs = a_true * (1 + rng.normal(0, PHOT_NOISE_FRAC, size=time.size))
    a_err = PHOT_NOISE_FRAC * a_true

    # Step 1: fit the single-lens model to the full light curve, seeded
    # only from what the data itself shows.
    p0 = data_driven_initial_guess(time, a_obs)
    popt0, _ = curve_fit(pspl_magnification, time, a_obs, p0=p0, sigma=a_err, maxfev=20000)
    trial_model = pspl_magnification(time, *popt0)
    trial_residual = a_obs / trial_model - 1.0
    trial_resid_err = a_err / trial_model

    # Step 2: locate the perturbation from the residuals of that fit --
    # no injected time or mass ratio is referenced.
    window = find_perturbation_window(time, trial_residual, trial_resid_err)
    if window is None:
        raise RuntimeError("No perturbation window detected above the noise threshold.")
    win_start, win_end = window
    pad = 0.5 * (win_end - win_start) + 0.5

    # Step 3: re-fit the single-lens model with the detected window masked
    # out, to get a baseline uncontaminated by the planetary signal.
    fit_mask = (time < win_start - pad) | (time > win_end + pad)
    popt, _ = curve_fit(pspl_magnification, time[fit_mask], a_obs[fit_mask], p0=popt0, sigma=a_err[fit_mask], maxfev=20000)
    fit_t0, fit_u0, fit_tE = popt
    baseline_model = pspl_magnification(time, fit_t0, fit_u0, fit_tE)
    residual_frac = a_obs / baseline_model - 1.0

    # Step 4: fit a Gaussian bump to the residual inside the detected
    # window, seeded from the window's own center and width.
    def gaussian_bump(t, amp, tp, sigma):
        return amp * np.exp(-0.5 * ((t - tp) / sigma) ** 2)

    pert_fit_mask = (time > win_start - pad) & (time < win_end + pad)
    p0_bump = [float(np.max(np.abs(residual_frac[pert_fit_mask]))), (win_start + win_end) / 2, max((win_end - win_start) / 2, 0.2)]
    bpopt, _ = curve_fit(
        gaussian_bump, time[pert_fit_mask], residual_frac[pert_fit_mask],
        p0=p0_bump, sigma=a_err[pert_fit_mask] / baseline_model[pert_fit_mask], maxfev=20000,
    )
    fit_amp, fit_tp, fit_sigma = bpopt
    fit_amp = abs(fit_amp)
    fit_sigma = abs(fit_sigma)

    # Step 5: invert the scaling relations for the mass ratio q, from the
    # perturbation's fitted amplitude and, independently, its duration.
    q_from_amp = (fit_amp / AMP_COEFF) ** 2
    q_from_duration = (fit_sigma / (DUR_COEFF * fit_tE)) ** 2
    q_recovered = np.sqrt(q_from_amp * q_from_duration)  # geometric mean of two independent estimates

    tE_error_pct = abs(fit_tE - TRUE_TE_DAYS) / TRUE_TE_DAYS * 100
    q_error_pct = abs(q_recovered - TRUE_Q) / TRUE_Q * 100
    tp_error_days = abs(fit_tp - TRUE_TP_DAYS)

    summary_path = FIG_DIR / "summary_statistics.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["quantity", "injected", "recovered", "error_pct"])
        writer.writerow(["detected_window_days", "-", f"[{win_start:.2f}, {win_end:.2f}]", "-"])
        writer.writerow(["tE_days", TRUE_TE_DAYS, f"{fit_tE:.3f}", f"{tE_error_pct:.2f}"])
        writer.writerow(["u0", TRUE_U0, f"{fit_u0:.4f}", f"{abs(fit_u0-TRUE_U0)/TRUE_U0*100:.2f}"])
        writer.writerow(["tp_days", TRUE_TP_DAYS, f"{fit_tp:.3f}", f"{tp_error_days:.3f} days abs. error"])
        writer.writerow(["mass_ratio_q_from_amplitude", TRUE_Q, f"{q_from_amp:.2e}", f"{abs(q_from_amp-TRUE_Q)/TRUE_Q*100:.2f}"])
        writer.writerow(["mass_ratio_q_from_duration", TRUE_Q, f"{q_from_duration:.2e}", f"{abs(q_from_duration-TRUE_Q)/TRUE_Q*100:.2f}"])
        writer.writerow(["mass_ratio_q_combined", TRUE_Q, f"{q_recovered:.2e}", f"{q_error_pct:.2f}"])

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    axes[0].errorbar(time, a_obs, yerr=a_err, fmt=".", ms=2, color="#9fb3a8", alpha=0.4, label="Simulated survey data")
    axes[0].plot(time, a_single, color="#a8431f", lw=1.2, ls="--", label="Injected single-lens model")
    axes[0].plot(time, a_true, color="#1f4e79", lw=1.0, label="Injected lens+planet model")
    axes[0].axvspan(win_start, win_end, color="#5cbf8a", alpha=0.15, label="Detected perturbation window")
    axes[0].set_xlabel("Time [days]")
    axes[0].set_ylabel("Magnification A(t)")
    axes[0].set_title("Microlensing event with planetary perturbation")
    axes[0].legend(fontsize=7)
    axes[0].grid(alpha=0.25)

    zoom_mask = np.abs(time - fit_tp) < 4
    axes[1].errorbar(time[zoom_mask], (a_obs[zoom_mask] / baseline_model[zoom_mask] - 1) * 100, yerr=(a_err[zoom_mask] / baseline_model[zoom_mask]) * 100, fmt="o", ms=4, color="#1f4e79", capsize=2, label="Residual against re-fit baseline")
    t_fine = np.linspace(time[zoom_mask].min(), time[zoom_mask].max(), 300)
    axes[1].plot(t_fine, gaussian_bump(t_fine, fit_amp, fit_tp, fit_sigma) * 100, color="#a8431f", lw=1.5, label="Fitted perturbation")
    axes[1].set_xlabel("Time [days]")
    axes[1].set_ylabel("Fractional deviation [%]")
    axes[1].set_title(f"Detected perturbation (q = {q_recovered:.2e})")
    axes[1].legend(fontsize=7)
    axes[1].grid(alpha=0.25)

    fig.suptitle("Gravitational microlensing: blind detection and mass-ratio recovery")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "microlensing_recovery.png", dpi=200)

    print(f"Wrote {summary_path}")
    print(f"Wrote {FIG_DIR / 'microlensing_recovery.png'}")
    print(f"Detected perturbation window: [{win_start:.2f}, {win_end:.2f}] days (injected tp was {TRUE_TP_DAYS} days)")
    print(f"Injected tE {TRUE_TE_DAYS} d -> recovered {fit_tE:.3f} d ({tE_error_pct:.2f}% error)")
    print(f"Injected q {TRUE_Q:.2e} -> recovered {q_recovered:.2e} ({q_error_pct:.2f}% error)")


if __name__ == "__main__":
    main()
