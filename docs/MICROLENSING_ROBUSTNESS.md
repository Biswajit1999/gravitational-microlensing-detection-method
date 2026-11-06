# Synthetic microlensing robustness diagnostics

Two independent numerical questions are addressed without modifying the existing demonstration.

1. **Blended PSPL Fisher information.** The normalized model is F(t)=1+f_s[A(t)-1]. Finite differences construct the Jacobian in (u0, tE, f_s), holding t0 fixed. The inverse Fisher matrix estimates local conditional uncertainty for independent, known Gaussian flux errors. Singular values and condition number expose parameter degeneracy. This is not a full posterior or a binary-lens inference. The pseudo-inverse is descriptive if the design is nearly singular.
2. **Cadence-phase completeness.** A uniformly distributed phase offset is scanned for a finite-duration top-hat perturbation window. The reported fraction is the proportion of phases containing at least N observations. It is *geometric sampling coverage*, not detection efficiency: no noise, weather, telescope gaps, selection effects, finite-source structure or caustics are included.

Reproduce: `pytest tests/test_microlensing_robustness.py -q`. This standalone module requires NumPy. The tests check PSPL symmetry, the unblended limit, Fisher scaling with noise, and exact cadence coverage cases. Synthetic calculations only; no observational planet claim. Prepared with AI coding assistance; Git author attribution follows the authenticated account and repository policy.
