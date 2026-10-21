# Module 02 — Joint trajectories and outcome-blind transport

This module contains the computational core used for Figure 2 and class assignment in subsequent analyses. Study-specific file loading, source extraction, monitoring, and fitted model objects are omitted. It is not a rerun of the submitted study.

## Inputs

- `MAP`: an `N × 48` array of reconciled hourly MAP in mmHg.
- `NEE`: an `N × 48` array in microgram/kg/min, with source-supported zeros and genuinely missing hours already distinguished.
- An initial parameter dictionary: `beta` (2 × 3 × 4), `var` (2 × 2), `pi` (3), and `rho` (2). The second random-intercept variance is zero.
- A 31 × 31 preconditioner: the Cholesky factor of the inverse per-stay observed information from the initialization reference.
- Patient identifiers, fixed patient-level training/holdout assignments, and eICU hospital identifiers for clustered resampling. Identifiers are used only as grouping variables.

Dependencies: Python, NumPy, SciPy, and PyTorch. The study fit used float64 arithmetic and PyTorch 2.11.0 with CUDA; CPU application uses the same formulas.

## Main calls

```python
from trajectory_model import transform_hourly, tensors, accepted_model, posterior

Z = transform_hourly(MAP, NEE)
model, H = accepted_model(
    tensors(Z[training]), weights[training], initial, reference,
    preconditioner, precise=True,
)
assert model["stage1_gate_passed"]
probability, observed_loglik = posterior(Z, model)
class_id = probability.argmax(axis=1) + 1
```

`initial`, `reference`, `preconditioner`, and all input arrays must be supplied by the authorized analyst. No study data or fitted parameter file is included with this code-only release. In external application, use the finalized development parameters without fitting in the external cohort.

## Locked specification

- Three latent classes, two markers: MAP and `ln(NEE + 0.01)`.
- Hourly midpoints 0.5–47.5 h. A natural cubic cardinal basis has nodes 0, 16, 32, and 48 h; no extra intercept is added.
- Class-specific marker means; marker-specific AR(1) residuals shared across classes.
- MAP random intercept, integrated analytically; logNEE random-intercept variance fixed at zero.
- Conditional independence of the two markers given class.
- Missing observations are omitted from their marker likelihood; correlation is `rho**gap` for the actual elapsed hourly gap.
- Thirty-one estimated parameters: 24 mean coefficients, three variance parameters, two free class-probability logits, and two AR coefficients.
- MAP clipping bounds and marker centers/scales are inherited from earlier MIMIC training preprocessing and are fixed in the canonical cohort and all bootstrap/external applications. They are not re-estimated from external data.

The reference fit used eight starts: the earlier 48 h fit and seven perturbations of its means/variances. Starts 2–4 perturbed mean coefficients with SD 0.15; starts 5–8 used SD 0.4. Variances were perturbed multiplicatively with log-SD 0.1. The accepted candidate with the highest likelihood was selected. BFGS is followed by up to four Newton refinement steps. Acceptance requires maximum preconditioned gradient < 2e-5 and positive-definite per-stay information. AR coefficients are constrained to (1e-5, 0.995).

`align` uses the Hungarian algorithm to minimize mean squared distances between standardized two-marker mean curves; no outcome is used. `observed_class_quality` reports entropy, APP, and OCC. Back-transforming a fitted logNEE location gives a log-location, not an arithmetic mean NEE.

## Paired bootstrap interface

`paired_cluster_weights(..., bootstrap_index=b)` returns multiplicities. MIMIC patients are sampled separately within training and holdout strata. eICU receives the product of independently sampled patient and hospital multiplicities. The exact same integers must be used for the trajectory likelihood and subsequent association resample.

For each of 200 replicates, fit on weighted MIMIC training stays, align to the reference, and assign all MIMIC/eICU stays. The bootstrap retries a failed gate at most three times using a precise fit (including small mean perturbations on later retries). Stop rather than release an unsuccessful replicate. Pass current classes/posteriors and expanded multiplicities to Module 03. The original pipeline also saves hashes and diagnostic logs; the scheduler itself is outside this core release.

K, covariate design, preprocessing, eligibility, and source measurement rules are fixed in this bootstrap. It therefore does not include uncertainty in selecting these design choices.
