# Module 03 — Mortality association and uncertainty

`association_core.R` contains the essential M0–M3 specification, fresh baseline multiple imputation, clustered covariance, scalar Rubin pooling, and the boot-then-impute variance calculation used for Table 2. `three_step_ml.py` contains the conditional three-step classification-error sensitivity. There is no database connection, dataset, fitted object, or batch scheduler in this release.

## Analysis frame and transformations

One row represents an eligible ICU stay at the 48 h landmark. Required columns are listed in `analysis_input`. Main outcome is binary hospital death after the landmark. `Jclass` is the current outcome-blind assigned class (1–3), with C1 as reference. `patient_id` and, for eICU, `hospital_id` remain outside the imputation matrix and identify sampling/covariance clusters.

The caller must supply the already constructed values with these exact scales:

| Column | Scale |
|---|---|
| `age_per10` | Age in years / 10 |
| `female` | Indicator for female sex |
| `charlson_retrospective` | Complete 17-domain, age-free Charlson index |
| `infection_site` | Recorded infection-site category |
| `noncv_sofa` | Pre-T0 SOFA excluding the cardiovascular domain, 0–20 |
| `baseline_lactate` | ln(1 + lactate in mmol/L) |
| `baseline_MAP` | Pre-T0 MAP in mmHg / 10 |
| `heart_rate10` | Pre-T0 heart rate in beats/min / 10 |
| `baseline_logNEE` | log2(initial NEE in microgram/kg/min + 0.01) |
| `icu_group` | Source-harmonized ICU type |

The base-2 baseline NEE transformation in association analysis differs from the natural-log trajectory transformation. This changes covariate units, not the trajectory definition.

M0 includes class alone; M1 adds age and sex; M2 adds Charlson and infection site; M3 adds non-CV SOFA, lactate, MAP, heart rate, initial NEE, and ICU type. The local association model omits the near-constant ICU-type term and retains the other terms.

Continuous observed 0–48 h marker summaries and the current posteriors are imputation auxiliaries in MIMIC/eICU; the local imputation uses the observed marker summaries. Mean-marker spline columns must already be prepared using the fixed development basis. Outcome, class, posteriors, and hourly trajectories are not imputed. Outcome use in association imputation is separate from outcome-blind trajectory fitting and from frozen prediction-model preprocessing.

## Reference analysis

```r
source("association_core.R")
x <- analysis_input(analysis_frame, local=FALSE)
imp <- fresh_imputation(x, m=40, seed=2026109000L)
result <- fit_nested(imp$mi, patient_id, hospital_id=NULL)
```

Use hospital identifiers for eICU. Numeric missing baseline covariates use predictive mean matching with five donors. Iterations begin at 20 and increase by 20 to at most 120; release requires all applicable final mean/SD PSRF values ≤1.1. If a resample has exactly one missing value for a covariate, its chain sample-SD diagnostic is inapplicable and its mean diagnostic is retained. Retain logged imputation events and actual effective predictor matrices for audit.

Reference/local intervals use 40 imputations, HC1 patient clustering (additional crossed hospital clustering in eICU), Rubin pooling, and finite cluster degrees of freedom. A negative-eigenvalue sandwich covariance triggers the explicitly recorded `sandwich::vcovCL(..., fix=TRUE)` projection, rather than an unreported correction.

## Full-process MIMIC/eICU bootstrap

1. Obtain replicate-specific trajectory assignments and integer multiplicities from Module 02.
2. Repeat rows by those multiplicities, keeping current class/posteriors paired with the correct stay.
3. Impute anew with `m=5`. Seed is `2026109000 + b` for MIMIC and `2026109000 + b + 1000` for eICU.
4. Call `fit_nested(..., outer_bootstrap=TRUE)` and save its individual log-ORs.
5. For each cohort/model/contrast, construct a **200 × 5** log-OR matrix `Q` in bootstrap/imputation order and call `boot_then_impute(Q)`.

This is boot-then-impute ANOVA pooling. It separates between-bootstrap from finite-imputation variation and uses a t interval. The individual model-based `within_U` values in outer replicates are diagnostic and are **not added again**. The final point estimate is the grand mean of bootstrap/imputation log-ORs; the original 40-imputation estimate is a separate reference. No percentile interval across the 1,000 individual values is used.

## Three-step sensitivity

Supply `X` with an intercept and completed covariates; binary `y`; zero-based assigned labels `w`; and the error matrix `D = error_matrix(posterior)`. The likelihood is

`sum_c P(C=c|X) D(W|C=c) P(Y|C=c,X)`.

The multinomial class gate and logistic outcome model are fitted jointly, while the measurement model and posterior-derived D remain fixed. For each imputation use `robust_fit`, then `profile_target_covariance` for the two class contrasts, and `pool` across imputations. These intervals condition on D and the measurement model. They do not replace the primary MIMIC/eICU full-process bootstrap.

R dependencies used in the study: R 4.5.2, mice 3.19.0, sandwich 3.1-1. The original orchestration additionally used data.table/jsonlite. Python sensitivity dependencies: NumPy and SciPy. No local trajectory bootstrap was performed; local intervals are conditional on the supplied frozen classifications.
