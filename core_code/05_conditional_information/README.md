# Module 05 — Conditional component information (Figure 4)

This module asks whether MAP summaries contribute conditional information given NEE summaries, and conversely. It does not compare the AUC of joint and single-marker latent classes.

## Files and use

1. `summary_features.py`: observed hourly summaries, using natural log(NEE + 0.01). Matrix shape is stays × 48; missing hours remain missing.
2. `robust_MI_core.R`: audited multivariate MI D1-type Wald pooling and cluster sandwich covariance.
3. `conditional_core.R`: training-only spline specification, common full regression, two whole-block tests, Holm adjustment, optional scalar contrasts.

In R, source `robust_MI_core.R` then `conditional_core.R`. Build `basis_spec <- conditional_basis(mimic_training_features)` once, then reuse unchanged across all cohorts. Call `conditional_analysis(completed, features, basis_spec, patient_id, hospital_id)` for eICU and omit `hospital_id` for MIMIC. Use `local=TRUE` for Local. No file access occurs automatically.

## Required inputs and scales

`completed` is a list of 40 completed current-cohort association data frames, identically ordered. Obtain it from the shared baseline-MI analysis; include the outcome and component summaries as auxiliaries in that imputation model. This module does not start a separate MI. Row alignment must be checked by ICU-stay identifier before calling.

Each data frame contains hospital_mortality (0/1), age_per10, female (0/1), noncv_sofa, baseline_lactate=ln(1+lactate), baseline_MAP=MAP/10, heart_rate10=HR/10, baseline_logNEE=log2(NEE+0.01), charlson_retrospective (complete age-free CCI), infection_site, and icu_group. MIMIC/eICU additionally contain ICU_fraction, MAP_coverage_ICU and NEE_coverage_ICU. These are precomputed with the study's eligible ICU observation denominator, not whole-48-hour missingness fractions. Local omits ICU type (near-constant) and the non-estimable observation covariates.

`features` contains MAP_mean, MAP_slope24, MAP_SD, logNEE_mean, logNEE_slope24, logNEE_SD. Slopes are OLS change per 24 h; SD uses sample denominator n−1. MAP slope and SD enter per 10 mmHg. The marker location uses a three-column natural cubic spline, knots at observed training-cohort 1/3 and 2/3 quantiles and boundary knots at training range. Each primary component block has five terms (three location columns, slope, SD).

## Inference

- Fit one common full logistic model in each imputation; test either entire component block while retaining the other and all adjustment variables.
- Patient HC1 sandwich covariance for MIMIC/Local; crossed patient and hospital HC1 sandwich for eICU.
- Negative-eigenvalue corrections are explicitly audited. The D1 pooled block covariance must be positive definite; no silent generalized inverse.
- Complete-data df: unique patients minus 1; eICU unique hospitals minus 1. D1 uses the finite-df MI approximation implemented in the supplied numerical core.
- Holm correction covers the two directional primary tests within cohort. The optional NEE-positive-fraction variant is separate sensitivity analysis, available for MIMIC/eICU.
- This uncertainty is conditional on current observed summaries and the frozen training spline basis, using the shared 40 MI; it is not the 200-replicate full trajectory-refitting bootstrap.

Sources: canonical `fig45_conditional.R`, `local_21_local_association_prepare.py`, `local_22_local_association.R`; `00_inventory_and_features.py` and audited `robust_D1.R` from the conditional-SOFA revision.

`reported_component_contrasts(result)` returns the six displayed contrasts: MAP location 80 vs 70 mmHg, MAP slope +10 mmHg/day, MAP SD +5 mmHg, log-NEE location corresponding to 0.10 vs 0.05 µg/kg/min, log-NEE slope +ln(2)/day, and log-NEE SD +0.5. The log-NEE location is not an arithmetic dose mean; the slope contrast concerns NEE+0.01 and should not be described as literal doubling of NEE at every baseline dose.
