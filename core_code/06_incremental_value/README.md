# Module 06 — Frozen risk comparisons and incremental value (Figure 5)

`risk_numeric_core.R` contains the exact ridge fitting, weighted AUC, Brier and calibration helpers. `risk_models.R` provides reusable training, frozen application and paired evaluation interfaces. Dependencies: R, glmnet, splines (base recommended package).

```r
library(glmnet)
source('risk_numeric_core.R')
source('risk_models.R')
# Objects below are prepared by the authorized analyst; none are distributed.
fit <- fit_risk_models(training)
p <- predict_risk_models(fit, evaluation)
result <- evaluate_pairs(evaluation$hospital_mortality, p,
                         evaluation$patient_id,
                         hospital_id = NULL) # supply eICU hospital_id
```

## Input schema

One row per stay. Raw fields: age (years), sex (same development coding), icu (harmonized development levels), MAP (mmHg), HR (/min), lactate (mmol/L), NEE (µg/kg/min), joint_p2, joint_p3. SOFA fields are `S0_full_proxy`, `S1_full_proxy`, `S2_full_proxy`, and their `nonCV` counterparts. Training additionally supplies hospital_mortality (0/1). Missing continuous observations are NA; do not substitute zero for an unobserved score. Evaluation supplies patient_id and eICU hospital_id for inference.

S0 denotes the pre-T0 score; S1 and S2 denote 0–24 and 24–48 h scores. Dynamic SOFA is the joint inclusion of these three summaries, not a separately fitted latent SOFA trajectory. Their source algorithms and missing-domain conventions are documented in the supplement.

## Model specification

Five references, each fitted alone (`a`) and plus two frozen joint-class posteriors (`b`):

| Key | Reference predictors |
|---|---|
| T0 | Age, sex, ICU type, baseline MAP, HR, ln(1+lactate), ln(NEE+0.01), non-CV S0 SOFA |
| static_full | Full S0 SOFA |
| static_nonCV | Non-CV S0 SOFA |
| dynamic_full | Full S0, S1, S2 SOFA |
| dynamic_nonCV | Non-CV S0, S1, S2 SOFA |

The main figure displays six models: T0, static full SOFA, dynamic full SOFA, each ± joint posteriors. All ten models enter the simultaneous-inference family. Chronic diagnoses/CCI are adjustment covariates in the association analysis but are not retrospective predictors in this frozen risk analysis.

Each continuous predictor is completed with its MIMIC training median and adds a missing indicator. Natural cubic splines use two training observed-value quantiles (1/3, 2/3); tied/invalid interior knots trigger a linear term. SOFA boundary knots are 0–24 (full) or 0–20 (non-CV); other boundaries are the training observed range. Categorical levels and modal fallback for unseen levels are trained in MIMIC. Centering, scaling, zero-variance removal and sequential QR alias removal are training-only and frozen. Joint posteriors add two standardized columns.

The ridge objective is the **summed** negative logistic log-likelihood plus 0.5 × sum(non-intercept coefficients²), equivalent to glmnet `alpha=0, lambda=1/n` after standardization. The intercept is unpenalized. A Newton refinement verifies the score equation and reports effective df. No penalty tuning uses evaluation outcomes.

## Evaluation and uncertainty

- No risk fitting or recalibration in holdout, eICU or Local. Calibration regressions describe performance and are not applied to predictions.
- 1,000 common paired bootstrap draws per cohort, by patient; eICU uses crossed patient × hospital multiplicities. Local seed 2026100931, MIMIC/eICU production stream seed 2026100911. Isolated cohort calls need not reproduce archived random-number consumption exactly.
- Percentile 95% intervals for AUC, Brier, calibration-in-the-large and calibration slope. AUC explicitly handles tied scores.
- ΔAUC = enhanced − reference; ΔBrier = reference − enhanced. Positive values indicate improvement.
- Simultaneous intervals use the 95th percentile of the maximum absolute standardized bootstrap deviation over 5 pairs × 2 increment metrics (10 comparisons) within cohort.
- These are conditional evaluation intervals for fixed preprocessing, trajectory and risk models. They do not propagate risk-training/trajectory-estimation uncertainty and must not be described as the full-pipeline bootstrap.
- T0 measurements are event-time baseline features assessed at the 48 h landmark. Some source entries were accepted up to 48 h; this is not a model deployable prospectively at T0.

Sources: canonical `fig45_risk_models.R`, `local_31_apply_frozen_risk_models.R`; numerical helpers from `30_prediction_revision.R`. This package does not contain fitted coefficients, patient-level predicted risks, ROC thresholds or any input dataset.
