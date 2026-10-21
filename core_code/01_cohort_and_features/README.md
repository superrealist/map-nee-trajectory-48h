# Cohort and feature core

This module contains reusable calculation kernels extracted or factorized from the locked 48-hour analysis. It is not the complete extraction pipeline and does not contain data or credentials.

- `cohort_rules.py`: official MIMIC Sepsis-3 event-time restriction; shared clinical eligibility contract; exact scope of local final-frame verification; paired 48-hour quality rule; eICU sustained-start T0.
- `NEE_source_state.py`: conversion and interval/point state reconstruction. Main eICU call uses `cap=120, stop_forever=True`. Returned arrays have historical allocation length 72; only the first 48 entries are used. Pass observation length no greater than 2880 minutes. Dose conversion expects explicit source units; unsupported units remain missing.
- `baseline_and_SOFA.py`: bounded pre-T0 laboratory selection, retrospectively ascertained initial-state vital signs, domain scoring, and explicitly labelled normal-default SOFA proxies.

## Units and time

Source-state functions use minutes relative to T0. Baseline selectors use hours relative to T0. `known_minutes` is the duration with quantifiable **total** NEE; a single unknown active drug makes total NEE unknown. Source NEE=0 requires an accepted stop-state rule. Clinical event time and electronic availability time are separate input fields.

T0 baseline NEE is the instantaneous source state, not the first hourly average. Baseline association uses `log2(NEE+0.01)`; trajectory measurement uses `ln(NEE+0.01)`.

## Local cohort

`common_eligibility()` is the common protocol contract. Its local flags require case-level source adjudication. `local_final_frame_verification()` documents the checks actually repeated on the supplied local final analysis frame; passing that function does not imply that the entire upstream selection has been independently repeated. No cohort is recalculated simply by importing this module.

## Integration

Prepare schema-matched source events privately, then call these functions as pure transformations. Unknown evidence is never substituted by a known positive flag. Figures and models use the already approved cohort memberships; changes in source eligibility require a separately documented rerun. Dependencies: Python, NumPy, pandas.
