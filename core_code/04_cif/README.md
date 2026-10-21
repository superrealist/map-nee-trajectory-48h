# Module 04 — Hospital-death cumulative incidence (Figure 3)

`cif_core.py` contains the numerical estimator and cluster bootstrap. Import its functions; the module has no automatic file I/O. Dependencies: Python 3, NumPy.

## Required input

One aligned row per eligible ICU stay: `time` (days from T0 + 48 h to the validated hospital endpoint; strictly positive), `death` (1 hospital death, 0 live hospital discharge), `class_id` (1–3), `patient_id`; eICU also requires `hospital_id`. Identifiers are used only to define resampling clusters. The package contains no input records.

Call `class_curves(time, death, class_id, patient_id, hospital_id=None)` for MIMIC. Supply hospital IDs for eICU. For Local set `B=500, seed=20261009071`; MIMIC/eICU used B=200. Each class is evaluated using the same cohort-level resampling draws. The original MIMIC/eICU script consumed one shared RNG stream successively; isolated calls preserve the algorithm but not necessarily the archived random draws.

## Endpoint and uncertainty contract

- Time zero is the independent 48 h landmark; the graph shows 28 subsequent days.
- Hospital death is the event of interest; live hospital discharge is competing. Follow-up for the endpoint is the remaining index hospitalization, not a universal 28-day mortality endpoint.
- The supplied estimator assumes complete ascertainment of death or live discharge. `death=0` must not contain administrative loss to follow-up.
- Tied deaths/discharges use a common pre-event risk set. Terminal CIF equals the observed hospital-death proportion.
- The supplementary KM curve censors live discharge and estimates a different quantity from CIF.
- Intervals are pointwise, conditional on the frozen classification. They do not propagate trajectory refitting uncertainty. eICU uses independent patient and hospital multinomial multiplicities multiplied at stay level.

Source: canonical `code/cohortpub_prepare.py` (`cluster_weights`, `estimators`, `curves`) and `code/local_30_local_curves.py`. Timestamp validation occurs upstream and is deliberately not replaced by this small analysis module.
