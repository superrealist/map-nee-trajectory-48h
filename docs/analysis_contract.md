# Analysis contract

## Time and eligibility

1. **T0** is the first source-supported vasopressor initiation in the index ICU
   stay. MIMIC uses qualifying continuous-infusion records; eICU requires the
   first positive record with a same-agent positive confirmation within 120
   minutes. It returns the first record's time, not the confirmation time.
2. **Landmark** is T0 + 48 hours. Patients must be alive and hospitalized then.
   Earlier ICU discharge does not itself exclude a patient; subsequent
   trajectory values remain structurally unobserved.
3. **Trajectory window** is [T0, T0 + 48 h), with hourly aggregation. Exposure
   minutes are truncated at ICU exit. Apply the observation-quality gate in
   `quality48()` before analysis.
4. **Baseline variables** follow their documented event/availability windows.
   Retrospective T0 vital-sign summaries and electronically available T0
   laboratory results are different contracts. Chronic history and weight
   require their source-specific preparation described in the study methods.
5. **Outcome** is hospital death after the landmark. The competing-risk
   analysis treats live discharge as a competing event. A 28-day CIF display
   horizon is distinct from the full-hospital-stay binary outcome.

## Preprocessing

- MAP is measured in mmHg. NEE is expressed in micrograms/kg/min-equivalent.
- The five drug components are norepinephrine, epinephrine, vasopressin,
  dopamine and phenylephrine. Conversion depends on the source unit and, where
  required, weight and explicitly documented concentration.
- A documented zero, an expired dose update and an unquantifiable positive
  infusion are separate states. State assumptions and the bounded-zero
  alternative are in `NEE_source_state.py`.
- The trajectory transform is natural log(NEE + 0.01). Other analysis modules
  have their own documented transforms; do not substitute them interchangeably.
- SOFA proxy sums retain support counts and default-normal counts. Missing
  organ-domain evidence must remain identifiable even when a proxy convention
  assigns a zero contribution.

## Modeling and external evaluation

Patient-level training/holdout allocation must be retained across all stays
for the same individual. Fit development models using training data; apply
their frozen parameters and preprocessing to holdout and external records.
Do not refit an external trajectory model and label it frozen-model validation.

The trajectory likelihood accepts prepared arrays, missingness masks,
initialization values and numerical preconditioning inputs. These must be
constructed by the analyst; this release does not distribute fitted model
objects. Mortality models use their specified M0–M3 adjustment sets. Multiple
imputation, pooled covariance, cluster resampling and classification-error
calculations are documented in their module README files.

## Scope of the public extraction layer

SQL queries export raw source records and reference derived tables. They do
not encode every hospital-specific data-cleaning decision, infection
adjudication rule, variable mapping or full fluid-balance reconciliation.
Raw exports must therefore be transformed and audited before use as model
inputs. In particular, directly summing all fluid rows, carrying forward all
drug records or treating each SOFA source row as a T0 score would not follow
the study's analysis contract.

The repository deliberately contains no patient-level demonstration dataset.
Module checks should use synthetic inputs or authorized data in private storage.
