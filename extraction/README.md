# Core database extraction

These PostgreSQL queries extract the source records used to build MAP–NEE
features from **MIMIC-IV v3.1** and **eICU-CRD v2.0**. They retain original
values, units, identifiers, event times and available recording times. They
are a portable raw-source layer, not a replacement for the study's full
screening and adjudication workflow. No hospital-specific extraction adapter
is included.

## Prerequisites and execution

Obtain authorized access to each database. Load its PostgreSQL tables using
the respective database documentation. MIMIC queries additionally require
`mimiciv_derived.age`, `sepsis3`, `charlson`, `sofa` and `ventilation` from the
MIMIC Code repository. These derived SQL definitions are external dependencies
and are not copied into this repository.

Install `extraction/requirements.txt`. Configure the connection locally with
standard libpq settings (`PGSERVICE`, or `PGHOST`, `PGPORT`, `PGDATABASE`,
`PGUSER` and a password file). Use a read-only database account. Do not commit
connection settings, password files, exported records or identifiers.

```bash
python extraction/extract.py --dataset mimiciv --list
python extraction/extract.py --dataset mimiciv --query cohort_candidates --output data/mimiciv_candidates.csv
python extraction/extract.py --dataset mimiciv --query vasopressors_raw --stay-ids data/mimiciv_candidates.csv --output data/mimiciv_vasopressors.csv
python extraction/extract.py --dataset eicu --query cohort_candidates --output data/eicu_candidates.csv
python extraction/extract.py --dataset eicu --query vasopressors_raw --stay-ids data/eicu_candidates.csv --output data/eicu_vasopressors.csv
```

Only `cohort_candidates` scans all adult ICU stays. Every other query requires
a private CSV containing an `icu_stay_id` column. Outputs are local CSV files;
an existing output is not overwritten. Large extracts can be divided into
disjoint stay-ID batches. The SQL files contain named blocks, so execute them
through the driver rather than sending the whole file to PostgreSQL.

## Query map

| Source | Content |
|---|---|
| MIMIC `cohort_candidates`; eICU `cohort_candidates` | Age, sex, admission/ICU identifiers, stay boundaries and hospital outcomes |
| `vasopressors_raw` | All candidate drug records, including zeros and invalid or unconvertible rates for downstream source-state handling |
| MIMIC `chart_events`; eICU `vital_periodic` | MAP and other physiological measurements; MIMIC includes weight and selected respiratory/GCS items |
| `laboratory_events` | Raw baseline and organ-function laboratory measurements |
| Fluid input/output blocks | Original amounts, units and timestamps; cumulative and component records require reconciliation |
| Diagnoses/history blocks | Recorded diagnoses and chronic disease history |
| MIMIC derived blocks | Official Sepsis-3, Charlson, rolling SOFA and ventilation source tables |
| eICU medication/microbiology/respiratory/treatment blocks | Source evidence for infection and organ-support preparation |

## Downstream temporal and source rules

- MIMIC uses chart/infusion timestamps; eICU uses **minutes from ICU admission**.
- Candidate extraction is broader than final eligibility. Use source-supported
  T0, adult ICU eligibility, survival and hospitalization at T0 + 48 h, source
  quality and valid hospital outcome information before fitting models.
- For MIMIC, `official_mimic_by_landmark()` requires an official positive
  Sepsis-3 record with antibiotic, culture, suspected-infection and SOFA event
  times all present and no later than T0 + 48 h. eICU uses its documented
  infection/organ-dysfunction adjudication, not a claim of official MIMIC Sepsis-3.
- Restrict trajectory events to **[T0, T0 + 48 h)** and the observed index ICU
  stay. Do not fill the post-ICU portion of the window with normal values.
- Apply baseline event and availability-time rules in
  `core_code/01_cohort_and_features/baseline_and_SOFA.py`. Raw exports can contain
  later records; their presence does not make them baseline measurements.
- Drug records are raw evidence. Resolve duplicate/conflicting updates, reject
  invalid orders, standardize units and weight, and apply the published source
  state rules before calculating hourly NEE. Do not treat absent records as a
  documented zero. `NEE_source_state.py` contains the core conversions and state
  integration; explicit bag-concentration preprocessing for mL/h rates is an
  upstream requirement, not implemented by its `dose()` function.
- MIMIC rolling SOFA intervals are source data, not automatically valid T0
  values. Respect each interval's endpoint and the timing of its inputs.
- MIMIC laboratory records may have a missing admission identifier. The query
  retains the original linkage type; resolve duplicates and validate specimen
  type, units and linkage before deriving organ-domain scores.
- eICU laboratory names and infusion documentation vary across hospitals.
  The explicit field/name lists must be checked against each authorized local
  installation. Unknown values remain unknown during source preparation.

The extraction queries were reorganized for this public code release. They
do not include patient lists, investigator adjudication files, preprocessing
artifacts or a one-command reconstruction of the published analysis cohort.
