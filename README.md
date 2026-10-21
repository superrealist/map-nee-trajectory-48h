# MAP–NEE joint trajectories over 48 hours

Core database extraction and statistical analysis code for a retrospective
study of joint mean arterial pressure (MAP) and norepinephrine-equivalent dose
(NEE) trajectories in adults with sepsis receiving vasopressors.

The analysis uses a **48-hour landmark** among patients alive and hospitalized
at that landmark. The primary model has three joint trajectory classes, with
serial correlation accounted for in the trajectory model. External application
uses the fitted development-model parameters.

## Repository contents

| Directory | Purpose |
|---|---|
| [`extraction/`](extraction/README.md) | Read-only PostgreSQL source queries for MIMIC-IV v3.1 and eICU-CRD v2.0 |
| [`core_code/01_cohort_and_features/`](core_code/01_cohort_and_features/README.md) | Eligibility, source-state NEE integration, bounded baseline measurements and SOFA proxies |
| [`core_code/02_trajectory/`](core_code/02_trajectory/README.md) | Joint trajectory likelihood, fitting and posterior classification |
| [`core_code/03_association/`](core_code/03_association/README.md) | Sequentially adjusted mortality associations, multiple imputation and classification uncertainty |
| [`core_code/04_cif/`](core_code/04_cif/README.md) | Cumulative incidence with live discharge as a competing event |
| [`core_code/05_conditional_information/`](core_code/05_conditional_information/README.md) | Conditional MAP/NEE information and pooled inference |
| [`core_code/06_incremental_value/`](core_code/06_incremental_value/README.md) | Risk models and incremental performance comparisons |
| [`docs/analysis_contract.md`](docs/analysis_contract.md) | Time definitions, inputs and reproducibility scope |

## Getting started

Python dependencies for the analysis are listed in `core_code/requirements.txt`;
the extraction driver has a separate `extraction/requirements.txt`. R packages
and input contracts are documented beside the relevant modules. These are core
functions: use the module README files to prepare their inputs and assemble the
analysis appropriate to an authorized dataset.

```bash
python -m pip install -r core_code/requirements.txt
python -m pip install -r extraction/requirements.txt
python extraction/extract.py --dataset mimiciv --list
python extraction/extract.py --dataset eicu --list
```

## Data access and reproducibility

**This repository contains code and documentation only.** It includes no
patient records, hospital extracts, credentials, private identifiers, fitted
patient-level outputs or unpublished manuscript files. Researchers must obtain
their own authorized access to [MIMIC-IV](https://physionet.org/content/mimiciv/)
and [eICU-CRD](https://physionet.org/content/eicu-crd/). The local hospital cohort
is not distributed, and its database extraction code is not included.

The analysis modules are the core-code release prepared with the manuscript.
The extraction directory provides a portable source-query layer. Cohort
adjudication, full feature engineering, initialization/preconditioning objects
and fitted model parameters are not included. Consequently, this is a core
methods release, not a complete executable reproduction archive. See the
[analysis code guide](core_code/README.md) for module-specific boundaries.

## Provenance and licensing

MIMIC derived concepts are external dependencies from
[MIT-LCP/mimic-code](https://github.com/MIT-LCP/mimic-code). Their licenses and
attribution requirements remain applicable to those external resources.
Database access agreements apply independently of this code repository.

No project-wide open-source license is assigned in this release. Public
availability alone does not grant additional reuse rights; contact the
repository owner for permissions beyond applicable law and GitHub's terms.

When referring to this code, cite the repository URL and the specific commit.
A publication citation can be added when bibliographic details are available.
