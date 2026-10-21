# Core analysis code for the 48 hour MAP and NEE study

This computational supplement contains the central algorithms supporting the final three-cohort study. It is organized by result, not by the historical project directory. It contains no participant records, hourly datasets, fitted model objects, credentials or database export pipeline.

## Modules

| Module | Result supported | Main functions |
| --- | --- | --- |
| 01_cohort_and_features | Figure 1 and Table 1 definitions | Eligibility flags, source-supported T0, NEE states, baseline windows and SOFA proxy aggregation |
| 02_trajectory | Figure 2 and class assignment | Observed-data joint Gaussian mixture likelihood, natural spline means, AR covariance, fitting, posteriors and label matching |
| 03_association | Table 2 | Nested M0–M3 logistic models, MI, cluster covariance, boot-then-impute pooling and conditional three-step sensitivity |
| 04_cif | Figure 3 and Figure S1 | Aalen–Johansen hospital-death CIF, discharge-censored KM and fixed-class cluster bootstrap |
| 05_conditional_information | Figure 4 | MAP and logNEE dynamic summaries, common-model conditional block tests and MI pooling |
| 06_incremental_value | Figure 5 and performance tables | Training-only preprocessing, frozen ridge risk models, paired AUC/Brier increments and calibration |

## Running the core

Read each module's README for input columns and scales. Functions consume authorized analysis arrays or tables supplied by the analyst. No example rows are represented as real study observations. An initial trajectory parameter object and information preconditioner are required for the original fitting strategy; neither is included in this code-only release. External application requires the finalized development parameters. The package deliberately does not initialize or locate private study data automatically.

1. Establish source eligibility and hourly measurement states before calling a model.
2. Keep the fixed patient-level MIMIC training/holdout split; do not split hourly records independently.
3. Fit or apply Module 02 using MAP and NEE only. Preserve class posterior probabilities.
4. For a full-process association bootstrap, pass the same original-unit sampling multiplicities and newly assigned classes from Module 02 to Module 03; use five fresh imputations per resample.
5. Use fixed-model predictions, separate from association MI, for Module 06 evaluation. Evaluation outcomes must not enter training preprocessing or external completion.

The trajectory marker uses natural log(NEE + 0.01). The initial-NEE covariate in the association design uses log base 2; the risk preprocessor has its own documented natural-log transform. Unit and transform contracts are intentional.

## Software

Python numerical dependencies are NumPy, SciPy, pandas and PyTorch. The trajectory implementation uses double precision. R dependencies are mice, sandwich, glmnet and their normal dependencies; R 4.5.2, mice 3.19.0 and sandwich 3.1-1 were used in the analysis. Use a compatible, isolated numerical environment. Database extraction and report-generation dependencies are outside this core release.

## What is reproducible from this release

The formulas and central computational procedures can be inspected and run with conforming authorized inputs. This is not a complete end-to-end extraction archive. Source abstraction, vocabulary mapping, source adjudication and the full batch scheduler remain outside the supplied core. The numerical algorithms were checked using synthetic inputs and, for frozen risk application, against the archived calculations without changing study estimates. Reorganized wrappers may consume random-number streams differently from the historical driver; exact archived draws are not promised by merely reusing a seed.

## Data access

No data are attached. Additional verification materials can be considered on reviewer request through the study team and the applicable source access arrangements. MIMIC-IV and eICU datasets are accessed through PhysioNet. A reviewer request does not itself replace database authorization or local approvals.
