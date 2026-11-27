# Experiment-to-manuscript mapping

Section numbers can change; the descriptive manuscript titles below are used
to identify experiments without depending on a journal or section numbering.

| Manuscript result | Fresh-computation entry point | Post-processing / supplied reference |
|---|---|---|
| Patients' data | `pkpd/src/patient_profiles.m` | Thirteen fixed rows; no regeneration of profile 13 by averaging |
| `PatientsBIS.pdf` | `generate_pkpd_profiles` | Nominal PK, zero measurement noise, 150 samples, 1-s period |
| Approximate-likelihood ARX grid tuning and screening comparison | `run_pkpd_arx` | `export_pkpd_results`; cohort MAT must be computed or supplied |
| `tab_lin_tuning`, `tab_lin_results`, `lin_numbers` | `run_pkpd_arx` | `pkpd/src/paper_lin_outputs.m` |
| `lin_errors_order.pdf`, `lin_time_screening.pdf` | `run_pkpd_arx` | Same exporter |
| Exact-ML double-pole and one-compartment screening comparison | `run_pkpd_ml` | `export_pkpd_results`; cohort MAT must be computed or supplied |
| `tab_ml_results`, `ml_numbers` | `run_pkpd_ml` | `pkpd/src/paper_ml_outputs.m` |
| Patient 13: six-model reconstruction, sigma=1 | `run_pkpd_patient13` | `export_pkpd_reference` uses the supplied fitted predictions |
| `patient13_six_models_BIS.pdf`, `patient13_six_models_errors.pdf` | Same | Saved fits, signals and metrics match the manuscript's illustrated RMSE values |
| NIST global certification against SCIP | `run_comparison_bnb`, `run_comparison_scip` after `setup_experiments` | NIST reference CSVs and `nist/scripts/build_tables.py` |
| NIST matched-grid screening ablation | `run_ablation` | Same, with final gaps and cheap/costly node partitions |

## Preserved settings and interpretation

The original final-batch ARX and ML scripts default to four noise levels,
`[0,0.01,0.1,1]`, although one manuscript sentence still describes two.
Those defaults are preserved, yielding 416 ARX main runs and 208 ML runs.
The illustrated six-model case uses only sigma=1, as in its caption.
No extrapolated or invented cohort results have been substituted for the
missing full-cohort MAT files.

ARX minimizes the approximate criterion after clipping above E0; its linear
domain is expanded according to the original boundary test. Exact ML uses
unclipped data and estimates E0. Both use norm-gap tolerances, whereas NIST
uses normalized-RSS tolerances. Queue iteration counts in PK/PD are kept
distinct from NIST's unique visited regions.

Adaptive tuning, capping, domain expansion and local searches can make
selected grids, search trajectories and timings machine-dependent. The
recorded algorithmic policies are preserved; exact replication of every
timing or fitted parameter is not guaranteed by packaging alone.

## Latest manuscript clarifications

The profile table is identified in the manuscript as Padula et al. (2015),
Table 1, with the literature patient set traced to Ionescu et al. (2008).
Demographics determine the nominal Schnider PK model; Ce50, gamma, E0 and
Emax are tabulated Hill-map parameters kept fixed during response generation.
Input rates are in mg/s; volumes and clearances are in L and L/min, with
rate constants converted to 1/s. Profile 13 remains exactly as tabulated.

The exact-ML formulations jointly estimate the baseline and maximum effect
as two linear parameters. This retains the nonlinear dimension used by the
existing experiments. No known-baseline ML experiment has been added.

The ARX error-plot exporter now relocates only exact zeros below the smallest
positive error, as stated in the current caption; positive errors, fitted
parameters, tables and statistics are not truncated or changed.
