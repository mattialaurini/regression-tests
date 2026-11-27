# Cleanup and source changes

The PK/PD directory is a dependency-based selection from the supplied archive.
It retains the two cohort scripts, the six-model illustration, their exporters,
the shared simulation/identification/bound routines, and numerical regression
checks relevant to the retained models. Earlier Wiener-model experiments,
real-data/Simulink scripts, exploratory benchmarks, unused plotting scripts,
`uniform_patient.mat`, duplicate `.fig`/PDF outputs, desktop metadata and the
development diary were excluded.

English replacements were made consistently at declarations and call sites:

| Original name | Distributed name |
|---|---|
| `fd_identifica` | `identify_arx` |
| `ml_identifica` | `identify_ml` |
| `sim_Propofol_optactiva` | `propofol_patient_parameters` |
| `dataPatients` | `patient_profiles` |
| `salva_figura` | `save_figure_pdf` |
| Output field `tempo_totale` | `total_time` |
| Output field `numero_nodi` | `iterations` |
| Output directory `risultati` | `results` |

`pkpd/source_mapping.json` records original/distributed filenames and checksums.
Non-numerical comments were updated, including explanatory notes on residual
scaling, seed order, screening validity and hybrid QP decisions. The private
patient-13 plotting routine was extracted to `export_patient13_outputs` so
saved predictions can be exported without invoking any identification.
Algorithms, grids, tolerances, data transformations and optimization constants
in the supplied PK/PD routines were not retuned. A later plotting-only
correction relocates exact zeros without flooring small positive errors,
matching the updated manuscript caption; statistics are computed from the
original values.

The patient-13 reference MAT contains only the data, model parameters,
predictions, metrics and selected numerical statistics needed to reconstruct
the illustration. Serialized function handles and MATLAB function-workspace
objects were removed. It is an export reference, not a fit-resume checkpoint.

The NIST numerical functions retained in this distribution are unchanged.
Obsolete examples, development-only functions, superseded standalone
launchers, generic old manuscript templates, and regenerable prebuilt tables
were excluded. Shared regression fixtures and source snapshots required by
the tests remain. Retained numerical-source hashes are verified by
`nist/scripts/verify_package.py`.

The machine-readable inventory `cleanup_manifest.json` lists excluded PK/PD
MATLAB files and removed NIST paths. File omissions never substitute a
different method for the recorded experiment.
