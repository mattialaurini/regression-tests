# Provenance and scope

The frozen numerical folders originate from:

- `nist_separable_benchmark_v1_0_6.zip`;
- `nist_tolerance_study.zip`, numerical recovery/instrumentation version 2.5.0;
- `nist_diagonal_grid_v250.zip`, matched diagonal-grid study.

The base MATLAB/Python files and tolerance-extension MATLAB files were
compared with the source texts embedded in the completed comparison MAT
file dated 2026-11-11. The diagonal package's numerical sources and its
frozen per-problem base options were retained. SHA-256 values for 73 frozen
MATLAB/Python files are in `reference/numerical_source_sha256.json`.
The JSON checksums in the original algorithm folders additionally document
their numerical-core snapshots.

`data/frozen_inputs.mat` is a portable MATLAB v5 file. Numerical arrays were
extracted from the saved tuning manifest, restoring MATLAB dimensions from
HDF5 storage order. It contains all eight problems, including their shared
feasible initial points, per-problem ablation base options, and the comparison
configuration. Machine-specific directories and MATLAB's serialized string
objects were removed. The loader restores string metadata from the unchanged
problem catalog and option strings explicitly. The saved observations and
boxes are not estimated from rounded manuscript tables.

The CSV files in `data/` expose original predictors and untransformed response
values for inspection. The frozen MAT additionally includes the transformed
Nelson response and normalization used by all solvers. NIST prescribes the
data/models/reference values; the finite parameter boxes are study choices.
The original catalog and supplied base package provide the data URLs.

The published comparison CSV combines the final 24 BnB runs (7/8, 8/8,
8/8 certified cases at the three tolerances) and the 24 previously completed
SCIP runs (1/8, 0/8, 0/8). The published ablation CSV contains 180 observations;
the progress file additionally includes 12 planned memory exclusions.
Source filenames are retained without machine-specific parent paths.

The new top-level portability drivers and `scripts/build_tables.py` are
packaging additions. They are not the original benchmark launch scripts.
The numerical core is unchanged. The ablation driver executes fresh trials
under the per-trial budgets recorded in the reference CSV by default, and
writes checkpoints using the original record/export helpers.
The alternative uniform-budget mode is explicitly a new experiment.

Published `wall_time`, native statuses, gap values, counters and actual
budgets are preserved in CSV. Displaying a configured timeout in a table
does not overwrite a recorded time. QP dimension counts are analytical
pre-presolve model sizes, not measurements of Gurobi's reduced models.

Original optional variable-projection routines remain in the frozen base
source tree for completeness, but no such runs are launched by this package.
Only BnB, SCIP and the matched-grid ablation are reproduction targets here.
