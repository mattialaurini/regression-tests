# Validation scope

## Checks performed during package preparation

- Retained NIST numerical files match their frozen SHA-256 values: 73 files.
- The 48 reference BnB/SCIP comparison records retain their success counts,
  objective-gap classifications, timing data and counters.
- The 180 NIST ablation records and 12 memory exclusions retain their gaps,
  grids, budgets, and cheap/costly node partitions.
- All 67 mapped original PK/PD MATLAB files were compared after the documented
  English identifier/path substitutions, removal of comments and the separate
  plot-export extraction. Their optimization and model code is otherwise unchanged; the later plotting-only
  update is described in `LATEST_ALIGNMENT.md`.
- PK/PD function declarations match their distributed filenames.
- The plotting-only update preserves every strictly positive error and all
  source statistics; a MATLAB unit test for this behavior is supplied. It has
  not been executed in the packaging environment.
- Six supplied patient-13 reconstructions were recomputed independently in
  Python/SciPy from the stored ARX coefficients/bases and exact-ML pole/gain
  parameters. Maximum output discrepancy: 2.85e-14 BIS units.
- The stored patient-13 concentration was independently simulated from the
  saved continuous PK parameters using a matrix exponential and 1-s ZOH.
  Maximum discrepancy: 9.8e-15 concentration units.
- MATLAB column-vector dimensions were checked in the portable patient-13
  MAT file, including the ARX coefficient vectors and all sampled outputs.
- Three NIST table-export tests and three patient-13 CSV/metric tests pass
  without MATLAB or solver licenses.
- English documentation and retained runtime code were checked for leftover
  Italian instructions, obsolete journal labels and machine-specific paths.

## Previous local NIST smoke tests

The NIST package's preceding version was exercised by the user in MATLAB:
portable inputs, instrumentation tests, five real-Gurobi QP recovery cases,
one BnB BoxBOD run, one SCIP BoxBOD run, and eight BoxBOD ablation runs passed.
The retained NIST numerical sources are unchanged; this distribution removes
unused files and places the suite under `nist/`.

## Checks still requiring local execution

The packaging environment has no MATLAB, Control System Toolbox, or licensed
Gurobi installation. The new unified MATLAB launchers, figure exporters and
the complete PK/PD numerical-check scripts have therefore not been executed
here. Run the commands below locally before publishing the unified release:

```matlab
setup_experiments;
export_pkpd_reference;
run_pkpd_checks;
test_portable_inputs;
```

Then use small, separately saved fit subsets to validate fresh PK/PD execution
before launching a full cohort. The supplied source's final-batch notes are
provenance, not a substitute for executing the unified launchers locally.

No complete original ARX/ML cohort MAT files were supplied. No full-cohort
statistics or figures have been reconstructed from prose or rounded tables.
No new full benchmark battery has been run during packaging.

Numerical certificates use the implementations' floating-point solver
safeguards. These tests are consistency and regression checks, not formal
directed-rounding proofs or guarantees of identical timings on other hardware.
