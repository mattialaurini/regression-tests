# NIST numerical experiments: reproducibility package

This page describes the `nist/` subtree. Shell paths below are relative to
that directory. In the unified repository, run `setup_experiments` from
the top level before invoking the MATLAB entry points.

Code and portable inputs for the NIST experiments in **An interpolation-based
branch-and-bound method for separable nonlinear regression: an application to
PK/PD model identification**.

This repository covers the comparison against SCIP and the matched-grid
ablation. It does not include the PK/PD experiments.

## What is included

- The numerical algorithms used for the experiments, with their QP recovery,
  branching, stopping and pruning instrumentation preserved.
- All eight frozen problem instances: observations, study parameter boxes,
  normalization and the common feasible initial points used in the experiments.
- Separate entry points for BnB and SCIP; no file-selection dialogs.
- Matched-grid ablation with and without the entire screening component.
- Published results as CSV files and a script to regenerate LaTeX tables.
- Checkpoints, partial reports, environment metadata and regression tests.

The numerical files are checked against saved experiment source snapshots.
New top-level drivers handle portability, published per-trial budgets and
reporting; they do not introduce lower bounds or change the numerical core.

## Dependencies

The experiment environment was:

| Component | Experiment version |
|---|---|
| CPU / RAM | Apple M2 Pro / 16 GiB |
| OS | macOS 26.6.2 |
| MATLAB | R2026a Update 3 |
| Optimization Toolbox | 26.1 |
| Gurobi | 12.0.2, MATLAB interface |
| Python | 3.12.15, out-of-process MATLAB integration |
| NumPy | 2.5.3 |
| PySCIPOpt | 6.2.1 |
| Native SCIP | 10.0.2 |

MATLAB, Optimization Toolbox and Gurobi must be installed and licensed
separately. This repository distributes no solver binaries or licenses.
The BnB experiments require `fmincon` and Gurobi. SCIP also requires Gurobi
because the shared initialization and preflight use it; the saved initial
points are loaded rather than recomputed during timed runs.

Install the SCIP Python dependencies in an interpreter supported by your
MATLAB release:

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-scip.txt
```

On Windows use `.venv\Scripts\python.exe`. Verify the **native SCIP** version
recorded in the resulting session metadata; a PySCIPOpt version alone does
not identify every possible native-library installation.

## Start here

Extract or clone the repository. Open MATLAB and change its current folder
to the repository root. Add your installed Gurobi MATLAB interface to the
MATLAB path, then execute:

```matlab
setup_reproduction;
test_portable_inputs;
run_comparison_bnb(DryRun=true);
run_ablation(DryRun=true);
```

These dry runs do not invoke Gurobi, SCIP or local optimization. They list
the selected experiments and memory exclusions. Inputs are fully bundled;
no NIST download or previous benchmark MAT file is required.

For SCIP, configure the interpreter before its first use:

```matlab
setup_reproduction(fullfile(pwd,'.venv','bin','python'));
```

On Windows substitute `fullfile(pwd,'.venv','Scripts','python.exe')`.
If MATLAB already loaded a different interpreter, restart MATLAB first.

Run the regression tests with licensed Gurobi available:

```matlab
run_regression_tests;
```

## 1. Comparison with SCIP

Run the two methods independently:

```matlab
[resultsBnB, summaryBnB] = run_comparison_bnb;
```

Then, after configuring Python:

```matlab
[resultsSCIP, summarySCIP] = run_comparison_scip;
export_comparison;
```

Each method runs 8 problems at each of the following settings:

| Relative tolerance | Nominal time limit |
|---|---|
| 0.1 | 100 s |
| 0.01 | 600 s |
| 0.001 | 3600 s |

The original BnB comparison uses initial `M=5`, screening `N=2`, and fine
grid cap `Nhat=3`, including its original grid coarsening/reinsertion policy.
It is distinct from the later ablation, where `Nhat=M`.
SCIP has no such grid parameters. It uses the same parameter boxes and
saved feasible initial point, one thread, feasibility tolerance `1e-6`,
fixed seeds and disabled presolve, as in the reported experiment.

Outputs: `outputs/bnb/`, `outputs/scip/`; merged comparison:
`outputs/comparison/`. Each includes CSV/MAT records. The legacy reporting
routine also writes LaTeX tables; the paper-table script below applies the
paper's convention for displayed timeout times explicitly.

## 2. Matched-grid ablation

```matlab
[resultsAblation, summaryAblation] = run_ablation;
```

All 8 problems, all 3 relative tolerances, and both variants are selected:

- `with_m2`: the full screening algorithm;
- `without_m2`: the frozen solver with the entire `m2` component removed.

The grids are `M in {3,5,9,17}`, `Nhat=M`, and `N=2` when screening is active.
Local searches remain enabled in both variants. Removing screening also
removes its queue-priority information and its poor-upper-bound local-search
trigger; the variants can consequently visit different regions. Other
local-search rules and budgets are retained. `N` is unused without `m2`.

There are 192 planned cases. The fixed 8-GiB planning guard excludes 12
Gauss1 cases (`M=9,17`, both variants and all tolerances) before allocation.
The memory estimate is **24 dense residual-by-grid arrays**, not a measured
peak RAM value. These exclusions are not timeouts or successful solves.

### Budget options

The nominal ablation schedule is 100/1000/3600 s. By default,
`BudgetMode="published"` starts fresh computations using the per-trial
limits recorded in `reference/ablation_results.csv`. The `time_limit`
column records the actual configured limit; `target_time_limit` records
its nominal schedule value. Both fields are preserved in new outputs.

To use the nominal schedule uniformly for a new experiment, select a
separate output directory:

```matlab
run_ablation(BudgetMode="uniform", ...
    OutputDirectory=fullfile(pwd,'outputs','ablation_uniform'));
```

Neither mode imports published timing results into new computations.

## Select, resume or restart

`CaseIndices=1`, `2`, `3` select relative tolerances 0.1, 0.01, 0.001.
For example, start with BoxBOD at tolerance 0.1:

```matlab
run_comparison_bnb(CaseIndices=1,ProblemIds="BoxBOD");
run_comparison_scip(CaseIndices=1,ProblemIds="BoxBOD");
run_ablation(CaseIndices=1,ProblemIds="BoxBOD");
```

After completing that subset, run the full command to add missing trials.
For the legacy comparison runner, changing a subset changes its frozen
manifest: use a separate `OutputDirectory` for a trial subset. The ablation
driver permits extending a subset in the same directory.

Repeat exactly the same command after an interruption to resume completed
checkpoints. Results are saved after each trial. A source/configuration
change requires a new output directory. To restart all computations without
deleting anything:

```matlab
run_comparison_bnb(OutputDirectory=fullfile(pwd,'outputs','bnb_fresh'));
run_ablation(OutputDirectory=fullfile(pwd,'outputs','ablation_fresh'));
```

An error is saved as an error and stops execution; it is never converted
to a timeout. Resolve the cause and inspect/remove only that failed
checkpoint before resuming. Keep the computer connected to power and its
power mode consistent; wall-time-dependent local searches can otherwise
change the search trajectory.

## Results and interpretation

- All bounds and global tolerances refer to `RSS/s`, with
  `s=max(NIST certified RSS,1e-8)`. Nelson uses logarithmic responses.
- The common success test is `L >= (1-epsilon_r)*U - epsilon_a`, with
  `epsilon_a=1e-10`. Raw RSS is also retained, using the recorded scale.
- `nodes_visited` counts distinct BnB regions, not repeated queue extractions.
  SCIP nodes are its processed-node count; these measure different work.
- For the ablation, **costly** means a visited region with at least one
  `m1` evaluation; **cheap** means another visited region. Their sum equals
  `nodes_visited`. Cheap does not mean no QP: screening and feasible-solution
  updates can still call Gurobi.
- `pruned_tolerance_only` is `Ntol`: discarded regions satisfying
  `(1-epsilon_r)*U-epsilon_a <= L_node < U`, at discard time. It includes
  eligible queued regions closed at successful global termination, and
  excludes regions left open at a time limit. No analogous SCIP count is
  claimed.
- Final absolute/relative objective gaps and actual/target budgets are kept
  separately. A timeout's raw `wall_time` is **not overwritten**. Tables mark
  timeouts and display the configured limit. MATLAB wall time may exceed
  the solver's internal time budget because a limit is not a hard external
  process deadline.
- PAR2 assigns measured time to a certified run and twice the comparison
  budget (or ablation target budget) otherwise. It is not a speedup ratio.

Gurobi solves the convex relaxation/screening QPs and constrained linear
least-squares subproblems used for feasible solutions and local searches.
QP numerical tolerances/safety margins are distinct from the global gap
tolerances. Unavailable QP bounds do not authorize pruning; the numerical
core retains valid inherited bounds and its original recovery checks.

## Regenerate paper tables without MATLAB or solvers

From the repository root, using Python's standard library only:

```bash
python3 scripts/verify_package.py
python3 -m unittest discover -s tests -v
python3 scripts/build_tables.py
```

This writes 3 comparison tables, 1 summary table, and 8 ablation tables to
`reports/`. They use normalized RSS, timeout marks, final timeout gaps,
node partitions and QP dimensions. The manuscript's visual formatting can
be applied around these generated tables. Load `amsmath` for the table
notation; no custom table macros are required.

To use newly computed outputs after `export_comparison`:

```bash
python3 scripts/build_tables.py \
  --comparison outputs/comparison/nist_tolerance_results.csv \
  --ablation outputs/ablation/all_trials.csv \
  --output reports/fresh
```

The QP dimension counts are original submitted models before presolve.
For `p` linear coefficients, `d` nonlinear variables, `n` residuals, and
`v=2^p` linear-box vertices, `m1` has `v*M^d+n` variables and `2*n+1`
linear rows. The screening QP uses `v*2^d+n` variables and the same rows;
its size is independent of `M`. Variable bounds are not constraint rows.

## Validation and reproducibility limits

`scripts/verify_package.py` checks frozen-source hashes, all published gap
certificates, pruning partitions, node partitions and memory exclusions.
GitHub Actions runs those checks and the table-export tests without solver
licenses. `test_portable_inputs` checks the MATLAB reconstruction without
calling an optimizer. `run_regression_tests` and full numerical runs require
MATLAB and licensed Gurobi, with Python/SCIP additionally required for SCIP.

The packaging checks and CSV/table tests were executed during preparation.
MATLAB/Gurobi/SCIP executions of these new portable entry points could not be
run in the packaging environment; run the supplied preflight and one-problem
commands locally before launching the entire battery. The frozen numerical
core has been used for the supplied experiment results. Exact timings,
node counts and successful completion can depend on machine performance,
solver builds and numerical behavior; the supplied observations are single
runs, not timing confidence intervals or formal interval-arithmetic proofs.

The original input generation routines and regression fixtures remain in
the versioned numerical folders for provenance. For this release, use the
portable commands above; obsolete examples and interactive launchers have
been excluded from this unified distribution.

## Citation, license and release

This package is distributed under the **PolyForm Noncommercial License
1.0.0**. See `../LICENSE` for its full terms and permitted purposes.
Uses outside those permissions require a separate written agreement;
commercial terms, including a license fee, may be negotiated with the rights
holders as described in `../COMMERCIAL_LICENSING.md`.

This is source-available software with restrictions on commercial use;
it is not distributed under an OSI-approved open-source license.
Third-party solvers, libraries and datasets retain their own terms.
See also `../CITATION.md` and `docs/PROVENANCE.md`.

For GitHub, upload the **contents** of this directory. Do not upload
`outputs/`, local Python environments, solver installations or license
files. `.gitignore` covers these generated/private files. A tagged release
and an archival DOI can be added once the authors publish the repository.
