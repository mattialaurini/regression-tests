# Quick start

Open MATLAB in the repository root and add your installed Gurobi MATLAB
interface to the MATLAB path. Then run:

```matlab
setup_reproduction;
test_portable_inputs;
run_regression_tests;
[resultsBnB, summaryBnB] = run_comparison_bnb;
```

To run SCIP, first select your Python interpreter:

```matlab
setup_reproduction(fullfile(pwd,'.venv','bin','python'));
[resultsSCIP, summarySCIP] = run_comparison_scip;
export_comparison;
```

On Windows use `fullfile(pwd,'.venv','Scripts','python.exe')` instead.
See README.md for Python dependency installation.

For the ablation study:

```matlab
[resultsAblation, summaryAblation] = run_ablation;
```

No file-selection dialogs are used. All eight problem instances are bundled.
New results are written to `outputs/`; published results are kept separately
in `reference/`. Repeat the same command to resume completed checkpoints.

Before launching the complete battery, use `DryRun=true` or select one
problem with `CaseIndices=1,ProblemIds="BoxBOD"` and a dedicated output
directory. The README explains dependencies, budget options, counters and
table generation. The software is distributed under PolyForm Noncommercial 1.0.0; see
`LICENSE` and `COMMERCIAL_LICENSING.md` for permitted uses and separate
commercial licensing.
