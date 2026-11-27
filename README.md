# Separable regression tests

Reproducible numerical experiments for **An interpolation-based branch-and-bound
method for separable nonlinear regression: an application to PK/PD model
identification**.

This repository contains the manuscript's PK/PD identification experiments,
the NIST comparison against SCIP, and the matched-grid NIST ablation.

## Layout

| Folder | Contents |
|---|---|
| `pkpd/src/` | Patient simulation, ARX/ML identification, bounds, grid tuning and manuscript exporters |
| `pkpd/reference/` | The six saved patient-13 fits, metrics and sampled signals |
| `nist/` | Portable NIST input data, BnB/SCIP drivers, ablation and reference results |
| `docs/` | Experiment-to-paper mapping, cleanup provenance and validation scope |
| `licenses/` | Third-party license notices |

Development benchmarks, obsolete identification experiments, real-data
scripts, Simulink-era inputs, OS metadata and duplicate binary figure files
have been excluded. The numerical algorithms are retained; English API names,
isolated launchers and explanatory comments are packaging changes.

## Requirements

- MATLAB; the recorded experiment environment was R2026a Update 3 on an
  Apple M2 Pro with 16 GiB RAM, using Gurobi 12.0.2.
- **Control System Toolbox** for PK/PD simulation and the nominal population
  model (`ss`, `c2d`, `ssdata`). No Simulink model is required.
- **Gurobi's MATLAB interface**, installed and licensed separately.
- **Optimization Toolbox** for the NIST BnB local searches (`fmincon`).
- For the SCIP comparison only: Python 3.12.15, NumPy 2.5.3,
  PySCIPOpt 6.2.1 and native SCIP 10.0.2, matching the recorded environment.

PK/PD uses the original hybrid projected-gradient/Gurobi solver, sequential
execution, automatic Gurobi thread selection and automatic presolve. NIST
uses its separately recorded one-thread configuration and disables SCIP
presolve. Keep these configurations separate when interpreting times.
For timing comparisons, run each suite in a fresh MATLAB session without
other concurrent numerical jobs, using a consistent power setting.

## Set up

Extract the repository into a new directory, make that directory MATLAB's
current folder, and add your installed Gurobi MATLAB interface to the path.
Then run:

```matlab
setup_experiments;
run_pkpd_arx(DryRun=true);
run_pkpd_ml(DryRun=true);
run_pkpd_patient13(DryRun=true);
```

Dry runs list the planned experiments without calling an optimizer.
No old result files or file-selection dialogs are required.

## PK/PD experiments

### Saved six-model illustration: no optimization required

```matlab
export_pkpd_reference;
```

This regenerates `patient13_six_models_BIS.pdf`,
`patient13_six_models_errors.pdf`, and their CSV data under
`outputs/pkpd_reference/` from the supplied fits. It needs MATLAB plotting
only; it does not run identification or require Control System Toolbox.

### Recompute the six-model illustration

```matlab
R = run_pkpd_patient13;
```

This fits ARX orders 1--4, exact-ML double pole, and exact-ML one compartment
for patient 13 with `sigma=1`, seed 13 and PK mismatch `omega=0.2`.
The relative/absolute gap tolerances are `0.05` and `1e-3`, on the residual
norm. ARX uses `(N,M,Nhat)=(3,5,5)` and a 1800-s limit per fit; ML uses a
3600-s limit, with `(3,3,3)` without screening for the double pole and
`(2,3,3)` with screening for one compartment. E0 is measured/fixed for ARX
and estimated for ML. Outputs are saved to `outputs/pkpd_patient13/`.
Every fit is checkpointed before figures are exported.

### ARX and exact-ML cohort experiments

```matlab
arx = run_pkpd_arx;
ml  = run_pkpd_ml;
export_pkpd_results;
generate_pkpd_profiles;
```

The ARX driver first tunes the grids independently with/without screening,
then runs the selected settings for all 13 profiles and orders 1--4.
The ML driver compares screening on/off for the double-pole and
one-compartment models on the same patient/noise realizations.

The defaults retain the supplied final-batch configuration:

| Setting | ARX | Exact ML |
|---|---|---|
| Patients | 1--13 | 1--13 |
| Noise standard deviations | 0, 0.01, 0.1, 1 | 0, 0.01, 0.1, 1 |
| PK mismatch | 0.2 | 0.2 |
| Relative/absolute norm tolerance | 0.05 / 0.001 | 0.05 / 0.001 |
| Time limit per main run | 1800 s | 3600 s |
| Main runs | 416 | 208 |
| Tuning | Patients 1, 6, 9; sigma 0.1; adaptive grids; nominal limit 600 s | Fixed grids from the manuscript |

The source and manuscript descriptions have a noise-level wording mismatch:
the prose mentions two levels, but the supplied final batch defaults contain
four. The four-level source configuration is preserved rather than silently
changing it. To select only the positive levels highlighted in the discussion:

```matlab
run_pkpd_arx(Sigmas=[0.1,1],OutputDirectory=fullfile(pwd,'outputs','arx_two_levels'));
run_pkpd_ml(Sigmas=[0.1,1],OutputDirectory=fullfile(pwd,'outputs','ml_two_levels'));
```

This produces 208 ARX and 104 ML main runs; it is a selected experiment and
changes the corresponding counts and generated tables.

The complete cohort result MAT files are **not present in the supplied
archive**. The package includes all code and patient inputs needed to
recompute them, but cannot regenerate their published tables without either
running the batches or supplying those original MAT files to
`export_pkpd_results(arxFile,mlFile,outputDirectory)`.

The exporters write `Tables/tab_lin_tuning.tex`, `tab_lin_results.tex`,
`tab_ml_results.tex`, `lin_numbers.tex`, `ml_numbers.tex`, and the ARX
figures under `outputs/paper/`. For logarithmic ARX error plots, only exact zeros are moved below the
smallest positive error; positive errors and statistical summaries retain
their original values. The generated LaTeX uses the manuscript's
notation/macros, such as `\C`; include it in the manuscript, not as a
standalone document without that preamble.

### Select or resume a PK/PD batch

```matlab
run_pkpd_ml(PatientIds=13,Sigmas=1, ...
    OutputDirectory=fullfile(pwd,'outputs','ml_patient13_check'),Export=false);
```

Each cohort batch runs inside its own function workspace, avoiding leakage
of `cfg`, `settings` or `settings_ml` from the caller. Settings and source
snapshots are saved; changing them requires a new output directory.
Repeating the same command resumes completed results. `Resume=false`
archives the old result file before starting again. An interrupted ARX
tuning phase without a completed grid selection must be restarted with
`Resume=false`; no incomplete selection is silently accepted.

Randomness is generated with MATLAB's Mersenne Twister, seed=patient index.
The seven PK perturbations are drawn before the measurement noise. ARX clips
samples above E0 before inverse-Hill evaluation; exact ML uses raw samples.
Profile 13 is used exactly as tabulated by Padula et al. (2015), Table 1,
as identified in the accompanying manuscript, rather than recomputed by
averaging. The tabulated Hill-map parameters are kept fixed during data
generation; only PK parameters are perturbed. The default input is in mg/s.
PK/PD objectives and gap tolerances refer to a **residual norm**, not RSS/s.
The reported iteration count is the original main-loop count and is not
relabelled as the distinct-region count used by the NIST suite.

## NIST experiments

Run `setup_experiments` first, then use the already validated commands:

```matlab
test_portable_inputs;
run_regression_tests;
[resultsBnB, summaryBnB] = run_comparison_bnb;
```

For SCIP, install the Python requirements in
`nist/requirements-scip.txt`, select that interpreter, then run:

```matlab
setup_experiments(fullfile(pwd,'.venv','bin','python'));
[resultsSCIP, summarySCIP] = run_comparison_scip;
export_comparison;
```

On Windows substitute `.venv/Scripts/python.exe`.
For the matched-grid screening ablation:

```matlab
[resultsAblation, summaryAblation] = run_ablation;
```

The comparison uses relative tolerances 0.1/0.01/0.001 and budgets
100/600/3600 s. The ablation uses `N=2`, `M=Nhat in {3,5,9,17}` and the
per-trial limits recorded in its reference CSV. The 8-GiB planning guard
excludes Gauss1 at M=9 and 17 before allocation. Raw wall times, actual and
target budgets, final gaps, Ntol and cheap/costly node counts are preserved.
All NIST objectives/tolerances are on **RSS/s**, with
`s=max(NIST certified RSS,1e-8)`.

New NIST results remain under `nist/outputs/`. For detailed configuration,
report fields and table regeneration, see `nist/README.md`.

## Checks and validation scope

Without MATLAB or solver licenses:

```bash
python3 scripts/verify_distribution.py
python3 nist/scripts/verify_package.py
python3 -m unittest discover -s nist/tests -v
python3 -m unittest discover -s tests -v
```

In MATLAB, the additional PK/PD numerical checks are:

```matlab
setup_experiments;
run_pkpd_checks;
```

They check the six saved outputs plus derivative/containment routines.
See `docs/VALIDATION.md` for what has been executed and what still needs
local MATLAB execution. The unified PK/PD launchers have not been executed
in the packaging environment, which lacks MATLAB and licensed solvers.
Successful smoke tests do not establish that every full-battery case or
every run time reproduces exactly on another machine.

## License and citation

Project-owned code is distributed under **PolyForm Noncommercial 1.0.0**;
see `LICENSE` and `COMMERCIAL_LICENSING.md`. Third-party components retain
their own terms, notably the BSD-licensed priority queue listed in
`THIRD_PARTY_NOTICES.md`. See `CITATION.md` for the accompanying manuscript.

The suggested repository name is `regression-tests`, covering both suites.
Copy the contents of this folder to GitHub. Do not publish local solver
installations, virtual environments, license keys, or newly generated run
folders; `.gitignore` excludes those paths. Add the repository URL and the
Zenodo release DOI once the authors publish and archive the release.

## Archived release

Version 1.1.1 is archived on Zenodo:
[https://doi.org/10.5281/zenodo.24583728](https://doi.org/10.5281/zenodo.24583728).
