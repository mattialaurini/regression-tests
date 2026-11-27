# Quick start

Open MATLAB in this directory, add your licensed Gurobi MATLAB interface
to the path, and run:

```matlab
setup_experiments;
run_pkpd_arx(DryRun=true);
run_pkpd_ml(DryRun=true);
run_pkpd_patient13(DryRun=true);
```

Regenerate the supplied patient-13 illustration without fitting:

```matlab
export_pkpd_reference;
```

Recompute the PK/PD experiments when ready:

```matlab
run_pkpd_patient13;
run_pkpd_arx;
run_pkpd_ml;
export_pkpd_results;
generate_pkpd_profiles;
```

The full cohort batches can take many hours. Their original result MAT
files are not bundled; patient-13 fits are supplied separately.

The NIST commands are unchanged after `setup_experiments`:

```matlab
run_comparison_bnb;
run_comparison_scip;
run_ablation;
```

Configure Python before the SCIP command as explained in README.md.
Use a fresh MATLAB session for timings in each suite. See the README for
dependencies, subsets, checkpoints, budgets, licensing and validation limits.
