# Alignment with the latest supplied manuscript

The manuscript update changes explanation and figure interpretation, not
the scope or numerical settings of the reported experiments. The following
experiment families remain in the package:

- ARX approximate-likelihood identification, grid tuning and screening comparison;
- exact-ML double-pole and one-compartment identification, screening on/off;
- the six-model patient-13 illustration at sigma=1;
- NIST BnB/SCIP comparison and matched-grid NIST screening ablation.

No additional numerical battery or known-baseline ML comparison is required.
The existing derivative, containment, QP and instrumentation regression checks
remain useful for validating these implementations.

## Documentation and comments

- Identify the profile table as the Padula et al. (2015), Table 1, source cited
  by the manuscript, with the literature set traced to Ionescu et al. (2008).
- Specify demographic, concentration, infusion and Hill-map parameter units.
- Distinguish demographic Schnider PK parameters from tabulated Hill-map
  parameters. Mismatch perturbs the former; the latter stay fixed during
  response generation.
- Retain profile 13 exactly as tabulated.
- Explain why exact ML estimates E0 and Emax together, while inverse-Hill ARX
  uses a fixed measured baseline. Both ML coefficients remain linear.
- Correct clearance comments to L/min; the existing division by 60 already
  computes rate constants in 1/s. No numerical conversion was changed.

## Plotting-only correction

The updated ARX figure caption states that exact zeros are displayed below
the smallest positive error. The earlier exporter instead floored every
error below 1e-4. The exporter now preserves all positive errors and moves
only exact zeros to a shared positive plotting level. This affects figure
display only; fits, certificates, tables and source statistics are unchanged.
`test_pkpd_plot_errors` checks this behavior in MATLAB.
The patient-13 input-axis label now includes mg/s.

## Outstanding source limitations

The supplied final-batch defaults still include four noise levels although
one manuscript sentence says two; the defaults have not been silently
changed. Full original ARX/ML cohort MAT files remain absent from the supplied
archive. These two points do not change the retained implementation set.
