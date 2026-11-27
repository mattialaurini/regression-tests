# Changelog

## 1.0.6
- Replaced node-count-based BnB progress printing with wall-clock throttling.
- Added `ProgressIntervalSeconds` to `nist_bnb_options` (default: 20 s; minimum finite value: 15 s).
- `Verbosity >= 2` now emits at most one periodic BnB progress line per configured interval, independent of node throughput.
- Updated the MGH10 local-search example to use a 20 s progress interval.


## 1.0.5
- Added a budgeted reduced-space local-search heuristic inside BnB.
- The heuristic is always invoked at the root by default and can be invoked periodically on promising nodes.
- Periodic invocation uses two adaptive triggers:
  - incumbent stagnation for a configurable number of processed nodes;
  - a poor-UB signal when the current node's `m2` value is much smaller than the incumbent objective.
- The poor-UB signal is used only to schedule a primal heuristic; `m2` is never treated as a valid lower bound for certification or pruning.
- Local search optimizes only the nonlinear block with bounded `fmincon`; the linear block is solved exactly at every trial point using `solve_nist_linear_block_qp`.
- Added explicit local-search budgets: cooldown, maximum calls, per-call function evaluations/iterations/time, and maximum fraction of total BnB wall time.
- Added local-search statistics including call counts, trigger counts, time, function evaluations, QP time, and incumbent improvements.
- Retains normalized longest-edge branching from v1.0.4.
- Added `matlab/examples/run_mgh10_bnb_local_search.m` for the MGH10 regression test.
- VarPro remains a completely independent baseline and is not used to initialize or warm-start BnB or SCIP.

## 1.0.4
- Added scale-invariant `normalized_longest_edge` branching for the nonlinear block.
- Kept `absolute_longest_edge` as an explicit compatibility/ablation option.
- Added BnB statistics `branchCountByDimension`, `maxBranchDepthByDimension`, `maxProcessedBranchDepthByDimension`, and `maxNodeDepth`.
- Added `matlab/examples/run_mgh10_normalized_branching.m` for an isolated MGH10 regression test.
- No changes to m1/m2 relaxations, QP settings, or stopping criteria.

## 1.0.3
- Keeps VarPro strictly as an independent local baseline; it is not used to initialize BnB or SCIP.
- Classifies paired outcomes by certification: if only one global solver certifies, that solver wins; if neither certifies, the outcome is `both_uncertified`.
- Adds `both_uncertified` to the paired summary.
- Adds `diagnose_nist_bnb_bounds.m`, a no-rerun diagnostic for near-zero global lower bounds using saved BnB histories.
- Includes the MATLAB R2026a-compatible `nist_benchmark_config.m` hotfix.

## 1.0.1
- Fixed `nist_basis_second_derivative_intervals.m` for `Gauss1`, nonlinear coordinate 1 (`b2`).
- Removed an early return that bypassed the common `residualScale` division.
- No changes to datasets, model formulas, parameter domains, SCIP model, or the other seven benchmark problems.
