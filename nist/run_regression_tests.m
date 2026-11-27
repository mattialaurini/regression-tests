function run_regression_tests()
%RUN_REGRESSION_TESTS Deterministic tests, then tests requiring licensed Gurobi.
state=prepare_nist_reproduction();
test_portable_inputs();
oldPath=path; cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(state.baseRoot,'matlab')); addpath(state.toleranceRoot,'-begin');
test_nist_tolerance_instrumentation(); test_nist_qp_recovery();
report=check_nist_qp_recovery(state.baseOptions{1}.QPOptions);
assert(report.passed,'Gurobi recovery preflight failed.');
addpath(state.ablationRoot,'-begin');
test_nist_diagonal_node_counts(); test_nist_fine_screening_v250(); test_nist_m2_removal_v250();
fprintf('Reproduction regression tests passed.\n');
end
