function preflight = check_nist_qp_recovery(qpOptions)
%CHECK_NIST_QP_RECOVERY Real-solver checks on five captured QPs; no BnB.
arguments
    qpOptions (1,1) struct = nist_qp_options()
end
root=fileparts(mfilename('fullpath'));
captured=load(fullfile(root,'test_fixtures','qp_recovery_examples.mat'));
labels=["Lanczos2_depth30","Lanczos2_depth45","Lanczos2_depth75","Hahn1_depth21","MGH10_failed_box"];
qpOptions.TimeLimit=min(double(qpOptions.TimeLimit),5);
qpOptions.Display=false;
rows=cell(5,1); details=cell(5,1);
for k=1:5
    if k<=4
        example=captured.examples{k}; minimum=example.expectedMinimum;
    else
        example=load(fullfile(root,'test_fixtures','mgh10_failed_qp.mat'));
        minimum=4.7e19;
    end
    result=solve_nist_outer_relaxation_qp_stable( ...
        example.G,example.cMinus,example.cPlus,qpOptions);
    dual=result.raw;
    if isfield(dual,'pi')
        dual.pi=double(dual.pi)*(result.nativeObjectiveFactor/result.nativeResidualScale);
    end
    [checked,validation]=nist_verify_qp_dual(example.G+example.cMinus, ...
        example.G+example.cPlus,dual,result.safeLowerBound,1);
    passed=result.boundAvailable && checked && result.safeLowerBound>minimum && ...
        result.safeLowerBound<=result.feasiblePrimalUpperBound;
    rows{k}=table(labels(k),string(result.status),result.safeLowerBound, ...
        minimum,result.qpRetryCount,logical(checked),logical(passed), ...
        'VariableNames',{'case_id','status','lower_bound','expected_minimum', ...
        'recovery_calls','dual_checked','passed'});
    details{k}=struct('result',result,'validation',validation);
end
report=vertcat(rows{:});
preflight=struct('report',report,'details',{details},'passed',all(report.passed), ...
    'time_limit_per_qp',qpOptions.TimeLimit);
disp(report);
if preflight.passed
    fprintf('Real Gurobi QP preflight passed (5 captured cases, version 2.5.0).\n');
else
    fprintf('Real Gurobi QP preflight FAILED: timed benchmarks must not start.\n');
end
end
