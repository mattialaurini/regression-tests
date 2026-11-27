function test_nist_qp_recovery()
%TEST_NIST_QP_RECOVERY Reproduce the real failed region without a QP solver.
root=fileparts(mfilename('fullpath')); temporary=tempname; mkdir(temporary);
previousPath=path; previousFolder=pwd;
cleanup=onCleanup(@()localRestore(previousPath,previousFolder,temporary)); %#ok<NASGU>
source=sprintf(['function result=solve_nist_qp(varargin)\n' ...
    'global NIST_RECOVERY_TEST\nresult=NIST_RECOVERY_TEST.raw;\nend\n']);
localWrite(fullfile(temporary,'solve_nist_qp.m'),source);
source=sprintf(['function data=nist_build_relaxation_data(varargin)\n' ...
    'global NIST_RECOVERY_TEST\ndata=NIST_RECOVERY_TEST.data;\nend\n']);
localWrite(fullfile(temporary,'nist_build_relaxation_data.m'),source);
source=sprintf(['function raw=gurobi(model,params)\n' ...
    'global NIST_RECOVERY_TEST\n' ...
    'NIST_RECOVERY_TEST.retryCalls=NIST_RECOVERY_TEST.retryCalls+1;\n' ...
    'assert(params.NumericFocus==3 && params.DualReductions==0);\n' ...
    'if isfield(params,''Method'') && params.Method==1\nname=''dualRaw'';\n' ...
    'elseif model.Q(end,end)>1\nname=''preserveRaw'';\n' ...
    'else\nname=''scaledRaw'';\nend\n' ...
    'if isfield(NIST_RECOVERY_TEST,name)\nraw=NIST_RECOVERY_TEST.(name);\n' ...
    'else\nraw=struct(''status'',''INFEASIBLE'');\nend\nend\n']);
localWrite(fullfile(temporary,'gurobi.m'),source);
addpath(root); addpath(temporary,'-begin'); cd(temporary);
clear solve_nist_qp nist_build_relaxation_data gurobi
opts=struct('BoundSafety',1e-8);
global NIST_RECOVERY_TEST
NIST_RECOVERY_TEST=struct('retryCalls',0);
d=load(fullfile(root,'test_fixtures','mgh10_failed_qp.mat'));
d.nColumns=size(d.G,2);
NIST_RECOVERY_TEST.data=d;
NIST_RECOVERY_TEST.raw=struct('hasSolution',false,'status',"INF_OR_UNBD", ...
    'x',[],'objective',Inf,'objectiveBound',-Inf,'wallTime',0,'raw',struct());
result=solve_nist_outer_relaxation_qp_stable(d.G,d.cMinus,d.cPlus,opts);
assert(~result.boundAvailable && ~result.hasSolution);
assert(result.safeLowerBound==0);
assert(isinf(result.feasiblePrimalUpperBound));
m1=compute_nist_m1_bound_stable(struct(),struct(),3,opts);
assert(m1.value==0 && ~m1.boundAvailable);
m2=compute_nist_m2_screening_stable(struct(),struct(),2,5,opts);
assert(isinf(m2.value) && ~m2.screeningAvailable);

% Repair the simplex weights and ignore an invalid raw residual/objective.
NIST_RECOVERY_TEST.raw.hasSolution=true;
NIST_RECOVERY_TEST.raw.status="NUMERIC";
NIST_RECOVERY_TEST.raw.x=[-0.1;1.1;0];
NIST_RECOVERY_TEST.raw.objective=0;
result=solve_nist_outer_relaxation_qp_stable([2,3],zeros(1,2),zeros(1,2),opts);
assert(all(result.weights>=0) && sum(result.weights)==1);
assert(result.feasiblePrimalUpperBound>=9 && result.feasiblePrimalUpperBound<9.00001);
assert(~result.boundAvailable && result.safeLowerBound==0);

% No usable dual bound: a feasible primal iterate must not create a lower bound.
result=solve_nist_outer_relaxation_qp_stable([-2,2],zeros(1,2),zeros(1,2),opts);
assert(result.safeLowerBound==0);

% Preserve a successful finite QP dual bound with the original safety margin.
NIST_RECOVERY_TEST.raw.objectiveBound=2;
NIST_RECOVERY_TEST.raw.status="OPTIMAL";
NIST_RECOVERY_TEST.raw.objective=9;
NIST_RECOVERY_TEST.raw.raw=struct('pi',[-4;0;8]);
result=solve_nist_outer_relaxation_qp_stable([2,3],zeros(1,2),zeros(1,2),opts);
assert(result.boundAvailable && abs(result.safeLowerBound-(2-3e-8))<1e-14);
% Regression: an ordinary successful solve does not require pi or a retry.
NIST_RECOVERY_TEST.raw.solver="gurobi";
NIST_RECOVERY_TEST.raw.raw=struct();
result=solve_nist_outer_relaxation_qp_stable([2,3],zeros(1,2),zeros(1,2),opts);
assert(result.boundAvailable && result.qpRetryCount==0);
assert(NIST_RECOVERY_TEST.retryCalls==0);
assert(result.boundValidation.reason=="original_optimal_solver_policy");
% A finite but inconsistent value is not accepted as a pruning bound.
inconsistent=NIST_RECOVERY_TEST.raw;
inconsistent=rmfield(inconsistent,'solver');
inconsistent.objectiveBound=12;
NIST_RECOVERY_TEST.raw=inconsistent;
result=solve_nist_outer_relaxation_qp_stable([2,3],zeros(1,2),zeros(1,2),opts);
assert(~result.boundAvailable && result.safeLowerBound==0);
NIST_RECOVERY_TEST.raw.objectiveBound=NaN;
result=solve_nist_outer_relaxation_qp_stable([2,3],zeros(1,2),zeros(1,2),opts);
assert(~result.boundAvailable && result.safeLowerBound==0);
% Regression: m2 needs only the repaired primal, even with no usable bound.
NIST_RECOVERY_TEST.data=struct('G',[2,3],'cMinus',zeros(1,2), ...
    'cPlus',zeros(1,2),'nColumns',2);
NIST_RECOVERY_TEST.raw.status="NUMERIC";
NIST_RECOVERY_TEST.raw.solver="gurobi";
NIST_RECOVERY_TEST.raw.objectiveBound=-Inf;
m2=compute_nist_m2_screening_stable(struct(),struct(),2,5,opts);
assert(m2.screeningAvailable && m2.value>=9 && m2.qpRetryCount==0);
assert(~m2.boundAvailable && NIST_RECOVERY_TEST.retryCalls==0);
NIST_RECOVERY_TEST.raw.status="OPTIMAL";
NIST_RECOVERY_TEST.raw.objectiveBound=NaN;
m2=compute_nist_m2_screening_stable(struct(),struct(),2,5,opts);
assert(m2.screeningAvailable && m2.qpRetryCount==0);
assert(NIST_RECOVERY_TEST.retryCalls==0);
% m1 does not turn that same NUMERIC primal into a pruning bound.
NIST_RECOVERY_TEST.raw.status="NUMERIC";
NIST_RECOVERY_TEST.raw.objectiveBound=-Inf;
noRetryRaw=NIST_RECOVERY_TEST.raw;
noRetryRaw=rmfield(noRetryRaw,'solver');
NIST_RECOVERY_TEST.raw=noRetryRaw;
m1=compute_nist_m1_bound_stable(struct(),struct(),3,opts);
assert(~m1.boundAvailable && m1.value==0);
% Replay the actual scaled solver response received from the user's Gurobi.
actual=load(fullfile(root,'test_fixtures','mgh10_scaled_qp_solution.mat'));
NIST_RECOVERY_TEST.scaledRaw=actual.nativeRaw;
NIST_RECOVERY_TEST.raw=struct('solver',"gurobi",'hasSolution',false, ...
    'status',"INF_OR_UNBD",'x',[],'objective',Inf,'objectiveBound',-Inf, ...
    'wallTime',0,'raw',struct());
retryOpts=struct('BoundSafety',1e-8,'TimeLimit',10,'Display',false, ...
    'Threads',1,'FeasibilityTolerance',1e-8,'OptimalityTolerance',1e-8);
result=solve_nist_outer_relaxation_qp_stable(d.G,d.cMinus,d.cPlus,retryOpts);
assert(result.qpRetryCount==1 && result.boundAvailable && result.boundValidation.verified);
assert(result.safeLowerBound>4.7e19 && result.safeLowerBound<4.73e19);
assert(result.originalQPStatus=="INF_OR_UNBD" && result.status=="OPTIMAL");
assert(result.objectiveBound<=result.feasiblePrimalUpperBound);
assert(NIST_RECOVERY_TEST.retryCalls==1);
% The recovered primal suffices for m2 even when no dual iterate is returned.
NIST_RECOVERY_TEST.scaledRaw=rmfield(actual.nativeRaw,'pi');
NIST_RECOVERY_TEST.data=d;
m2=compute_nist_m2_screening_stable(struct(),struct(),2,2,retryOpts);
assert(m2.screeningAvailable && ~m2.boundAvailable && m2.qpRetryCount==1);
% Replay the four ACTUAL numerical cases supplied by the user on Lanczos2/Hahn1.
captured=load(fullfile(root,'test_fixtures','qp_recovery_examples.mat'));
for k=1:numel(captured.examples)
    example=captured.examples{k};
    NIST_RECOVERY_TEST.raw=localWrapCaptured(example.originalRaw);
    NIST_RECOVERY_TEST.dualRaw=example.dualRaw;
    NIST_RECOVERY_TEST.preserveRaw=example.preserveRaw;
    NIST_RECOVERY_TEST.scaledRaw=example.scaledRaw;
    NIST_RECOVERY_TEST.retryCalls=0;
    bound=solve_nist_outer_relaxation_qp_stable( ...
        example.G,example.cMinus,example.cPlus,retryOpts);
    assert(bound.boundAvailable && bound.boundValidation.verified);
    assert(bound.safeLowerBound>example.expectedMinimum);
    assert(bound.qpRetryCount==example.expectedRetries);
    assert(NIST_RECOVERY_TEST.retryCalls==example.expectedRetries);
    assert(bound.safeLowerBound<=bound.feasiblePrimalUpperBound);
    if k==3
        % A slightly infeasible raw primal objective must not discard this
        % independently checked original bound or trigger any re-solve.
        assert(bound.safeLowerBound>bound.objective);
        assert(~bound.boundValidation.raw_primal_bound_consistent);
    end
end
% A recovery that yields only zero must not hide the later useful bound.
example=captured.examples{1};
NIST_RECOVERY_TEST.raw=localWrapCaptured(example.originalRaw);
NIST_RECOVERY_TEST.dualRaw=example.scaledRaw;
NIST_RECOVERY_TEST.preserveRaw=example.preserveRaw;
NIST_RECOVERY_TEST.scaledRaw=example.scaledRaw;
bound=solve_nist_outer_relaxation_qp_stable( ...
    example.G,example.cMinus,example.cPlus,retryOpts);
assert(bound.qpRetryCount==2 && bound.boundAvailable && bound.safeLowerBound>0.001);
fprintf('QP recovery tests passed (2.5.0): real Lanczos2/Hahn1 responses, feasible-primal consistency, useful recovery bounds, legacy MGH10 and separate screening.\n');
end
function native=localWrapCaptured(raw)
native=struct('solver',"gurobi",'status',string(raw.status),'hasSolution',isfield(raw,'x'), ...
    'x',[],'objective',Inf,'objectiveBound',-Inf,'wallTime',0,'raw',raw);
if native.hasSolution
    native.x=double(raw.x(:)); native.objective=double(raw.objval);
end
if isfield(raw,'objbound'), native.objectiveBound=double(raw.objbound); end
end
function localWrite(file,source)
fid=fopen(file,'w'); assert(fid>=0); fwrite(fid,source,'char'); fclose(fid);
end
function localRestore(previousPath,previousFolder,temporary)
cd(previousFolder); path(previousPath);
clear solve_nist_qp nist_build_relaxation_data gurobi
clear global NIST_RECOVERY_TEST
if isfolder(temporary), rmdir(temporary,'s'); end
end
