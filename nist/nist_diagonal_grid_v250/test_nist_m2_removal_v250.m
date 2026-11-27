function test_nist_m2_removal_v250()
%TEST_NIST_M2_ABLATION Deterministic tests against temporary subproblem doubles.
% Any m2 invocation raises an error. No optimization solver is needed here.
root=fileparts(mfilename('fullpath'));
previousPath=path; previousFolder=pwd; temporary=tempname; mkdir(temporary);
cleanup=onCleanup(@()localCleanup(previousPath,previousFolder,temporary)); %#ok<NASGU>
files=dir(fullfile(root,'test_fixtures','*.m.txt'));
for i=1:numel(files)
    copyfile(fullfile(files(i).folder,files(i).name),fullfile(temporary,erase(files(i).name,'.txt')));
end
addpath(root); addpath(temporary,'-begin'); cd(temporary); localClear();
global NIST_TOLERANCE_TEST_CASE
source=fileread(fullfile(root,'reference','solve_nist_bnb_tolerance_study.m.txt'));
tokens=regexp(source,'options\.([A-Za-z]\w*)','tokens');
names=unique(string([tokens{:}])); opts=struct();
for name=names, opts.(name)=0; end
opts.TimeLimit=Inf; opts.NodeLimit=Inf;
opts.EpsilonRelative=0.1; opts.EpsilonAbsolute=1e-10;
opts.CoarseGridPoints=2; opts.InitialApproxGridPoints=2; opts.MaximumFineGridPoints=2;
opts.BranchingRule="normalized_longest_edge"; opts.StrongBranchingTieBreak=false;
opts.TieTolerance=1e-12; opts.RecordRelativeTargets=[0.1,0.01,0.001];
opts.Verbosity=0; opts.ProgressIntervalSeconds=Inf; opts.QPOptions=struct('TimeLimit',Inf);
opts.EnableLocalSearch=false; opts.LocalSearchAtRoot=false;

% Tolerance-only discard at the inclusive threshold, with m2 forbidden.
r=localRun("constant",0.9-1e-10,opts,1);
assert(r.paper_certificate_satisfied && r.nodes==1);
assert(r.stats.prunedToleranceOnly==1 && r.stats.m1Evaluations==1);
assert(r.stats.m2Evaluations==0 && r.stats.m1SkippedByScreening==0);

% Independent no-m2 tuning starts directly on K, without grid reinsertion.
for K=[3,5]
    fixed=opts; fixed.InitialApproxGridPoints=K; fixed.MaximumFineGridPoints=K;
    r=localRun("constant",0.95,fixed,1);
    assert(r.paper_certificate_satisfied && r.stats.m1Evaluations==1);
    assert(r.stats.nodesReinserted==0 && r.stats.m2Evaluations==0);
    assert(isequal(NIST_TOLERANCE_TEST_CASE.fineGrids,K));
end

% Ordinary incumbent-based pruning stays distinct.
r=localRun("constant",1,opts,1);
assert(r.stats.prunedAtIncumbent==1 && r.stats.prunedToleranceOnly==0);

% The original fine-grid cap/reinsertion policy is retained.
reinsert=opts; reinsert.InitialApproxGridPoints=5;
r=localRun("constant",0.95,reinsert,1);
assert(r.nodes==1 && r.queue_extractions==3 && r.stats.nodesReinserted==2);
assert(r.stats.m1Evaluations==1 && r.stats.m2Evaluations==0);

% A child objective improvement closes the remaining inherited-bound region.
r=localRun("global_queue",0.24,opts,1);
assert(r.nodes==2 && r.stats.nodesCreated==3);
assert(r.stats.prunedAtIncumbentExplicit==1);
assert(r.stats.prunedToleranceOnlyGlobal==1 && r.stats.prunedTotal==2);

% A node limit must leave unresolved queue regions unclassified.
limited=opts; limited.NodeLimit=1;
r=localRun("global_queue",0.24,limited,1);
assert(~r.paper_certificate_satisfied && r.stats.prunedTotal==0);
assert(r.stats.m2Evaluations==0);

% Initial certification need not visit or evaluate any node.
r=localRun("constant",0.95,opts,1e-12);
assert(r.nodes==0 && r.stats.m1Evaluations==0 && r.stats.prunedToleranceOnlyGlobal==1);

% Local search remains enabled, including root and stagnation calls.
local=opts; local.NodeLimit=3; local.EnableLocalSearch=true; local.LocalSearchAtRoot=true;
local.LocalSearchMaxCalls=3; local.LocalSearchMaxTimePerCall=2;
local.LocalSearchMaxTotalTimeFraction=1; local.LocalSearchCooldownNodes=1;
local.LocalSearchStagnationNodes=0; local.LocalSearchPoorUBTrigger=false;
r=localRun("constant",0.2,local,1);
assert(r.stats.localSearchRootCalls==1 && r.stats.localSearchPeriodicCalls==2);
assert(r.stats.localSearchStagnationTriggerCalls==2 && NIST_TOLERANCE_TEST_CASE.localCalls==3);
assert(r.stats.m1Evaluations==3 && r.stats.m2Evaluations==0);

% The poor-UB trigger uses an inherited valid bound when m2 is absent.
local.LocalSearchStagnationNodes=Inf; local.LocalSearchPoorUBTrigger=true;
local.LocalSearchPoorUBMinNodes=1; local.LocalSearchPoorUBRatio=0.3;
r=localRun("constant",0.2,local,1);
assert(r.stats.localSearchRootCalls==1 && r.stats.localSearchPoorUBTriggerCalls==2);
assert(r.stats.localSearchCalls==3 && r.stats.m2Evaluations==0);
% Actual baseline grid controls 5 -> 3 are retained, not replaced by a fixed grid.
actualGrid=opts; actualGrid.InitialApproxGridPoints=5; actualGrid.MaximumFineGridPoints=3;
r=localRun("constant",0.95,actualGrid,1);
assert(r.nodes==1 && r.queue_extractions==2 && r.stats.nodesReinserted==1);
assert(isequal(NIST_TOLERANCE_TEST_CASE.fineGrids,3));
assert(r.options.InitialApproxGridPoints==5 && r.options.MaximumFineGridPoints==3);

% A failed m1 must not prune the root; descendants may still be certified.
r=localRun("failed_root",0,opts,1);
assert(r.paper_certificate_satisfied && r.nodes==3);
assert(r.stats.m1UnavailableBounds==1 && r.stats.nodesPruned==2);
assert(r.stats.m2Evaluations==0);

% Adjacent floating-point endpoints must not produce non-progressing splits.
precision=struct('zLower',1,'zUpper',1+eps(1));
r=localRun("constant",0.2,opts,1,precision);
assert(r.status=="numerical_resolution_limit" && ~r.paper_certificate_satisfied);
assert(r.stats.precisionLimitedNodes==1 && r.stats.nodesBranched==0);
assert(r.stats.m2Evaluations==0);
fprintf('NIST m2 removal v2.5.0 self-tests passed: m2 forbidden, original grids, root/periodic/poor-UB search, failures and precision limits.\n');
end

function result=localRun(mode,lower,options,initialUpper,region)
if nargin<5, region=struct('zLower',0,'zUpper',1); end
global NIST_TOLERANCE_TEST_CASE
NIST_TOLERANCE_TEST_CASE=struct('mode',mode,'lower',lower,'options',options,'localCalls',0,'fineGrids',zeros(1,0));
problem=struct('id',"test",'instance_id',"test_study", ...
    'nLinear',1,'nNonlinear',1,'nResidual',1,'domainProfile',"study", ...
    'zLower',region.zLower,'zUpper',region.zUpper,'objectiveScale',1,'residualScale',1,'certifiedRSS',1, ...
    'initial',struct('linear',sqrt(initialUpper),'z',region.zLower+(region.zUpper-region.zLower)/2));
result=solve_nist_bnb_m2_removed(problem,options);
end
function localClear()
clear nist_bnb_options validate_nist_problem evaluate_nist_objective
clear solve_nist_linear_block_qp compute_nist_m1_bound_stable compute_nist_m2_screening_stable
clear solve_nist_reduced_local_search
end
function localCleanup(previousPath,previousFolder,temporary)
cd(previousFolder); path(previousPath); localClear();
clear global NIST_TOLERANCE_TEST_CASE
if isfolder(temporary), rmdir(temporary,'s'); end
end
