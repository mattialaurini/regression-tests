function test_nist_fine_screening_v250()
%TEST_NIST_FINE_SCREENING_V250 Exercise unchanged solver with isolated doubles.
root=fileparts(mfilename('fullpath'));
previousPath=path; previousFolder=pwd; temporary=tempname; mkdir(temporary);
cleanup=onCleanup(@()localCleanup(previousPath,previousFolder,temporary)); %#ok<NASGU>
files=dir(fullfile(root,'screening_test_fixtures','*.m.txt'));
for i=1:numel(files)
    copyfile(fullfile(files(i).folder,files(i).name),fullfile(temporary,erase(files(i).name,'.txt')));
end
addpath(root); addpath(temporary,'-begin'); cd(temporary); localClear();
global NIST_TOLERANCE_TEST_CASE
source=fileread(fullfile(root,'solve_nist_bnb_tolerance_study.m'));
tokens=regexp(source,'options\.([A-Za-z]\w*)','tokens');
names=unique(string([tokens{:}])); opts=struct();
for name=names, opts.(name)=0; end
opts.TimeLimit=Inf; opts.NodeLimit=Inf;
opts.EpsilonRelative=0.1; opts.EpsilonAbsolute=1e-10;
opts.CoarseGridPoints=2; opts.InitialApproxGridPoints=641; opts.MaximumFineGridPoints=3;
opts.BranchingRule="normalized_longest_edge"; opts.StrongBranchingTieBreak=false;
opts.TieTolerance=1e-12; opts.RecordRelativeTargets=[0.1,0.01,0.001];
opts.Verbosity=0; opts.ProgressIntervalSeconds=Inf; opts.QPOptions=struct('TimeLimit',Inf);
opts.EnableLocalSearch=false; opts.LocalSearchAtRoot=false;

% A failed screening test coarsens M but never constructs a large m1 QP.
r=localRun("constant",0.95,1.1,opts);
assert(r.paper_certificate_satisfied && r.nodes==1 && r.queue_extractions==9);
assert(r.stats.nodesReinserted==8 && r.stats.m1Evaluations==1);
assert(isequal(NIST_TOLERANCE_TEST_CASE.fineGrids,3));
assert(isequal(NIST_TOLERANCE_TEST_CASE.screenGrids,[641,321,161,81,41,21,11,6,3]));
assert(r.stats.m2Evaluations==9 && r.stats.prunedToleranceOnly==1);

% Every sweep grid reaches exactly 3 with the original coarsening rule.
for M=[5,17,65,257]
    sweep=opts; sweep.InitialApproxGridPoints=M;
    r=localRun("constant",0.95,1.1,sweep);
    assert(r.paper_certificate_satisfied && r.stats.m1Evaluations==1);
    assert(isequal(NIST_TOLERANCE_TEST_CASE.fineGrids,3));
    assert(NIST_TOLERANCE_TEST_CASE.screenGrids(1)==M);
    assert(NIST_TOLERANCE_TEST_CASE.screenGrids(end)==3);
end

% Joint study: exercise the actual solver for every M/N/Nhat tuple.
for M=[5,17,65]
    for N=[2,3]
        for H=[3,5]
            joint=opts; joint.InitialApproxGridPoints=M;
            joint.CoarseGridPoints=N; joint.MaximumFineGridPoints=H;
            r=localRun("constant",0.95,1.1,joint);
            assert(r.paper_certificate_satisfied && r.stats.m1Evaluations==1);
            assert(isequal(NIST_TOLERANCE_TEST_CASE.fineGrids,H));
            assert(all(mod(NIST_TOLERANCE_TEST_CASE.screenGrids-1,N-1)==0));
        end
    end
end

% Successful screening branches before constructing any m1 relaxation.
limited=opts; limited.NodeLimit=1;
r=localRun("constant",0.2,0.5,limited);
assert(~r.paper_certificate_satisfied && r.stats.m1Evaluations==0);
assert(r.stats.m1SkippedByScreening==1 && r.stats.nodesBranched==1);
assert(isempty(NIST_TOLERANCE_TEST_CASE.fineGrids));
assert(all(NIST_TOLERANCE_TEST_CASE.screenGrids==641));

% An unavailable m2 forces the m1 path and cannot justify discarding a node.
r=localRun("constant",0.95,Inf,opts);
assert(r.paper_certificate_satisfied && r.stats.m2UnavailableEvaluations==9);
assert(r.stats.m1Evaluations==1 && r.stats.m1SkippedByScreening==0);
assert(isequal(NIST_TOLERANCE_TEST_CASE.fineGrids,3));

% Failed root lower bound branches; child bounds, not failure, certify.
r=localRun("failed_root",0,1.1,opts);
assert(r.paper_certificate_satisfied && r.nodes==3 && r.stats.m1UnavailableBounds==1);
assert(r.stats.nodesPruned==2 && all(NIST_TOLERANCE_TEST_CASE.fineGrids==3));

% Root and periodic local searches remain available with large initial M.
local=opts; local.NodeLimit=3; local.EnableLocalSearch=true; local.LocalSearchAtRoot=true;
local.LocalSearchMaxCalls=3; local.LocalSearchMaxTimePerCall=2;
local.LocalSearchMaxTotalTimeFraction=1; local.LocalSearchCooldownNodes=1;
local.LocalSearchStagnationNodes=0; local.LocalSearchPoorUBTrigger=false;
r=localRun("constant",0.2,1.1,local);
assert(r.stats.localSearchRootCalls==1 && r.stats.localSearchPeriodicCalls==2);
assert(NIST_TOLERANCE_TEST_CASE.localCalls==3);
fprintf('Fine-screening orchestration tests passed: M641 to 3, m1 grid cap, screening, unavailable bounds and local search.\n');
end

function result=localRun(mode,lower,screen,options)
global NIST_TOLERANCE_TEST_CASE
NIST_TOLERANCE_TEST_CASE=struct('mode',mode,'lower',lower,'screen',screen, ...
    'options',options,'localCalls',0,'fineGrids',zeros(1,0),'screenGrids',zeros(1,0));
problem=struct('id',"test",'instance_id',"test_study", ...
    'nLinear',1,'nNonlinear',1,'nResidual',1,'domainProfile',"study", ...
    'zLower',0,'zUpper',1,'objectiveScale',1,'residualScale',1,'certifiedRSS',1, ...
    'initial',struct('linear',1,'z',0.5));
result=solve_nist_bnb_tolerance_study(problem,options);
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
