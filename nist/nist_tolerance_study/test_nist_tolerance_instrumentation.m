function test_nist_tolerance_instrumentation()
%TEST_NIST_TOLERANCE_INSTRUMENTATION Deterministic orchestration tests, no QPs.
% The actual instrumented solver runs against temporary, isolated test doubles.
% No original package file is changed. Temporary paths are restored on error.
root=fileparts(mfilename('fullpath'));
assert(isequal(nist_tolerance_solver_selection("bnb"),"BnB"));
assert(isequal(nist_tolerance_solver_selection("scip"),"SCIP"));
assert(isequal(nist_tolerance_solver_selection(["scip","bnb"]),["SCIP","BnB"]));
rejected=false;
try
    nist_tolerance_solver_selection(["BnB","bnb"]);
catch exception
    rejected=strcmp(exception.identifier,'NISTStudy:SolverSelection');
end
assert(rejected,'Duplicate solvers must be rejected.');
previousPath=path; previousFolder=pwd; temporary=tempname; mkdir(temporary);
cleanup=onCleanup(@()localCleanup(previousPath,previousFolder,temporary)); %#ok<NASGU>
files=dir(fullfile(root,'test_fixtures','*.m.txt'));
for i=1:numel(files)
    copyfile(fullfile(files(i).folder,files(i).name), ...
        fullfile(temporary,erase(files(i).name,'.txt')));
end
addpath(root); addpath(temporary,'-begin');
% MATLAB searches the current folder before entries added to the path.
% Changing into the private folder guarantees that the doubles take priority.
cd(temporary); localClear();
global NIST_TOLERANCE_TEST_CASE
source=fileread(fullfile(root,'reviewed_source','solve_nist_bnb.m.txt'));
tokens=regexp(source,'options\.([A-Za-z]\w*)','tokens');
names=unique(string([tokens{:}])); opts=struct();
for name=names, opts.(name)=0; end
opts.TimeLimit=Inf; opts.NodeLimit=Inf;
opts.EpsilonRelative=0.1; opts.EpsilonAbsolute=1e-10;
opts.CoarseGridPoints=2; opts.InitialApproxGridPoints=2;
opts.MaximumFineGridPoints=2; opts.StrongBranchingTieBreak=false;
opts.BranchingRule="normalized_longest_edge"; opts.TieTolerance=1e-12;
opts.RecordRelativeTargets=[0.1,0.01,0.001]; opts.Verbosity=0;
opts.ProgressIntervalSeconds=Inf; opts.QPOptions=struct('TimeLimit',Inf);
opts.EnableLocalSearch=false; opts.LocalSearchAtRoot=false;

% Regression: a tiny but longer side must not tie with a zero-width side.
bp=struct('zLower',[0;0],'zUpper',[1;1]);
br=struct('zLower',[1;200],'zUpper',[1;200+300*2^-40]);
[candidates,midpoints]=nist_branch_candidates(br,bp,opts);
assert(isequal(candidates,2));
assert(br.zLower(2)<midpoints(2)&&midpoints(2)<br.zUpper(2));
br=struct('zLower',[0;0],'zUpper',[1e-14;1e-13]);
assert(isequal(nist_branch_candidates(br,bp,opts),2));
% Ordinary near ties keep their deterministic first-coordinate tie break.
br=struct('zLower',[0;0],'zUpper',[1;1+5e-13]);
assert(isequal(nist_branch_candidates(br,bp,opts),[1;2]));
% Adjacent floating-point endpoints have no representable interior midpoint.
br=struct('zLower',1,'zUpper',1+eps(1));
assert(isempty(nist_branch_candidates(br,struct('zLower',1,'zUpper',1+eps(1)),opts)));

tau=0.9-1e-10;
[relaxed,ordinary,eligible]=nist_tolerance_prune_masks( ...
    [tau-eps(1),tau,0.95,1,1.1,NaN],1,0.1,1e-10);
assert(isequal(relaxed,[false,true,true,false,false,false]));
assert(isequal(ordinary,[false,false,false,true,true,false]));
assert(isequal(eligible,relaxed|ordinary));

r=localRun("constant",0.95,1.1,opts,1);
assert(r.nodes==1&&r.stats.prunedToleranceOnlyExplicit==1);
assert(r.stats.prunedAtIncumbent==0&&r.stats.prunedTotal==1);

r=localRun("constant",1,1.1,opts,1);
assert(r.stats.prunedAtIncumbentExplicit==1&&r.stats.prunedToleranceOnly==0);

% Equality belongs to the inclusive pruning criterion, not to branching.
r=localRun("constant",tau,tau,opts,1);
assert(r.stats.prunedToleranceOnlyExplicit==1&&r.stats.nodesBranched==0);

% A region extracted three times after grid reductions is one visited node.
reinsert=opts; reinsert.InitialApproxGridPoints=5;
r=localRun("constant",0.95,1.1,reinsert,1);
assert(r.nodes==1&&r.queue_extractions==3&&r.stats.nodesReinserted==2);

% Precision fallback may prune only using a valid m1 certificate.
precision=struct('zLower',1,'zUpper',1+eps(1));
r=localRun("constant",0.95,0,opts,1,precision);
assert(r.paper_certificate_satisfied && r.nodes==1 && r.stats.nodesCreated==1);
assert(r.stats.precisionFallbackM1Evaluations==1 && r.stats.nodesBranched==0);
assert(r.stats.m1SkippedByScreening==0 && r.stats.prunedToleranceOnly==1);
% Otherwise keep the region unresolved; precision is never an optimality proof.
r=localRun("constant",0.2,0,opts,1,precision);
assert(r.status=="numerical_resolution_limit" && ~r.paper_certificate_satisfied);
assert(r.stats.precisionLimitedNodes==1 && r.stats.prunedTotal==0);
assert(r.stats.nodesBranched==0 && r.stats.nodesCreated==1);

% An unavailable m2 must force m1, never act as a zero screening value.
r=localRun("constant",0.95,Inf,opts,1);
assert(r.paper_certificate_satisfied && r.stats.m1Evaluations==1);
assert(r.stats.m2UnavailableEvaluations==1 && r.stats.m1SkippedByScreening==0);

% Initial certificate closes an existing, unvisited root region.
r=localRun("constant",0.95,1.1,opts,1e-12);
assert(r.nodes==0&&r.stats.closedByGlobalGap==1);
assert(r.stats.prunedToleranceOnlyGlobal==1&&r.stats.nodesPruned==0);

% Root LB=.24; a child improves UB to .25 and is pruned at LB=.25.
% The remaining child's inherited LB=.24 is closed by the global gap test.
r=localRun("global_queue",0.24,1.1,opts,1);
assert(r.nodes==2&&r.stats.nodesCreated==3);
assert(r.stats.prunedAtIncumbentExplicit==1);
assert(r.stats.prunedToleranceOnlyGlobal==1&&r.stats.prunedTotal==2);

% Unresolved queued nodes at a resource limit must not be counted as pruned.
limited=opts; limited.NodeLimit=1;
r=localRun("global_queue",0.24,1.1,limited,1);
assert(~r.paper_certificate_satisfied&&r.stats.prunedTotal==0);
assert(r.stats.closedByGlobalGap==0&&r.nodes==1);

% Historical classification uses the old incumbent, not a later improved UB.
[oldClass,~]=nist_tolerance_prune_masks(0.95,1,0.1,1e-10);
[newClass,~]=nist_tolerance_prune_masks(0.95,0.9,0.1,1e-10);
assert(oldClass&&~newClass);

% Regression check: manuscript tables must use scaled values and preserve
% unavailable SCIP pruning counters as NaN, rather than reporting zero.
folder=fullfile(temporary,'outputs'); mkdir(folder);
fixture=struct('case_id',"rel_0p1",'epsilon_relative',0.1,'time_limit',60, ...
    'solver',"BnB",'certified',false,'wall_time',60,'par2',120, ...
    'nodes_visited',17,'pruned_tolerance_only',3,'pruned_at_incumbent',1, ...
    'objective_scaled',0.00223,'dual_bound_scaled',0,'problem_id',"Lanczos2");
t=struct2table(fixture); other=t; other.solver="SCIP";
other.pruned_tolerance_only=NaN; other.pruned_at_incumbent=NaN;
summary=write_nist_tolerance_outputs([t;other],string(folder));
assert(isnan(summary.total_pruned_tolerance_only(summary.solver=="SCIP")));
latex=fileread(fullfile(folder,'tables','nist_rel_0p1.tex'));
assert(contains(latex,'0.00223'));
assert(~contains(latex,'10^{-11}'));
assert(count(string(latex),'\begin{table}')==2);
fprintf('NIST instrumentation self-tests passed (branching precision regression, 11 scenarios, boundary masks, scaled-table export).\n');
end

function result=localRun(mode,lower,screen,options,initialUpper,region)
if nargin<6, region=struct('zLower',0,'zUpper',1); end
global NIST_TOLERANCE_TEST_CASE
NIST_TOLERANCE_TEST_CASE=struct('mode',mode,'lower',lower, ...
    'screen',screen,'options',options);
problem=struct('id',"test",'instance_id',"test_study", ...
    'nLinear',1,'nNonlinear',1,'nResidual',1,'domainProfile',"study", ...
    'zLower',region.zLower,'zUpper',region.zUpper,'objectiveScale',1,'residualScale',1,'certifiedRSS',1, ...
    'initial',struct('linear',sqrt(initialUpper),'z',region.zLower+(region.zUpper-region.zLower)/2));
result=solve_nist_bnb_tolerance_study(problem,options);
end

function localClear()
clear nist_bnb_options validate_nist_problem evaluate_nist_objective
clear solve_nist_linear_block_qp compute_nist_m1_bound compute_nist_m2_screening
clear compute_nist_m1_bound_stable compute_nist_m2_screening_stable
end

function localCleanup(previousPath,previousFolder,temporary)
cd(previousFolder); path(previousPath); localClear();
clear global NIST_TOLERANCE_TEST_CASE
if isfolder(temporary), rmdir(temporary,'s'); end
end
