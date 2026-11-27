function test_nist_diagonal_node_counts()
% Exercise actual counted solver, including a double-m1 precision fallback.
root=fileparts(mfilename('fullpath')); previousPath=path; previousFolder=pwd;
temporary=tempname; mkdir(temporary);
cleanup=onCleanup(@()restore(previousPath,previousFolder,temporary)); %#ok<NASGU>
files=dir(fullfile(root,'screening_test_fixtures','*.m.txt'));
for j=1:numel(files), copyfile(fullfile(files(j).folder,files(j).name),fullfile(temporary,erase(files(j).name,'.txt'))); end
addpath(root); addpath(temporary,'-begin'); cd(temporary); clearFixtures();
source=fileread(fullfile(root,'solve_nist_bnb_diagonal_counted.m'));
tokens=regexp(source,'options\.([A-Za-z]\w*)','tokens');
names=unique(string([tokens{:}])); o=struct();
for n=names, o.(n)=0; end
o.TimeLimit=Inf; o.NodeLimit=Inf; o.EpsilonRelative=0.1; o.EpsilonAbsolute=1e-10;
o.CoarseGridPoints=2; o.InitialApproxGridPoints=3; o.MaximumFineGridPoints=3;
o.BranchingRule="normalized_longest_edge"; o.StrongBranchingTieBreak=false;
o.TieTolerance=1e-12; o.RecordRelativeTargets=[0.1,0.01,0.001];
o.Verbosity=0; o.ProgressIntervalSeconds=Inf; o.QPOptions=struct('TimeLimit',Inf);
o.EnableLocalSearch=false; o.LocalSearchAtRoot=false;
p=struct('id',"test",'instance_id',"test_study",'nLinear',1,'nNonlinear',1, ...
    'nResidual',1,'domainProfile',"study",'zLower',0,'zUpper',1, ...
    'objectiveScale',1,'residualScale',1,'certifiedRSS',1,'initial',struct('linear',1,'z',0.5));
limited=o; limited.NodeLimit=1;
r=run("constant",0.2,0.5,p,limited);
assert(r.nodes==1&&r.stats.nodesWithM1==0&&r.stats.nodesWithoutM1==1&&r.stats.nodesBranched==1);
r=run("constant",0.95,Inf,p,o);
assert(r.paper_certificate_satisfied&&r.stats.nodesWithM1==1&&r.stats.nodesWithoutM1==0);
r=run("failed_root",0,Inf,p,o);
assert(r.nodes==3&&r.stats.nodesWithM1==3&&r.stats.m1UnavailableBounds==1);
% Adjacent floating-point endpoints cannot be bisected. Both m1 evaluations
% concern the same region, which must be counted as expensive only once.
adjacent=p; adjacent.zLower=1; adjacent.zUpper=1+eps(1); adjacent.initial.z=1;
r=run("constant",0.2,Inf,adjacent,o);
assert(r.stats.m1Evaluations==2&&r.stats.nodesWithM1==1&&r.stats.nodesWithoutM1==0);
assert(r.status=="numerical_resolution_limit");
r=run("constant",0.2,0.5,adjacent,o);
assert(r.stats.m1Evaluations==1&&r.stats.nodesWithM1==1&&r.stats.nodesWithoutM1==0);
fprintf('Unique-node cost tests passed: screening branches, expensive nodes, unavailable bounds and repeated m1 on one region.\n');
end
function r=run(mode,lower,screen,p,o)
global NIST_TOLERANCE_TEST_CASE
NIST_TOLERANCE_TEST_CASE=struct('mode',mode,'lower',lower,'screen',screen, ...
    'options',o,'localCalls',0,'fineGrids',zeros(1,0),'screenGrids',zeros(1,0));
r=solve_nist_bnb_diagonal_counted(p,o);
assert(r.stats.nodesWithM1+r.stats.nodesWithoutM1==r.nodes);
end
function clearFixtures()
clear nist_bnb_options validate_nist_problem evaluate_nist_objective
clear solve_nist_linear_block_qp compute_nist_m1_bound_stable compute_nist_m2_screening_stable
clear solve_nist_reduced_local_search
end
function restore(oldPath,oldFolder,temporary)
cd(oldFolder); path(oldPath); clearFixtures(); clear global NIST_TOLERANCE_TEST_CASE
if isfolder(temporary), rmdir(temporary,'s'); end
end
