function test_nist_diagonal_grid(problems,bases,sources,matlabSources)
%TEST_NIST_DIAGONAL_GRID Planning, budget reuse, gaps, counters and exports.
plan=nist_diagonal_plan(problems,8);
assert(height(plan)==192 && all(plan.M==plan.Nhat));
assert(all(plan.N(plan.variant=="with_m2")==2));
assert(all(isnan(plan.N(plan.variant=="without_m2"))));
assert(all(plan.target_time_limit(plan.epsilon_relative==0.01)==1000));
gauss=plan.problem_id=="Gauss1";
assert(all(plan.resource_blocked(gauss & plan.M==17)));
assert(all(plan.minimum_dense_GiB(gauss & plan.M==17)>60));
assert(all(plan.resource_blocked(gauss & plan.M==9)));
assert(~any(plan.resource_blocked(~gauss)));

p=problems{1}; spec=plan(plan.problem_id==string(p.id)&plan.variant=="with_m2"& ...
    plan.epsilon_relative==0.01&plan.M==3,:);
expected=nist_diagonal_options(bases{1},spec); old=expected; old.TimeLimit=600;
raw=fixture(p); st=struct('nodesReinserted',0,'precisionFallbackM1Evaluations',0, ...
    'm1Evaluations',1,'m2Evaluations',1,'m1SkippedByScreening',0);
assert(nist_diagonal_reuse_eligible(raw,old,expected));
record=nist_diagonal_record(raw,st,old,spec,"fixture");
assert(record.row.time_limit==600&&record.row.target_time_limit==1000&&record.row.reuse_shorter_certified_budget);
assert(record.row.nodes_expensive==1&&record.row.nodes_cheap==0);
assert(abs(record.row.gap_relative_percent-0.5)<1e-10);
failed=raw; failed.status="timelimit"; failed.certified=false;
failed.dual_bound_scaled=0; failed.wall_time=600.1;
assert(~nist_diagonal_reuse_eligible(failed,old,expected),'600-s failures must be rerun.');
failed.time_limit=1000; old1000=expected;
assert(nist_diagonal_reuse_eligible(failed,old1000,expected),'A completed 1000-s timeout must not loop on resume.');
bad=old; bad.QPOptions.FeasibilityTolerance=1e-4;
assert(~nist_diagonal_reuse_eligible(raw,bad,expected),'QP options must remain identical.');

test=st; test.precisionFallbackM1Evaluations=1;
x=nist_diagonal_metrics(raw,test,old);
assert(~x.node_cost_counts_available&&isnan(x.nodes_expensive),'Ambiguous historical counts must not be invented.');
test.nodesWithM1=1; test.nodesWithoutM1=0;
x=nist_diagonal_metrics(raw,test,old); assert(x.nodes_expensive==1);
zero=raw; zero.objective_scaled=0; zero.dual_bound_scaled=0;
x=nist_diagonal_metrics(zero,st,old); assert(x.gap_relative_scaled==0&&x.gap_absolute_scaled==0);
gap=raw; gap.dual_bound_scaled=0;
x=nist_diagonal_metrics(gap,st,old); assert(x.gap_relative_percent==100);
bad=raw; bad.dual_bound_scaled=2;
x=nist_diagonal_metrics(bad,st,old); assert(~x.bounds_consistent&&isnan(x.gap_relative_scaled));

folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
results=raw; results.M=3; results.N=2; results.Nhat=3;
results.variant="with_m2"; results.phase="tune";
stats={st};
manifest=struct('version',"instance-tuning-test",'problems',{{p}}, ...
    'baseOptions',{{bases{1}}},'matlabSources',matlabSources,'experimentSources',sources);
file=string(fullfile(folder,'reuse.mat')); save(file,'results','stats','manifest');
[imported,~,history]=nist_diagonal_reuse(file,problems,bases,sources,matlabSources,plan);
assert(numel(imported)==1&&numel(history)==1,'The per-instance M/N/Nhat schema must import.');
assert(imported{1}.row.time_limit==600&&imported{1}.row.target_time_limit==1000);
results=failed; results.M=3; results.N=2; results.Nhat=3;
results.variant="with_m2"; results.phase="tune"; results.time_limit=600;
save(file,'results','stats','manifest');
[imported,~,history]=nist_diagonal_reuse(file,problems,bases,sources,matlabSources,plan);
assert(isempty(imported)&&numel(history)==1,'Old failed attempts must be retained but not reused.');

% Confirmation timings are a separate experiment, not replacement tuning data.
results=raw; results.M=3; results.N=2; results.Nhat=3;
results.variant="with_m2"; results.phase="confirm";
save(file,'results','stats','manifest');
[imported,~]=nist_diagonal_reuse(file,problems,bases,sources,matlabSources,plan);
assert(isempty(imported));

% The strict grid-5 study has initial_M/m1_grid_cap but no coarse_N.
strict=plan(plan.problem_id==string(p.id)&plan.variant=="with_m2"&plan.M==5& ...
    plan.epsilon_relative==0.001,:);
o=nist_diagonal_options(bases{1},strict); results=raw;
results.epsilon_relative=0.001; results.time_limit=3600; results.case_id="rel_0p001";
results.dual_bound_scaled=0.9995; results.initial_M=5; results.m1_grid_cap=5;
results.variant="with_m2";
manifest=rmfield(manifest,'baseOptions'); manifest.options={o};
save(file,'results','stats','manifest');
[imported,~]=nist_diagonal_reuse(file,problems,bases,sources,matlabSources,plan);
assert(numel(imported)==1&&imported{1}.row.M==5,'Strict .001 grid-5 results must be reused.');

[out,summary,progress]=nist_diagonal_write({record},plan,manifest,{},folder);
assert(height(out)==1&&height(progress)==192);
assert(sum(summary.completed)+sum(summary.pending)+sum(summary.resource_skipped)==192);
assert(~all(summary.all_planned_completed));
assert(isfile(fullfile(folder,'nist_diagonal_grid_results.mat')));
fprintf('Diagonal study self-tests passed: 192 cases, memory skips, 600/1000-s reuse, gaps, exact node partitions and partial reports.\n');
end
function t=fixture(p)
r=struct('instance_id',string(p.instance_id),'problem_id',string(p.id),'case_id',"rel_0p01", ...
    'solver',"BnB",'status',"paper_gap",'certified',true,'objective_scaled',1, ...
    'dual_bound_scaled',0.995,'rss',p.objectiveScale,'dual_bound_rss',0.995*p.objectiveScale, ...
    'objective_scale',p.objectiveScale,'wall_time',2,'time_limit',600,'par2',2, ...
    'epsilon_relative',0.01,'epsilon_absolute',1e-10,'nodes_visited',1, ...
    'pruned_tolerance_only',1,'pruned_at_incumbent',0,'pruned_total',1, ...
    'm1_evaluations',1,'m2_evaluations',1,'m1_qp_retries',0,'m2_qp_retries',0, ...
    'm1_unavailable_bounds',0,'m2_unavailable_evaluations',0,'error_message',"");
t=struct2table(r);
end
