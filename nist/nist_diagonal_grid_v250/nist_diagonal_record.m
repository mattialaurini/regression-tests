function record=nist_diagonal_record(raw,stats,options,spec,source)
%NIST_DIAGONAL_RECORD Preserve actual historical budget and time, unmodified.
keep={'instance_id','problem_id','case_id','solver','status','certified', ...
    'objective_scaled','dual_bound_scaled','rss','dual_bound_rss','objective_scale', ...
    'wall_time','time_limit','par2','nodes_visited','pruned_tolerance_only', ...
    'pruned_at_incumbent','pruned_total','m1_evaluations','m2_evaluations', ...
    'm1_qp_retries','m2_qp_retries','m1_unavailable_bounds','m2_unavailable_evaluations','error_message'};
assert(all(ismember(keep,raw.Properties.VariableNames)),'Incomplete result schema.');
row=raw(:,keep);
assert(string(raw.problem_id)==spec.problem_id && raw.epsilon_relative==spec.epsilon_relative && ...
    raw.epsilon_absolute==spec.epsilon_absolute);
metrics=nist_diagonal_metrics(raw,stats,options);
gap=metrics.bounds_consistent&&metrics.gap_margin_scaled>=0;
if strlength(raw.error_message)==0, assert(logical(raw.certified)==gap,'Gap certificate mismatch.'); end
row.variant=spec.variant; row.M=spec.M; row.N=spec.N; row.Nhat=spec.Nhat;
row.epsilon_relative=spec.epsilon_relative; row.epsilon_absolute=spec.epsilon_absolute;
row.config_id=spec.config_id; row.task_id=spec.task_id;
row.target_time_limit=spec.target_time_limit;
row.source_file=string(source); row.reused=string(source)~="new_run";
row.reuse_shorter_certified_budget=row.reused && raw.time_limit<spec.target_time_limit;
row.m1_columns=spec.m1_columns;
row.minimum_dense_GiB=spec.minimum_dense_GiB;
row.estimated_memory_GiB=spec.estimated_memory_GiB;
row.time_reported=raw.wall_time;
if ~raw.certified && string(raw.status)=="timelimit", row.time_reported=raw.time_limit; end
row.target_par2=raw.wall_time;
if ~raw.certified || ~isfinite(raw.wall_time), row.target_par2=2*spec.target_time_limit; end
for n=reshape(string(fieldnames(metrics)),1,[]), row.(n)=metrics.(n); end
record=struct('row',row,'rawRow',raw,'stats',stats,'options',options);
end
