function [results,summary,progress]=nist_diagonal_write(records,plan,manifest,priorAttempts,folder)
%NIST_DIAGONAL_WRITE Atomic partial reports after every completed trial.
results=table(); stats=cell(numel(records),1); rawRows=stats; options=stats;
for j=1:numel(records)
    results=[results;records{j}.row]; %#ok<AGROW>
    stats{j}=records{j}.stats; rawRows{j}=records{j}.rawRow; options{j}=records{j}.options;
end
if ~isempty(results), assert(numel(unique(results.task_id))==height(results)); end
progress=plan;
progress.execution_state=repmat("pending",height(plan),1);
progress.execution_state(plan.resource_blocked)="resource_limit";
progress.status=progress.execution_state; progress.certified=false(height(plan),1);
progress.reused=false(height(plan),1); progress.source_file=strings(height(plan),1);
progress.node_cost_count_source=strings(height(plan),1);
numeric=["wall_time","time_reported","time_limit","objective_scaled","dual_bound_scaled", ...
    "gap_absolute_scaled","gap_relative_scaled","gap_relative_percent", ...
    "nodes_visited","nodes_cheap","nodes_expensive","pruned_tolerance_only"];
for name=numeric, progress.(name)=nan(height(plan),1); end
progress.node_cost_counts_available=false(height(plan),1);
for j=1:height(results)
    i=find(progress.task_id==results.task_id(j)); assert(isscalar(i));
    progress.execution_state(i)="completed";
    for name=["status","certified","reused","source_file","node_cost_count_source", ...
            "node_cost_counts_available",numeric]
        progress.(name)(i)=results.(name)(j);
    end
end
summary=table();
for c=unique(plan.case_id,'stable').'
    for variant=["with_m2","without_m2"]
        q=progress(progress.case_id==c & progress.variant==variant,:);
        completed=q.execution_state=="completed";
        row=table(c,variant,height(q),sum(completed),sum(q.reused),sum(q.certified), ...
            sum(q.execution_state=="pending"),sum(q.execution_state=="resource_limit"), ...
            sum(completed & q.status=="timelimit"), ...
            ~any(q.execution_state=="pending"),all(completed), ...
            'VariableNames',{'case_id','variant','planned','completed','reused','certified', ...
            'pending','resource_skipped','timeouts','all_executable_completed','all_planned_completed'});
        summary=[summary;row]; %#ok<AGROW>
    end
end
comparison=table();
for c=unique(plan.case_id,'stable').'
    for id=unique(plan.problem_id,'stable').'
        for M=[3,5,9,17]
            q=progress(progress.case_id==c & progress.problem_id==id & progress.M==M,:);
            a=q(q.variant=="with_m2",:); b=q(q.variant=="without_m2",:);
            ratio=NaN;
            if a.certified&&b.certified&&a.wall_time>0, ratio=b.wall_time/a.wall_time; end
            row=table(c,id,M,a.execution_state,b.execution_state,a.certified,b.certified, ...
                a.time_reported,b.time_reported,a.gap_relative_percent,b.gap_relative_percent, ...
                a.nodes_visited,b.nodes_visited,a.nodes_cheap,a.nodes_expensive, ...
                a.pruned_tolerance_only,b.pruned_tolerance_only,ratio, ...
                'VariableNames',{'case_id','problem_id','M','with_m2_state','without_m2_state', ...
                'with_m2_certified','without_m2_certified','with_m2_time','without_m2_time', ...
                'with_m2_gap_percent','without_m2_gap_percent','with_m2_nodes','without_m2_nodes', ...
                'with_m2_cheap','with_m2_expensive','with_m2_Ntol','without_m2_Ntol', ...
                'certified_time_ratio_without_over_with'});
            comparison=[comparison;row]; %#ok<AGROW>
        end
    end
end
csv(results,fullfile(folder,'all_trials.csv'));
csv(progress,fullfile(folder,'study_progress.csv'));
csv(summary,fullfile(folder,'summary.csv'));
csv(comparison,fullfile(folder,'grid_comparison.csv'));
csv(plan,fullfile(folder,'search_plan.csv'));
csv(progress(progress.execution_state=="resource_limit",:),fullfile(folder,'resource_skips.csv'));
historical=table();
for j=1:numel(priorAttempts)
    h=priorAttempts{j}; r=h.rawRow; m=h.metrics;
    one=table(h.task_id,string(h.source_file),r.time_limit,r.wall_time,r.certified, ...
        r.objective_scaled,r.dual_bound_scaled,m.gap_absolute_scaled,m.gap_relative_percent, ...
        r.nodes_visited,m.nodes_cheap,m.nodes_expensive,h.budget_eligible, ...
        'VariableNames',{'task_id','source_file','time_limit','wall_time','certified', ...
        'objective_scaled','dual_bound_scaled','gap_absolute_scaled','gap_relative_percent', ...
        'nodes_visited','nodes_cheap','nodes_expensive','budget_eligible'});
    historical=[historical;one]; %#ok<AGROW>
end
csv(historical,fullfile(folder,'historical_attempts.csv'));
payload=struct('results',results,'summary',summary,'comparison',comparison,'progress',progress, ...
    'stats',{stats},'rawRows',{rawRows},'options',{options},'manifest',manifest, ...
    'priorAttempts',{priorAttempts},'historical',historical);
file=fullfile(folder,'nist_diagonal_grid_results.mat'); temp=string(file)+'.tmp.mat';
save(temp,'-struct','payload','-v7.3'); movefile(temp,file,'f');
end
function csv(t,file)
if width(t)==0, return, end
temp=string(file)+'.tmp.csv'; writetable(t,temp); movefile(temp,file,'f');
end
