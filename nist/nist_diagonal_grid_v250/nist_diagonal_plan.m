function plan=nist_diagonal_plan(problems,maxEstimatedMemoryGiB)
%NIST_DIAGONAL_PLAN Every problem/grid/variant/tolerance; no dimension cap.
plan=table(); ers=[0.1,0.01,0.001]; limits=[100,1000,3600];
cases=["rel_0p1","rel_0p01","rel_0p001"];
for c=1:3
    for j=1:numel(problems)
        p=problems{j};
        for M=[3,5,9,17]
            for variant=["with_m2","without_m2"]
                N=2; if variant=="without_m2", N=NaN; end
                key="m2_M"+M+"_N2_H"+M;
                if variant=="without_m2", key="no_m2_K"+M; end
                task=cases(c)+"_"+lower(string(p.id))+"_"+key;
                columns=size(p.linearVertices,2)*M^p.nNonlinear;
                % Three dense residual/interpolation arrays alone are a
                % hard lower bound, not a prediction of Gurobi memory.
                minGiB=3*8*p.nResidual*columns/2^30;
                % Dense copies, stacked constraints and sparse conversion
                % coexist in the existing backend. A planning estimate only.
                estimateGiB=24*8*p.nResidual*columns/2^30;
                blocked=estimateGiB>maxEstimatedMemoryGiB;
                reason="";
                if blocked
                    reason=sprintf('Memory estimate %.2f GiB exceeds configured %.2f GiB; three dense arrays alone %.2f GiB.', ...
                        estimateGiB,maxEstimatedMemoryGiB,minGiB);
                end
                row=table(task,cases(c),string(p.id),variant,M,N,M,ers(c),1e-10,limits(c), ...
                    key,columns,minGiB,estimateGiB,blocked,string(reason), ...
                    'VariableNames',{'task_id','case_id','problem_id','variant','M','N','Nhat', ...
                    'epsilon_relative','epsilon_absolute','target_time_limit','config_id', ...
                    'm1_columns','minimum_dense_GiB','estimated_memory_GiB','resource_blocked','resource_reason'});
                plan=[plan;row]; %#ok<AGROW>
            end
        end
    end
end
assert(height(plan)==192 && numel(unique(plan.task_id))==192);
end
