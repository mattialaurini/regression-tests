function [records,audit,priorAttempts]=nist_diagonal_reuse(files,problems,bases,core,matlabSources,plan)
%NIST_INSTANCE_REUSE Import only source-, problem- and option-compatible runs.
records=cell(0,1); audit=table(); priorAttempts=cell(0,1); seen=containers.Map('KeyType','char','ValueType','logical');
for file=reshape(string(files),1,[])
    before=numel(records); message="";
    try
        s=load(file);
        if isfield(s,'manifest')
            if isfield(s,'results'), visitTable(s.results,localField(s,'stats',{}),s.manifest,"",file); end
            visitManifest(s.manifest,file,0);
        end
    catch exception
        message=string(exception.message);
    end
    audit=[audit;table(file,numel(records)-before,message, ...
        'VariableNames',{'source_file','imported','note'})]; %#ok<AGROW>
end

    function visitManifest(m,file,depth)
        if depth>16 || ~isstruct(m), return, end
        if isfield(m,'referenceResults') && istable(m.referenceResults)
            context=localField(m,'referenceManifest',m);
            visitTable(m.referenceResults,localField(m,'referenceStats',{}),context,"",file);
        end
        if isfield(m,'baselineResults') && istable(m.baselineResults)
            visitTable(m.baselineResults,{},m,"with_m2",file);
        end
        if isfield(m,'baselineRuns') && iscell(m.baselineRuns)
            for j=1:numel(m.baselineRuns)
                run=m.baselineRuns{j};
                if isstruct(run) && all(isfield(run,{'row','options','result'}))
                    one(run.row,localField(run.result,'stats',struct()),m,"with_m2",file,run.options);
                end
            end
        end
        if isfield(m,'referenceManifest'), visitManifest(m.referenceManifest,file,depth+1); end
    end

    function visitTable(t,stats,m,hint,file)
        if ~istable(t) || ~all(ismember({'problem_id','epsilon_relative','error_message'},t.Properties.VariableNames))
            return
        end
        for j=1:height(t)
            st=struct();
            if iscell(stats) && numel(stats)==height(t), st=stats{j}; end
            one(t(j,:),st,m,hint,file,[]);
        end
    end

    function one(raw,st,m,hint,file,explicitOptions)
        % Reuse is an optimization only. Missing metadata schedules a fresh
        % trial; it can never authorize a certificate or silently change options.
        try
            if ~ismember(raw.epsilon_relative,[0.1,0.01,0.001]) || strlength(raw.error_message)>0, return, end
            if ismember('solver',raw.Properties.VariableNames) && raw.solver~="BnB", return, end
            if ismember('phase',raw.Properties.VariableNames) && raw.phase~="tune", return, end
            id=string(raw.problem_id); k=find(cellfun(@(p)string(p.id)==id,problems));
            if ~isscalar(k), return, end
            sourceProblems=localField(m,'problems',{});
            if ~iscell(sourceProblems), return, end
            j=find(cellfun(@(p)string(p.id)==id,sourceProblems));
            if ~isscalar(j) || ~isequaln(validate_nist_problem(sourceProblems{j}),problems{k}), return, end
            variant=hint;
            if ismember('variant',raw.Properties.VariableNames), variant=string(raw.variant); end
            if strlength(variant)==0
                version=string(localField(m,'version',""));
                if contains(version,"m1-grid-tuning") || contains(version,"m2-removal")
                    variant="without_m2";
                else, variant="with_m2"; end
            end
            if ~any(variant==["with_m2","without_m2"]), return, end
            if ~sourceMatches(m,variant), return, end
            o=explicitOptions;
            if isempty(o)
                if isfield(m,'options'), o=m.options{j};
                elseif isfield(m,'originalOptions'), o=m.originalOptions{j};
                elseif isfield(m,'baseOptions'), o=m.baseOptions{j};
                elseif isfield(m,'baselineRuns') && isfield(m,'schedule')
                    c=find(m.schedule.epsilon_relative==raw.epsilon_relative);
                    if ~isscalar(c), return, end
                    o=m.baselineRuns{c,j}.options;
                else, return, end
                if ismember('initial_M',raw.Properties.VariableNames)
                    o.InitialApproxGridPoints=raw.initial_M; o.MaximumFineGridPoints=raw.m1_grid_cap;
                    if variant=="with_m2" && ismember('coarse_N',raw.Properties.VariableNames), o.CoarseGridPoints=raw.coarse_N; end
                end
                if all(ismember({'M','N','Nhat'},raw.Properties.VariableNames))
                    o.InitialApproxGridPoints=raw.M; o.MaximumFineGridPoints=raw.Nhat;
                    if variant=="with_m2", o.CoarseGridPoints=raw.N; end
                end
                o.EpsilonRelative=raw.epsilon_relative; o.EpsilonAbsolute=raw.epsilon_absolute;
                o.TimeLimit=raw.time_limit;
            end
            M=o.InitialApproxGridPoints; H=o.MaximumFineGridPoints; N=o.CoarseGridPoints;
            if variant=="without_m2"
                if M~=H || ~isfield(st,'m2Evaluations') || st.m2Evaluations~=0 || st.nodesReinserted~=0, return, end
                N=NaN;
            elseif N~=2 || M~=H
                return
            end
            if M~=H || ~ismember(M,[3,5,9,17]), return, end
            er=raw.epsilon_relative;
            selected=find(plan.problem_id==id & plan.variant==variant & ...
                plan.M==M & plan.epsilon_relative==er);
            if ~isscalar(selected), return, end
            spec=plan(selected,:); task=spec.task_id;
            expected=nist_diagonal_options(bases{k},spec);
            % CoarseGridPoints is never inspected by this no-m2 solver.
            if variant=="without_m2", o.CoarseGridPoints=2; end
            if isempty(fieldnames(st)), return, end
            U=raw.objective_scaled; L=raw.dual_bound_scaled;
            gap=isfinite(U) && isfinite(L) && U>=0 && L<=U && L>=(1-er)*U-1e-10;
            if logical(raw.certified)~=gap, return, end
            budgetEligible=nist_diagonal_reuse_eligible(raw,o,expected);
            priorAttempts{end+1,1}=struct('task_id',task,'source_file',file, ...
                'rawRow',raw,'stats',st,'options',o, ...
                'metrics',nist_diagonal_metrics(raw,st,o),'budget_eligible',budgetEligible);
            if ~budgetEligible, return, end
            if isKey(seen,char(task)), return, end
            record=nist_diagonal_record(raw,st,o,spec,file);
            records{end+1,1}=record; seen(char(task))=true;
        catch
            % Metadata-incomplete historical runs are not reused. Their
            % planned configurations remain pending and are measured anew.
        end
    end

    function yes=sourceMatches(m,variant)
        yes=false;
        if ~isfield(m,'matlabSources') || ~isequaln(m.matlabSources,matlabSources), return, end
        names=["nist_branch_candidates.m","compute_nist_m1_bound_stable.m", ...
            "compute_nist_m2_screening_stable.m","solve_nist_outer_relaxation_qp_stable.m", ...
            "nist_solve_outer_qp_verified.m","nist_outer_qp_primal_upper.m", ...
            "nist_verify_qp_dual.m","nist_tolerance_prune_masks.m","nist_tolerance_result_row.m"];
        if variant=="with_m2", names=[names,"solve_nist_bnb_tolerance_study.m"];
        else, names=[names,"solve_nist_bnb_m2_removed.m"]; end
        if startsWith(string(localField(m,'version',"")),"diagonal-grid-v250") && variant=="with_m2"
            names=[names,"solve_nist_bnb_diagonal_counted.m"];
        end
        for name=names
            text=findSource(m,name,0); index=find(core.names==name);
            if ~isscalar(index) || strlength(string(text))==0 || ~strcmp(text,core.contents{index}), return, end
        end
        yes=true;
    end
end

function text=findSource(m,name,depth)
text=''; if depth>16 || ~isstruct(m), return, end
for field=["experimentSources","extensionSources","ablationSources"]
    if ~isfield(m,field), continue, end
    s=m.(field);
    index=find(s.names==name);
    if isempty(index), index=find(s.names=="reference/"+name+".txt"); end
    if isscalar(index), text=s.contents{index}; return, end
end
if isfield(m,'referenceManifest'), text=findSource(m.referenceManifest,name,depth+1); end
end
function value=localField(s,name,default)
if isstruct(s)&&isfield(s,name), value=s.(name); else, value=default; end
end
