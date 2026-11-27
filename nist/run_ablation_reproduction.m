function [results,summary]=run_ablation_reproduction(state,output,args,ids)
%RUN_ABLATION_REPRODUCTION Fresh trials with published per-trial budgets.
% Numerical algorithms and their stopping criteria are the frozen 2.5.0 code.
plan=nist_diagonal_plan(state.problems,8);
executionLimits=plan.target_time_limit;
if args.BudgetMode=="published"
    published=readtable(fullfile(state.root,'reference','ablation_results.csv'),TextType='string');
    for k=1:height(plan)
        j=find(published.task_id==plan.task_id(k));
        if ~isempty(j)
            assert(isscalar(j)); executionLimits(k)=published.time_limit(j);
        else
            assert(plan.resource_blocked(k),'Published observation is missing.');
        end
    end
end
if ~isfolder(output), mkdir(output); end
manifest=state.manifest;
manifest.version="reproduction-ablation-1.0.0";
manifest.plan=plan; manifest.executionLimits=executionLimits; manifest.budgetMode=args.BudgetMode;
manifest.driverSources=driverSources(state.root);
configFile=fullfile(output,'reproduction_config.mat');
if isfile(configFile)
    old=load(configFile,'manifest');
    assert(isequaln(old.manifest,manifest),'Inputs/source/budgets changed: use a new output directory.');
else
    atomicSave(configFile,struct('manifest',manifest));
end
records=cell(0,1); completed=strings(0,1);
files=dir(fullfile(output,'runs','*.mat'));
for j=1:numel(files)
    if endsWith(files(j).name,'.tmp.mat'), continue, end
    d=load(fullfile(files(j).folder,files(j).name),'record');
    assert(strlength(d.record.row.error_message)==0,'An error checkpoint exists; inspect it before resuming.');
    records{end+1,1}=d.record; completed(end+1,1)=d.record.row.task_id; %#ok<AGROW>
end
assert(numel(unique(completed))==numel(completed),'Duplicate checkpoint.');
[results,summary,progress]=nist_diagonal_write(records,plan,manifest,{},output);
cases=["rel_0p1","rel_0p01","rel_0p001"];
chosen=ismember(plan.case_id,cases(args.CaseIndices))&ismember(plan.problem_id,ids);
todo=find(chosen&~plan.resource_blocked&~ismember(plan.task_id,completed));
preview=plan(chosen,:); preview.execution_time_limit=executionLimits(chosen);
writetable(preview,fullfile(output,'execution_plan.csv'));
if args.DryRun
    disp(preview(:,{'problem_id','case_id','variant','M','target_time_limit','execution_time_limit','resource_blocked'}));
    fprintf('Dry run: %d pending solves; %d memory exclusions. No optimizer called.\n', ...
        numel(todo),sum(chosen&plan.resource_blocked)); return
end
if isempty(todo), fprintf('All selected trials are already complete.\n'); return, end
if args.RunSelfTests
    test_nist_diagonal_grid(state.problems,state.baseOptions, ...
        manifest.experimentSources,manifest.matlabSources);
    test_nist_diagonal_node_counts(); test_nist_fine_screening_v250();
    test_nist_m2_removal_v250(); test_nist_qp_recovery();
    report=check_nist_qp_recovery(state.baseOptions{1}.QPOptions);
    atomicSave(fullfile(output,'qp_recovery_preflight.mat'),struct('report',report));
    assert(report.passed,'Gurobi QP recovery preflight failed.');
end
assert(~isempty(which('fmincon')),'Optimization Toolbox (fmincon) is required.');
warmup_nist_solvers(state.baseOptions{1}.QPOptions);
for j=reshape(todo,1,[])
    spec=plan(j,:); k=find(cellfun(@(p)string(p.id)==spec.problem_id,state.problems));
    p=state.problems{k}; options=nist_diagonal_options(state.baseOptions{k},spec);
    options.TimeLimit=executionLimits(j);
    fprintf('\nRunning %s / %s / %s: M=Nhat=%g, limit=%g s\n', ...
        spec.case_id,spec.problem_id,spec.variant,spec.M,options.TimeLimit);
    message="";
    try
        if spec.variant=="with_m2", result=solve_nist_bnb_diagonal_counted(p,options);
        else, result=solve_nist_bnb_m2_removed(p,options); end
    catch exception
        message=string(getReport(exception,'extended','hyperlinks','off'));
        result=struct('solver',"BnB",'status',"error",'problem_id',string(p.id), ...
            'instance_id',string(p.instance_id));
    end
    raw=nist_tolerance_result_row(result,p,spec.case_id,spec.epsilon_relative,1e-10,options.TimeLimit,message);
    stats=struct(); if isfield(result,'stats'), stats=result.stats; end
    record=nist_diagonal_record(raw,stats,options,spec,"new_run");
    if strlength(message)==0
        assert(stats.nodesReinserted==0);
        assert(record.row.nodes_cheap+record.row.nodes_expensive==record.row.nodes_visited);
        if spec.variant=="without_m2", assert(stats.m2Evaluations==0); end
    end
    atomicSave(fullfile(output,'runs',spec.task_id+'.mat'),struct('record',record,'result',result));
    records{end+1,1}=record; %#ok<AGROW>
    [results,summary]=nist_diagonal_write(records,plan,manifest,{},output);
    r=record.row;
    fprintf('U=%.9g L=%.9g (RSS/s), t=%.3f s, nodes=%.0f (%g cheap, %g costly), gap=%.5g%%, Ntol=%g\n', ...
        r.objective_scaled,r.dual_bound_scaled,r.wall_time,r.nodes_visited,r.nodes_cheap, ...
        r.nodes_expensive,r.gap_relative_percent,r.pruned_tolerance_only);
    if strlength(message)>0
        error('NISTReproduction:Failure','Error checkpoint saved. This is not a timeout.');
    end
end
disp(summary);
end
function atomicSave(file,payload)
if ~isfolder(fileparts(file)), mkdir(fileparts(file)); end
temporary=file+'.tmp.mat'; save(temporary,'-struct','payload','-v7.3'); movefile(temporary,file,'f');
end
function s=driverSources(root)
files=dir(fullfile(root,'*.m')); names=string({files.name})'; contents=cell(numel(files),1);
for j=1:numel(files), contents{j}=fileread(fullfile(root,files(j).name)); end
s=struct('names',names,'contents',{contents});
end
