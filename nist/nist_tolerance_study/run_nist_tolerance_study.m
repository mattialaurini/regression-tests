function [results,summary] = run_nist_tolerance_study(benchmarkMat,outputDirectory,settings)
%RUN_NIST_TOLERANCE_STUDY Paired, resumable NIST tolerance/budget study.
% Select the ORIGINAL benchmark_results.mat when prompted.
% Requires the original nist_separable_benchmark_v1_0_6 package and solvers.
arguments
    benchmarkMat (1,1) string = ""
    outputDirectory (1,1) string = ""
    settings (1,1) struct = struct()
end
defaults=struct('ProblemIds',strings(1,0),'CaseIndices',1:3,'Solvers',["BnB","SCIP"], ...
    'RunVarPro',false,'RunSelfTests',true,'StopOnError',true);
names=fieldnames(settings);
assert(all(ismember(names,fieldnames(defaults))),'Unknown study setting.');
for k=1:numel(names), defaults.(names{k})=settings.(names{k}); end
settings=defaults;
settings.Solvers=nist_tolerance_solver_selection(settings.Solvers);
runBnB=any(settings.Solvers=="BnB");
runSCIP=any(settings.Solvers=="SCIP");
if strlength(benchmarkMat)==0
    parent=fileparts(fileparts(mfilename('fullpath')));
    candidate=string(fullfile(parent,'nist_separable_benchmark_v1_0_6', ...
        'output','nist_full_v106_t1200','benchmark_results.mat'));
    if isfile(candidate), benchmarkMat=candidate; end
end
if strlength(benchmarkMat)==0
    [name,folder]=uigetfile('*.mat','Select the original benchmark_results.mat');
    if isequal(name,0), results=table(); summary=table(); return, end
    benchmarkMat=string(fullfile(folder,name));
end
stored=load(benchmarkMat,'config');
assert(isfield(stored,'config'),'The selected MAT file must contain config.');
config=stored.config;
config.QPOptions.Threads=1;
if isfield(config,'VarProOptions'), config.VarProOptions.QPOptions.Threads=1; end
% Prefer the package that actually contains the selected original MAT.
% Stored absolute paths may refer to the package before it was moved.
packageRoot=fileparts(fileparts(fileparts(char(benchmarkMat))));
if ~isfolder(fullfile(packageRoot,'matlab'))
    packageRoot=fileparts(fileparts(char(config.OutputDirectory)));
end
matlabRoot=fullfile(packageRoot,'matlab');
if ~isfolder(matlabRoot)
    matlabRoot=uigetdir(pwd,'Select the original benchmark matlab directory');
    assert(~isequal(matlabRoot,0),'No benchmark directory selected.');
    packageRoot=fileparts(matlabRoot);
end
if ~isfolder(config.DataDirectory)
    dataRoot=fullfile(packageRoot,'data','raw');
    if ~isfolder(dataRoot)
        dataRoot=uigetdir(packageRoot,'Select the original data/raw directory');
        assert(~isequal(dataRoot,0),'No data directory selected.');
    end
    config.DataDirectory=string(dataRoot);
end
extensionRoot=fileparts(mfilename('fullpath'));
previousPath=path; restorePath=onCleanup(@()path(previousPath)); %#ok<NASGU>
addpath(genpath(matlabRoot)); addpath(extensionRoot,'-begin');
required=["load_nist_problem","nist_result_to_row","gurobi"];
if runBnB
    required=[required,"solve_nist_bnb","nist_bnb_options"];
    if config.BnBOptions.EnableLocalSearch, required=[required,"fmincon"]; end
end
if runSCIP
    required=[required,"solve_nist_scip","nist_scip_options","import_nist_scip_bridge"];
end
if settings.RunVarPro
    required=[required,"solve_nist_varpro_multistart","fmincon"];
end
for name=required
    assert(~isempty(which(char(name))),'Missing function: %s',name);
end
if runBnB
    reviewed=fileread(fullfile(extensionRoot,'reviewed_source','solve_nist_bnb.m.txt'));
    original=fileread(which('solve_nist_bnb'));
    assert(strcmp(localNormalize(reviewed),localNormalize(original)), ...
        ['The installed BnB source differs from the reviewed v1.0.6 source. ' ...
         'Do not run this private solver copy against a different version.']);
    if settings.RunSelfTests
        test_nist_tolerance_instrumentation();
        test_nist_qp_recovery();
    end
end

problemIds=string(settings.ProblemIds);
if isempty(problemIds), problemIds=string(config.ProblemIds); end
problemIds=reshape(problemIds,1,[]);
assert(numel(unique(lower(problemIds)))==numel(problemIds),'Duplicate problem IDs.');
assert(all(ismember(lower(problemIds),lower(string(config.ProblemIds)))), ...
    'Problem IDs must belong to the original configuration.');
caseIds=["rel_0p1";"rel_0p01";"rel_0p001"];
schedule=table(caseIds,[0.1;0.01;0.001],[100;600;3600], ...
    'VariableNames',{'case_id','epsilon_relative','time_limit'});
indices=settings.CaseIndices;
assert(isnumeric(indices)&&isvector(indices)&&~isempty(indices)&& ...
    all(ismember(indices,1:3))&&numel(unique(indices))==numel(indices), ...
    'CaseIndices must contain distinct members of 1:3.');
schedule=schedule(indices,:);
if strlength(outputDirectory)==0
    outputDirectory=string(config.OutputDirectory)+"_tolerance_study_qprecovery_v250_100_600_3600";
end
assert(~strcmp(char(outputDirectory),char(config.OutputDirectory)), ...
    'A separate output directory is required.');
for folder=[outputDirectory,fullfile(outputDirectory,'runs'), ...
        fullfile(outputDirectory,'problems'),fullfile(outputDirectory,'tables'), ...
        fullfile(outputDirectory,'sessions')]
    if ~isfolder(folder), mkdir(folder); end
end
manifest=struct('schemaVersion',"2.5.0",'config',config, ...
    'schedule',schedule,'problemIds',problemIds,'settings',settings, ...
    'matlabSources',localSources(matlabRoot), ...
    'pythonSources',localSources(fullfile(packageRoot,'python')), ...
    'extensionSources',localSources(extensionRoot));
manifestFile=fullfile(outputDirectory,'tolerance_study_config.mat');
if isfile(manifestFile)
    prior=load(manifestFile,'manifest');
    assert(isequaln(prior.manifest,manifest), ...
        'Configuration or source changed: use a new output directory.');
else
    localSave(manifestFile,struct('manifest',manifest));
end
writetable(schedule,fullfile(outputDirectory,'tolerance_schedule.csv'));

% Actual Gurobi responses must pass the captured regressions before a long run.
if runBnB && settings.RunSelfTests
    qpPreflight=check_nist_qp_recovery(config.QPOptions);
    preflightFile=fullfile(outputDirectory,'qp_recovery_preflight.mat');
    localSave(preflightFile,struct('qpPreflight',qpPreflight));
    assert(qpPreflight.passed,'QP preflight failed. Share this file: %s',preflightFile);
end

% All solver/module initialization occurs outside the timed benchmark calls.
if exist('warmup_nist_solvers','file')~=0
    warmup_nist_solvers(config.QPOptions);
end
bridgeInfo=struct(); pythonVersions=struct();
if runSCIP
    [module,bridgeInfo]=import_nist_scip_bridge(PackageRoot=packageRoot,Reload=true);
    smoke=jsondecode(string(module.smoke_test_json()));
    assert(smoke.has_solution,'SCIP smoke test did not return a solution.');
    pythonVersions=jsondecode(string(module.versions_json()));
end
environment=struct('matlab',version,'architecture',computer('arch'), ...
    'optimizationToolbox',ver('optim'),'gurobiInterface',which('gurobi'), ...
    'bridgeInfo',bridgeInfo,'pythonVersions',pythonVersions, ...
    'solvers',settings.Solvers,'threadCount',1,'startedAt',string(datetime('now')));
sessionTag=string(datetime('now','Format','yyyyMMdd_HHmmss'));
localSave(fullfile(outputDirectory,'sessions','environment_'+sessionTag+'.mat'), ...
    struct('environment',environment));

% Prepare each common initial point once; retain the same problem across cases.
problems=cell(numel(problemIds),1);
for i=1:numel(problemIds)
    problemFile=fullfile(outputDirectory,'problems',lower(problemIds(i))+'_study.mat');
    if isfile(problemFile)
        p=load(problemFile,'problem'); problem=p.problem;
    else
        originalProblem=fullfile(config.OutputDirectory,'problems',lower(problemIds(i))+'_study.mat');
        if isfile(originalProblem)
            p=load(originalProblem,'problem'); problem=validate_nist_problem(p.problem);
        else
            problem=load_nist_problem(problemIds(i), ...
                DataDirectory=config.DataDirectory,AutoDownload=false,QPOptions=config.QPOptions);
        end
        localSave(problemFile,struct('problem',problem));
    end
    problems{i}=problem;
end

rows=cell(0,1);
for c=1:height(schedule)
    caseId=schedule.case_id(c); relative=schedule.epsilon_relative(c);
    timeLimit=schedule.time_limit(c); absolute=1e-10;
    caseFolder=fullfile(outputDirectory,'runs',caseId);
    if ~isfolder(caseFolder), mkdir(caseFolder); end
    for i=1:numel(problemIds)
        problem=problems{i}; id=problemIds(i);
        order=settings.Solvers;
        if numel(order)>1&&mod(i+c,2)==0, order=fliplr(order); end
        for solver=order
            runFile=fullfile(caseFolder,lower(id)+'_'+lower(solver)+'.mat');
            failed=false;
            if isfile(runFile)
                cached=load(runFile,'row'); row=cached.row;
                fprintf('Resume: %s / %s / %s\n',caseId,id,solver);
            else
                fprintf('\nRunning %s / %s / %s: eps_r=%g, limit=%g s\n', ...
                    caseId,id,solver,relative,timeLimit);
                errorMessage=""; options=struct();
                try
                    switch solver
                        case "BnB"
                            options=config.BnBOptions;
                            options.TimeLimit=timeLimit;
                            options.EpsilonRelative=relative;
                            options.EpsilonAbsolute=absolute;
                            options.RecordRelativeTargets=[0.1,0.01,0.001];
                            options.QPOptions=config.QPOptions;
                            options.QPOptions.Threads=1;
                            result=solve_nist_bnb_tolerance_study(problem,options);
                        case "SCIP"
                            options=config.SCIPOptions;
                            options.time_limit=timeLimit;
                            options.paper_epsilon_relative=relative;
                            options.paper_epsilon_absolute=absolute;
                            options.threads=1;
                            options.reload_python_module=false;
                            options.output_directory='';
                            if options.save_artifacts
                                options.output_directory=char(fullfile(caseFolder,lower(id)+'_scip_artifacts'));
                            end
                            result=solve_nist_scip(problem,options);
                    end
                catch exception
                    failed=true;
                    errorMessage=string(getReport(exception,'extended','hyperlinks','off'));
                    result=struct('solver',solver,'status',"error", ...
                        'problem_id',id,'instance_id',string(problem.instance_id));
                end
                row=nist_tolerance_result_row(result,problem,caseId,relative,absolute,timeLimit,errorMessage);
                localSave(runFile,struct('result',result,'row',row, ...
                    'options',options,'environment',environment));
            end
            rows{end+1,1}=row; %#ok<AGROW>
            results=vertcat(rows{:});
            summary=write_nist_tolerance_outputs(results,outputDirectory);
            localSave(fullfile(outputDirectory,'nist_tolerance_results.mat'), ...
                struct('results',results,'summary',summary,'manifest',manifest));
            fprintf('%s / %s / %s: U=%.9g, L=%.9g (RSS/s), t=%.3f s, nodes=%.0f, tol_only=%.0f, gap_ok=%d\n', ...
                caseId,id,solver,row.objective_scaled,row.dual_bound_scaled, ...
                row.wall_time,row.nodes_visited,row.pruned_tolerance_only,row.certified);
            if row.solver_certified && ~row.certified
                warning('NISTStudy:GapMismatch', ...
                    'The solver flag for %s/%s differs from the common exported gap test.',id,solver);
            end
            if failed && settings.StopOnError
                error('NISTStudy:SolverFailure', ...
                    'Run saved with an error. See error_message in the CSV/MAT; fix and remove that run checkpoint before retrying.');
            end
        end
    end
end
if settings.RunVarPro
    varproRows=cell(numel(problemIds),1);
    for i=1:numel(problemIds)
        file=fullfile(outputDirectory,'runs',lower(problemIds(i))+'_varpro.mat');
        if isfile(file)
            cached=load(file,'row'); row=cached.row;
        else
            result=solve_nist_varpro_multistart(problems{i},config.VarProOptions);
            row=nist_result_to_row(result,config.VarProOptions.TimeLimit);
            row.objective_scale=problems{i}.objectiveScale;
            localSave(file,struct('result',result,'row',row,'environment',environment));
        end
        varproRows{i}=row;
    end
    writetable(vertcat(varproRows{:}),fullfile(outputDirectory,'varpro_results_scaled.csv'));
end
fprintf('\nMain results: %s\n',fullfile(outputDirectory,'nist_tolerance_results.csv'));
fprintf('LaTeX tables: %s\n',fullfile(outputDirectory,'tables'));
disp(summary)
end

function text=localNormalize(text)
text=strtrim(strrep(text,sprintf('\r\n'),sprintf('\n')));
end

function sources=localSources(folder)
% Relative filenames and contents only: timestamps do not invalidate resume.
files=[dir(fullfile(folder,'**','*.m'));dir(fullfile(folder,'**','*.py')); ...
    dir(fullfile(folder,'**','*.txt'));dir(fullfile(folder,'**','*.json'))];
names=strings(numel(files),1); contents=cell(numel(files),1);
for i=1:numel(files)
    file=fullfile(files(i).folder,files(i).name);
    names(i)=string(file(numel(char(folder))+2:end));
    contents{i}=localNormalize(fileread(file));
end
[names,index]=sort(names); contents=contents(index);
sources=struct('names',names,'contents',{contents});
end

function localSave(file,payload)
temporary=string(file)+'.tmp.mat';
save(temporary,'-struct','payload','-v7.3');
movefile(temporary,file,'f');
end
