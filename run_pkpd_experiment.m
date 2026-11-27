function result=run_pkpd_experiment(experiment,args)
%RUN_PKPD_EXPERIMENT Isolated, non-interactive drivers for manuscript experiments.
arguments
    experiment (1,1) string {mustBeMember(experiment,["arx","ml","patient13"])}
    args.Sigmas (1,:) double = []
    args.PatientIds (1,:) double = 1:13
    args.OutputDirectory (1,1) string = ""
    args.TimeLimit (1,1) double = NaN
    args.Resume (1,1) logical = true
    args.DryRun (1,1) logical = false
    args.Export (1,1) logical = true
end
root=string(fileparts(mfilename('fullpath'))); src=fullfile(root,'pkpd','src');
output=args.OutputDirectory;
if strlength(output)==0, output=fullfile(root,'outputs','pkpd_'+experiment); end
if ~isAbsolute(output), output=fullfile(pwd,output); end
sigmas=args.Sigmas;
if isempty(sigmas)
    if experiment=="patient13", sigmas=1; else, sigmas=[0,0.01,0.1,1]; end
end
validateattributes(sigmas,{'double'},{'real','finite','nonnegative','vector'});
validateattributes(args.PatientIds,{'double'},{'integer','>=',1,'<=',13,'nonempty','vector'});
assert(numel(unique(sigmas))==numel(sigmas),'Duplicate noise level.');
assert(numel(unique(args.PatientIds))==numel(args.PatientIds),'Duplicate patient id.');
timeLimit=args.TimeLimit;
if isnan(timeLimit)
    if experiment=="arx", timeLimit=1800; else, timeLimit=3600; end
end
validateattributes(timeLimit,{'double'},{'scalar','finite','positive'});
if experiment=="patient13", patients=13; models=6; else, patients=args.PatientIds;
    if experiment=="arx", models=8; else, models=4; end
end
plan=struct('experiment',experiment,'patients',patients,'sigmas',sigmas, ...
    'main_trials',numel(patients)*numel(sigmas)*models,'time_limit',timeLimit, ...
    'output_directory',output);
if experiment=="arx", plan.tuning="adaptive grids; patients [1 6 9], sigma=0.1, nominal cap 600 s"; end
if experiment=="patient13"
    plan.time_limit="ARX=1800 s, ML="+string(timeLimit)+" s";
end
if args.DryRun, disp(plan); result=plan; return, end
oldPath=path; oldFolder=pwd;
addpath(src,'-begin');
oldSolver=qp_solver(); oldThreads=gurobi_threads(); oldPresolve=gurobi_presolve();
cleanup=onCleanup(@()restoreSession(oldFolder,oldPath,oldSolver,oldThreads,oldPresolve)); %#ok<NASGU>
for name=["ss","c2d","gurobi"]
    assert(~isempty(which(char(name))),'Required function is unavailable: %s',name);
end
if ~isfolder(output), mkdir(output); end
cfg=struct('ids',args.PatientIds,'sigmas',sigmas,'omega',0.2, ...
    'timeout',timeLimit,'Reltol',0.05,'Abstol',1e-3,'parallel',false, ...
    'gurobi_threads',0,'gurobi_presolve',-1,'qp_solver','hybrid','progress',15);
manifest=struct('protocol',"pkpd-portable-1.0.0",'experiment',experiment, ...
    'cfg',cfg,'sources',sourceSnapshot(src));
manifestFile=fullfile(output,'experiment_config.mat');
if isfile(manifestFile)
    previous=load(manifestFile,'manifest');
    assert(isequaln(previous.manifest,manifest), ...
        'Settings or source changed. Use a new OutputDirectory.');
else
    save(manifestFile,'manifest','-v7');
end
record_pkpd_environment(output);
cd(output); % Legacy batch scripts write relative folders; keep them out of src.
if experiment=="patient13"
    file=fullfile(output,'patient13_checkpoint.mat');
    if ~args.Resume, archiveExisting(file); end
    result=paper_patient13_six_models('sigmas',sigmas,'outdir',char(output), ...
        'solver','hybrid','timeout_arx',1800,'timeout_ml',timeLimit);
    return
end
folder=fullfile(output,'results'); if ~isfolder(folder), mkdir(folder); end
file=fullfile(folder,'tests_'+string(experiment)+'.mat');
% Preserve the original ARX filename consumed by the manuscript output script.
if experiment=="arx", file=fullfile(folder,'tests_lin.mat'); end
cfg.out=char(file);
if experiment=="ml", cfg.estimate_e0=true; end
if isfile(file)
    if args.Resume
        previous=load(file);
        if experiment=="arx"
            assert(isfield(previous,'best'), ...
                ['The saved tuning phase is incomplete. Use Resume=false to archive it ' ...
                 'and restart tuning, or choose a new OutputDirectory.']);
        end
        cfg.resume=char(file);
    else
        archiveExisting(file);
    end
end
if experiment=="arx"
    runARX(cfg);
    if args.Export, paper_lin_outputs(char(file),char(fullfile(output,'paper'))); end
else
    runML(cfg);
    if args.Export, paper_ml_outputs(char(file),char(fullfile(output,'paper'))); end
end
result=load(file);
end
function runARX(cfg)
% A function workspace prevents stale cfg/settings variables from other batches.
paper_tests_lin
end
function runML(cfg)
paper_tests_ml
end
function archiveExisting(file)
if ~isfile(file), return, end
[folder,name,ext]=fileparts(file);
backup=fullfile(folder,string(name)+'_backup_'+string(datetime('now','Format','yyyyMMdd_HHmmssSSS'))+string(ext));
movefile(file,backup); fprintf('Previous checkpoint archived: %s\n',backup);
end
function restoreSession(folder,previousPath,solver,threads,presolve)
qp_solver(solver); gurobi_threads(threads); gurobi_presolve(presolve);
cd(folder); path(previousPath);
end
function s=sourceSnapshot(folder)
files=dir(fullfile(folder,'*.m')); names=string({files.name})'; contents=cell(numel(files),1);
for j=1:numel(files), contents{j}=fileread(fullfile(folder,files(j).name)); end
s=struct('names',names,'contents',{contents});
end
function yes=isAbsolute(p)
c=char(p); yes=startsWith(p,filesep)||(~isempty(regexp(c,'^[A-Za-z]:[\\/]','once')));
end
