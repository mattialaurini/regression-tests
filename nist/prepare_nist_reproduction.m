function state=prepare_nist_reproduction()
%PREPARE_NIST_REPRODUCTION Materialize portable frozen inputs, never solve a QP.
root=string(fileparts(mfilename('fullpath')));
state=struct('root',root,'baseRoot',fullfile(root,'nist_separable_benchmark_v1_0_6'), ...
    'toleranceRoot',fullfile(root,'nist_tolerance_study'), ...
    'ablationRoot',fullfile(root,'nist_diagonal_grid_v250'));
previousPath=path; cleanup=onCleanup(@()path(previousPath)); %#ok<NASGU>
addpath(fullfile(state.baseRoot,'matlab'),'-begin');
d=load(fullfile(root,'data','frozen_inputs.mat'));
catalog=nist_problem_catalog();
assert(iscell(d.problems)&&numel(d.problems)==8&&iscell(d.baseOptions));
problems=cell(1,8); bases=d.baseOptions;
for k=1:8
    frozen=d.problems{k}; j=find(strcmpi([catalog.id],string(frozen.id)));
    assert(isscalar(j),'Unknown frozen problem.');
    problem=catalog(j);
    for name=reshape(string(fieldnames(frozen)),1,[])
        if name~="id", problem.(name)=frozen.(name); end
    end
    problem.domainProfile="study";
    problem.instance_id=lower(string(problem.id))+"_study";
    problems{k}=validate_nist_problem(problem);
    bases{k}=restoreOptionStrings(bases{k});
end
state.problems=problems; state.baseOptions=bases;
config=d.config;
config.Profile="full"; config.ProblemIds=string(cellfun(@(p)char(p.id),problems,'UniformOutput',false));
config.BnBOptions=restoreOptionStrings(config.BnBOptions);
config.QPOptions.Solver=string(config.QPOptions.Solver);
config.VarProOptions.Display="off"; config.VarProOptions.QPOptions.Solver="gurobi";
config.SCIPOptions.output_directory='';
config.DataDirectory=fullfile(state.baseRoot,'data','raw');
config.OutputDirectory=fullfile(state.baseRoot,'output','reproduction_seed');
state.seedMat=fullfile(config.OutputDirectory,'benchmark_results.mat');
state.referenceMat=fullfile(root,'.reproduction','ablation_reference.mat');
% Persist exactly the same data/initial point for both independent solvers.
for k=1:8
    problem=problems{k};
    file=fullfile(config.OutputDirectory,'problems',lower(string(problem.id))+'_study.mat');
    if ~isfolder(fileparts(file)), mkdir(fileparts(file)); end
    save(file,'problem','-v7');
end
save(state.seedMat,'config','-v7');
manifest=struct('version',"portable-reference-v1",'problems',{problems}, ...
    'baseOptions',{bases},'matlabSources',sourceSnapshot(fullfile(state.baseRoot,'matlab')), ...
    'experimentSources',sourceSnapshot(state.ablationRoot));
state.manifest=manifest;
if ~isfolder(fileparts(state.referenceMat)), mkdir(fileparts(state.referenceMat)); end
save(state.referenceMat,'manifest','-v7');
end
function o=restoreOptionStrings(o)
o.BranchingRule=string(o.BranchingRule); o.QPOptions.Solver=string(o.QPOptions.Solver);
end
function s=sourceSnapshot(folder)
files=[dir(fullfile(folder,'**','*.m'));dir(fullfile(folder,'**','*.py')); ...
    dir(fullfile(folder,'**','*.txt'));dir(fullfile(folder,'**','*.json'))];
names=strings(numel(files),1); contents=cell(numel(files),1);
for j=1:numel(files)
    f=fullfile(files(j).folder,files(j).name); names(j)=string(f(numel(char(folder))+2:end));
    contents{j}=strtrim(strrep(fileread(f),sprintf('\r\n'),sprintf('\n')));
end
[names,ix]=sort(names); s=struct('names',names,'contents',{contents(ix)});
end
