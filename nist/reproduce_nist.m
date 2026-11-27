function [results,summary]=reproduce_nist(experiment,args)
%REPRODUCE_NIST Fresh or resumed NIST comparison/ablation, without GUI dialogs.
arguments
    experiment (1,1) string {mustBeMember(experiment,["bnb","scip","ablation"])}
    args.CaseIndices (1,:) double = 1:3
    args.ProblemIds (1,:) string = strings(1,0)
    args.OutputDirectory (1,1) string = ""
    args.DryRun (1,1) logical = false
    args.RunSelfTests (1,1) logical = true
    args.BudgetMode (1,1) string {mustBeMember(args.BudgetMode,["published","uniform"])} = "published"
end
setup_reproduction();
state=prepare_nist_reproduction();
previousPath=path; cleanup=onCleanup(@()path(previousPath)); %#ok<NASGU>
output=args.OutputDirectory;
if strlength(output)==0, output=fullfile(state.root,'outputs',experiment); end
assert(all(ismember(args.CaseIndices,1:3))&&~isempty(args.CaseIndices)&& ...
    numel(unique(args.CaseIndices))==numel(args.CaseIndices),'Invalid CaseIndices.');
allIds=string(cellfun(@(p)char(p.id),state.problems,'UniformOutput',false));
ids=args.ProblemIds; if isempty(ids), ids=allIds; end
assert(all(ismember(ids,allIds))&&numel(unique(ids))==numel(ids),'Invalid ProblemIds.');
if experiment=="ablation"
    addpath(state.ablationRoot,'-begin');
    if ~args.DryRun
        addpath(fullfile(state.baseRoot,'matlab'));
        record_reproduction_environment(fullfile(output,'sessions'),false);
    end
    [results,summary]=run_ablation_reproduction(state,output,args,ids);
else
    addpath(state.toleranceRoot,'-begin');
    if args.DryRun
        er=[0.1,0.01,0.001]; budgets=[100,600,3600]; rows=cell(0,1);
        for c=args.CaseIndices
            for id=ids
                rows{end+1,1}=table(id,upper(experiment),er(c),budgets(c),5,2,3, ...
                    'VariableNames',{'problem_id','solver','epsilon_relative','time_limit','M','N','Nhat'}); %#ok<AGROW>
            end
        end
        results=vertcat(rows{:}); summary=table(); disp(results);
        fprintf('Dry run: no optimizer called. Comparison grids apply to BnB only.\n');
        return
    end
    solver="BnB"; if experiment=="scip", solver="SCIP"; end
    addpath(fullfile(state.baseRoot,'matlab'));
    record_reproduction_environment(fullfile(output,'sessions'),experiment=="scip");
    settings=struct('CaseIndices',args.CaseIndices,'ProblemIds',ids, ...
        'Solvers',solver,'RunVarPro',false,'RunSelfTests',args.RunSelfTests,'StopOnError',true);
    [results,summary]=run_nist_tolerance_study(state.seedMat,output,settings);
end
end
