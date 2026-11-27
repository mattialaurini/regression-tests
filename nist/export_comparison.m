function summary=export_comparison()
%EXPORT_COMPARISON Merge the separately executed BnB and SCIP comparison runs.
root=string(fileparts(mfilename('fullpath')));
a=load(fullfile(root,'outputs','bnb','nist_tolerance_results.mat'),'results');
b=load(fullfile(root,'outputs','scip','nist_tolerance_results.mat'),'results');
results=[a.results;b.results];
assert(height(unique(results(:,{'case_id','problem_id','solver'})))==height(results), ...
    'Duplicate comparison records.');
folder=fullfile(root,'outputs','comparison'); if ~isfolder(folder), mkdir(folder); end
oldPath=path; cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'nist_tolerance_study'));
summary=write_nist_tolerance_outputs(results,folder);
save(fullfile(folder,'nist_tolerance_results.mat'),'results','summary','-v7.3');
end
