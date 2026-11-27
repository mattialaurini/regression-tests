function export_pkpd_results(arxFile,mlFile,outputDirectory)
%EXPORT_PKPD_RESULTS Generate the manuscript ARX/ML tables, macros and figures.
arguments
    arxFile (1,1) string = ""
    mlFile (1,1) string = ""
    outputDirectory (1,1) string = ""
end
root=string(fileparts(mfilename('fullpath')));
if strlength(arxFile)==0, arxFile=fullfile(root,'outputs','pkpd_arx','results','tests_lin.mat'); end
if strlength(mlFile)==0, mlFile=fullfile(root,'outputs','pkpd_ml','results','tests_ml.mat'); end
if strlength(outputDirectory)==0, outputDirectory=fullfile(root,'outputs','paper'); end
oldPath=path; cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'pkpd','src'),'-begin');
assert(isfile(arxFile)||isfile(mlFile),'No completed cohort result files were found.');
if isfile(arxFile), paper_lin_outputs(char(arxFile),char(outputDirectory)); end
if isfile(mlFile), paper_ml_outputs(char(mlFile),char(outputDirectory)); end
end
