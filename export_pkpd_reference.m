function export_pkpd_reference(outputDirectory)
%EXPORT_PKPD_REFERENCE Rebuild patient-13 figures from stored fits; no fitting.
arguments
    outputDirectory (1,1) string = ""
end
root=string(fileparts(mfilename('fullpath')));
if strlength(outputDirectory)==0, outputDirectory=fullfile(root,'outputs','pkpd_reference'); end
if ~isfolder(outputDirectory), mkdir(outputDirectory); end
oldPath=path; cleanup=onCleanup(@()path(oldPath)); %#ok<NASGU>
addpath(fullfile(root,'pkpd','src'),'-begin');
stored=load(fullfile(root,'pkpd','reference','patient13_fits.mat'),'R');
export_patient13_outputs(stored.R,char(outputDirectory));
end
