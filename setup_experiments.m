function setup_experiments(pythonExecutable)
%SETUP_EXPERIMENTS Add the unified entry points without starting any solver.
arguments
    pythonExecutable (1,1) string = ""
end
root=string(fileparts(mfilename('fullpath')));
addpath(root); addpath(fullfile(root,'nist'));
if strlength(pythonExecutable)>0, setup_reproduction(pythonExecutable); end
fprintf('Regression experiments ready: %s\n',root);
end
