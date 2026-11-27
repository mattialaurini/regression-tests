function setup_reproduction(pythonExecutable)
%SETUP_REPRODUCTION Add entry points; optionally select a Python interpreter.
arguments
    pythonExecutable (1,1) string = ""
end
root=string(fileparts(mfilename('fullpath')));
addpath(root);
% Set these before Python/BLAS initialization, not during timed calls.
for name=["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS", ...
        "VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS"]
    setenv(name,'1');
end
if strlength(pythonExecutable)>0
    assert(isfile(pythonExecutable),'Python executable not found: %s',pythonExecutable);
    pe=pyenv;
    if pe.Status=="Loaded"
        assert(string(pe.Executable)==pythonExecutable, ...
            'Python is already loaded from another interpreter. Restart MATLAB first.');
    else
        pyenv('Version',pythonExecutable,'ExecutionMode','OutOfProcess');
    end
end
fprintf('NIST reproduction entry points ready: %s\n',root);
end
