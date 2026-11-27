function options = nist_scip_options(varargin)
%NIST_SCIP_OPTIONS Default options for the SCIP comparison.
p = inputParser;
addParameter(p,"TimeLimit",600);
addParameter(p,"EpsilonRelative",0.01);
addParameter(p,"EpsilonAbsolute",1e-10);
addParameter(p,"Threads",1);
addParameter(p,"Verbosity",2);
addParameter(p,"FeasibilityTolerance",1e-6);
addParameter(p,"DisablePresolve",true);
addParameter(p,"SaveArtifacts",false);
addParameter(p,"OutputDirectory","");
addParameter(p,"ReloadPythonModule",false);
parse(p,varargin{:});
options = struct();
options.time_limit = double(p.Results.TimeLimit);
options.paper_epsilon_relative = double(p.Results.EpsilonRelative);
options.paper_epsilon_absolute = double(p.Results.EpsilonAbsolute);
options.threads = double(p.Results.Threads);
options.verbosity = double(p.Results.Verbosity);
options.feasibility_tolerance = double(p.Results.FeasibilityTolerance);
options.disable_presolve = logical(p.Results.DisablePresolve);
options.save_artifacts = logical(p.Results.SaveArtifacts);
options.output_directory = char(string(p.Results.OutputDirectory));
options.reload_python_module = logical(p.Results.ReloadPythonModule);
options.random_seed_shift = 0;
options.permutation_seed = 0;
options.lp_seed = 0;
end
