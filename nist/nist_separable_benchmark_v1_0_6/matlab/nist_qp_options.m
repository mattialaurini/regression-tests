function options = nist_qp_options(varargin)
%NIST_QP_OPTIONS Options for convex QP subproblems.
p = inputParser;
addParameter(p,"Solver","gurobi");
addParameter(p,"Threads",1);
addParameter(p,"Display",false);
addParameter(p,"FeasibilityTolerance",1e-9);
addParameter(p,"OptimalityTolerance",1e-9);
addParameter(p,"TimeLimit",Inf);
addParameter(p,"BoundSafety",1e-8);
parse(p,varargin{:});
options = p.Results;
options.Solver = string(options.Solver);
end
