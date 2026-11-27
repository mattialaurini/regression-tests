function s = qp_solver(new)
%QP_SOLVER Get or set the solver used by SOLVE_MODEL.
%   qp_solver('gurobi')  Gurobi (default when available)
%   qp_solver('matlab')  quadprog / linprog (Optimization Toolbox)
%   qp_solver('native')  structured first-order solver for the 2-norm bounds
%                        (NATIVE_L2: no external solver; LPs/1-norm still use linprog)
%   qp_solver('hybrid')  native solver for upper bounds and as a first pass for
%                        lower bounds, Gurobi only when the native bracket does not
%                        decide pruning against the current threshold (HYBRID_TAU)
%   qp_solver('auto')    Gurobi if found on the path, otherwise quadprog/linprog
%   s = qp_solver        current choice ('gurobi' or 'matlab')
persistent choice
if nargin > 0, choice = lower(new); end
if isempty(choice), choice = 'auto'; end
s = choice;
if strcmp(s, 'auto')
    if exist('gurobi', 'file') == 2 || exist('gurobi', 'file') == 3, s = 'gurobi'; else, s = 'matlab'; end
end
end
