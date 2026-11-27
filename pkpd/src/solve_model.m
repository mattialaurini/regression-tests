function result = solve_model(model, params)
%SOLVE_MODEL Solve a Gurobi-format LP/QP with the solver chosen by QP_SOLVER.
%   Fields used: A, rhs, sense ('<', '>', '='), obj, lb, ub, optional Q (objective
%   x'Qx + obj'x), modelsense 'min'. The result has the Gurobi fields status
%   ('OPTIMAL' or other), objval, x and, when available, pi (constraint duals in
%   Gurobi's sign convention; the certificates in NORM_MIN_LP try both signs).
%   With 'matlab': quadprog (interior-point-convex) for QPs, linprog (dual simplex)
%   for LPs. Quadratic constraints (quadcon) require Gurobi.
if any(strcmp(qp_solver(), {'gurobi', 'hybrid'})) || (strcmp(qp_solver(), 'native') && exist('gurobi', 'file') >= 2 && ~exist('linprog', 'file'))
    result = gurobi(model, params);
    return
end
if isfield(model, 'quadcon')
    error('solve_model:quadcon', 'Quadratic constraints require Gurobi (qp_solver(''gurobi'')).');
end
A = sparse(model.A); rhs = model.rhs(:); sense = model.sense(:);
eq = sense == '='; le = sense == '<'; ge = sense == '>';
Ain = [A(le, :); -A(ge, :)]; bin = [rhs(le); -rhs(ge)];
Aeq = A(eq, :); beq = rhs(eq);
lb = model.lb(:); ub = model.ub(:); f = model.obj(:);
isqp = isfield(model, 'Q') && nnz(model.Q) > 0;
result = struct('status', 'FAILED');
try
    if isqp
        H = 2 * sparse(model.Q); H = (H + H') / 2;
        opt = qp_options('quadprog');
        [x, fval, flag, ~, lam] = quadprog(H, f, Ain, bin, Aeq, beq, lb, ub, [], opt);
    else
        opt = qp_options('linprog');
        [x, fval, flag, ~, lam] = linprog(f, Ain, bin, Aeq, beq, lb, ub, opt);
    end
catch err
    warning('solve_model:failed', '%s', err.message); return
end
if isempty(x) || any(~isfinite(x)), return; end
result.x = x(:); result.objval = fval;
if flag > 0, result.status = 'OPTIMAL'; else, result.status = sprintf('EXITFLAG_%d', flag); end
if isstruct(lam) && isfield(lam, 'ineqlin') && isfield(lam, 'eqlin')
    pi = zeros(numel(rhs), 1); nl = sum(le);
    li = lam.ineqlin(:); le_ = find(le); ge_ = find(ge);
    if numel(li) == nl + numel(ge_) && numel(lam.eqlin) == sum(eq)
        pi(le_) = -li(1:nl); pi(ge_) = li(nl+1:end); pi(eq) = -lam.eqlin(:);
        result.pi = pi;
    end
end
end

function opt = qp_options(which)
try
    if strcmp(which, 'quadprog')
        opt = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex', ...
            'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-12);
    else
        opt = optimoptions('linprog', 'Display', 'off', 'Algorithm', 'dual-simplex');
    end
catch
    opt = optimset('Display', 'off');                % older MATLAB / Octave
end
end
