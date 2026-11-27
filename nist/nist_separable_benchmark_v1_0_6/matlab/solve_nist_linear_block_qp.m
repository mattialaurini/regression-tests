function result = solve_nist_linear_block_qp(problem,z,qpOptions)
%SOLVE_NIST_LINEAR_BLOCK_QP Bounded linear least squares at fixed z.
arguments
    problem (1,1) struct
    z double
    qpOptions (1,1) struct = nist_qp_options()
end
problem = validate_nist_problem(problem);
Phi = nist_basis_matrix(problem,z);
y = problem.response;
d = problem.nLinear;
qp.Q = sparse(Phi.'*Phi);
qp.f = -2*Phi.'*y;
qp.objcon = y.'*y;
qp.Aineq = sparse(0,d); qp.bineq = zeros(0,1);
qp.Aeq = sparse(0,d); qp.beq = zeros(0,1);
qp.lb = problem.linearLower; qp.ub = problem.linearUpper;
raw = solve_nist_qp(qp,qpOptions);
result = raw;
result.linear = [];
result.residual = [];
result.prediction = [];
if raw.hasSolution
    linear = min(max(raw.x,problem.linearLower),problem.linearUpper);
    [objective,residual,prediction] = evaluate_nist_objective(problem,linear,z);
    result.linear = linear;
    result.residual = residual;
    result.prediction = prediction;
    result.objective = objective;
end
end
