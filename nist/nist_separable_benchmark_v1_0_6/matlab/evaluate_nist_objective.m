function [objective,residual,prediction] = evaluate_nist_objective(problem,linear,z)
%EVALUATE_NIST_OBJECTIVE Evaluate residuals and residual sum of squares.
problem = validate_nist_problem(problem);
linear = double(linear(:));
z = double(z(:));
if numel(linear) ~= problem.nLinear
    error("The linear parameter vector has the wrong length.");
end
Phi = nist_basis_matrix(problem,z);
prediction = Phi*linear;
residual = prediction-problem.response;
objective = residual.'*residual;
end
