function result = solve_nist_outer_relaxation_qp_stable(G,cMinus,cPlus,qpOptions,purpose)
%SOLVE_NIST_OUTER_RELAXATION_QP Minimize squared two-norm over outer image.
arguments
    G double
    cMinus double
    cPlus double
    qpOptions (1,1) struct = nist_qp_options()
    purpose (1,1) string = "lower"
end
[m,r] = size(G);
if ~isequal(size(cMinus),[m,r]) || ~isequal(size(cPlus),[m,r])
    error("Interpolation matrices must match G.");
end
raw=nist_solve_outer_qp_verified(G,cMinus,cPlus,qpOptions,purpose);
result=raw;
if ~isfield(result,'feasiblePrimalUpperBound')
    [result.feasiblePrimalUpperBound,result.weights]= ...
        nist_outer_qp_primal_upper(G+cMinus,G+cPlus,raw,qpOptions);
end
end
