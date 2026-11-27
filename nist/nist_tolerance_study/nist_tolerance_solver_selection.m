function solvers = nist_tolerance_solver_selection(value)
%NIST_TOLERANCE_SOLVER_SELECTION Canonical names, preserving requested order.
solvers=reshape(string(value),1,[]);
canonical=["BnB","SCIP"];
[valid,index]=ismember(lower(solvers),lower(canonical));
assert(~isempty(solvers)&&all(valid)&&numel(unique(index))==numel(index), ...
    'NISTStudy:SolverSelection','Select BnB, SCIP, or both without duplicates.');
solvers=canonical(index);
end
