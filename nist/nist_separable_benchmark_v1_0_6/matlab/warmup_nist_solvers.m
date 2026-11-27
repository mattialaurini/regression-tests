function info = warmup_nist_solvers(qpOptions)
%WARMUP_NIST_SOLVERS Load Gurobi and SCIP before timed benchmark runs.
arguments
    qpOptions (1,1) struct = nist_qp_options()
end
qp=struct(); qp.Q=speye(1); qp.f=-2; qp.objcon=1;
qp.Aineq=sparse(0,1); qp.bineq=zeros(0,1); qp.Aeq=sparse(0,1); qp.beq=zeros(0,1);
qp.lb=0; qp.ub=2;
qpResult=solve_nist_qp(qp,qpOptions);
assert(qpResult.hasSolution,"QP warm-up failed.");
[module,bridgeInfo]=import_nist_scip_bridge();
scipResult=jsondecode(string(module.smoke_test_json()));
assert(logical(scipResult.has_solution),"SCIP warm-up failed.");
info=struct("qp_solver",string(qpResult.solver),"qp_objective",qpResult.objective, ...
    "scip_status",string(scipResult.status),"scip_objective",double(scipResult.objective), ...
    "bridge_file",string(bridgeInfo.module_file));
end
