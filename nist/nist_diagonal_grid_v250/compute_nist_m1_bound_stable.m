function result = compute_nist_m1_bound_stable(problem,region,nPoints,qpOptions)
%COMPUTE_NIST_M1_BOUND Valid interpolation lower bound for RSS.
startTime = tic;
data = nist_build_relaxation_data(problem,region,nPoints);
result = solve_nist_outer_relaxation_qp_stable(data.G,data.cMinus,data.cPlus,qpOptions);
result.value = result.safeLowerBound;
result.rawObjective = result.objective;
result.totalWallTime = toc(startTime);
result.nColumns = data.nColumns;
result.nPoints = nPoints;
end
