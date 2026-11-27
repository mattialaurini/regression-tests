function result = compute_nist_m2_screening_stable(problem,region,nPoints,M,qpOptions)
%COMPUTE_NIST_M2_SCREENING Coarse-grid screening quantity.
if M<nPoints, error("M must be at least nPoints."); end
startTime = tic;
data = nist_build_relaxation_data(problem,region,nPoints);
scale = ((nPoints-1)/(M-1))^2;
result = solve_nist_outer_relaxation_qp_stable( ...
    data.G,scale*data.cMinus,scale*data.cPlus,qpOptions,"screening");
result.screeningAvailable=isfinite(result.feasiblePrimalUpperBound);
if result.screeningAvailable
    result.value=result.feasiblePrimalUpperBound;
else
    % Unavailable screening MUST force the m1 path, never bypass it.
    result.value=Inf;
end
result.totalWallTime = toc(startTime);
result.nColumns = data.nColumns;
result.nPoints = nPoints;
result.M = M;
result.scale = scale;
end
