function problem = validate_nist_problem(problem)
%VALIDATE_NIST_PROBLEM Validate dimensions and domain assumptions.
required = ["id","modelId","response","predictors","linearLower", ...
    "linearUpper","zLower","zUpper","linearVertices","nLinear", ...
    "nNonlinear","nResidual","residualScale","objectiveScale"];
for name = required
    if ~isfield(problem,name)
        error("NIST problem is missing field '%s'.",name);
    end
end
problem.response = double(problem.response(:));
problem.predictors = double(problem.predictors);
problem.linearLower = double(problem.linearLower(:));
problem.linearUpper = double(problem.linearUpper(:));
problem.zLower = double(problem.zLower(:));
problem.zUpper = double(problem.zUpper(:));
problem.linearVertices = double(problem.linearVertices);
if size(problem.predictors,1) ~= problem.nResidual || ...
        numel(problem.response) ~= problem.nResidual
    error("Data dimensions are inconsistent.");
end
if numel(problem.linearLower) ~= problem.nLinear || ...
        numel(problem.linearUpper) ~= problem.nLinear
    error("Linear-bound dimensions are inconsistent.");
end
if numel(problem.zLower) ~= problem.nNonlinear || ...
        numel(problem.zUpper) ~= problem.nNonlinear
    error("Nonlinear-bound dimensions are inconsistent.");
end
if any(problem.linearLower>=problem.linearUpper) || any(problem.zLower>=problem.zUpper)
    error("Parameter bounds are invalid.");
end
if size(problem.linearVertices,1) ~= problem.nLinear
    error("linearVertices has an invalid row dimension.");
end
if any(~isfinite(problem.response)) || any(~isfinite(problem.predictors),'all')
    error("Data contain NaN or Inf.");
end

if ~isscalar(problem.residualScale) || problem.residualScale<=0 || ...
        ~isscalar(problem.objectiveScale) || problem.objectiveScale<=0
    error("Residual/objective scaling is invalid.");
end
% Verify that rational denominators remain separated from zero.
if any(problem.modelId == ["kirby2","hahn1"])
    [denLower,denUpper] = nist_denominator_interval(problem, ...
        problem.zLower,problem.zUpper);
    margin = min(abs([denLower;denUpper]),[],2);
    crosses = denLower<=0 & denUpper>=0;
    if any(crosses) || min(margin) < 1e-6
        error("The selected rational-model domain does not keep the denominator away from zero.");
    end
end
end
