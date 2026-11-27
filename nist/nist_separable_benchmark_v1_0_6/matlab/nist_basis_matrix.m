function Phi = nist_basis_matrix(problem,z)
%NIST_BASIS_MATRIX Evaluate the model basis matrix Phi(z).
problem = validate_nist_problem(problem);
z = double(z(:));
if numel(z) ~= problem.nNonlinear
    error("The nonlinear parameter vector has the wrong length.");
end
X = problem.predictors;
x = X(:,1);
switch problem.modelId
    case "boxbod"
        Phi = 1-exp(-z(1)*x);
    case "mgh17"
        Phi = [ones(size(x)),exp(-z(1)*x),exp(-z(2)*x)];
    case "lanczos2"
        Phi = [exp(-z(1)*x),exp(-z(2)*x),exp(-z(3)*x)];
    case "kirby2"
        D = 1+z(1)*x+z(2)*x.^2;
        Phi = [ones(size(x)),x,x.^2]./D;
    case "hahn1"
        D = 1+z(1)*x+z(2)*x.^2+z(3)*x.^3;
        Phi = [ones(size(x)),x,x.^2,x.^3]./D;
    case "gauss1"
        Phi = [exp(-z(1)*x), ...
            exp(-((x-z(2))/z(3)).^2), ...
            exp(-((x-z(4))/z(5)).^2)];
    case "nelson"
        x1 = X(:,1); x2 = X(:,2);
        Phi = [ones(size(x1)),-x1.*exp(-z(1)*x2)];
    case "mgh10"
        Phi = exp(z(1)./(x+z(2)));
    otherwise
        error("Unsupported NIST model '%s'.",problem.modelId);
end
if isfield(problem,"residualScale")
    Phi = Phi/double(problem.residualScale);
end
if ~isequal(size(Phi),[problem.nResidual,problem.nLinear]) || any(~isfinite(Phi),'all')
    error("Basis evaluation failed for %s.",problem.id);
end
end
