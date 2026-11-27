function [lower,upper] = nist_denominator_interval(problem,zLower,zUpper)
%NIST_DENOMINATOR_INTERVAL Interval enclosure of rational denominators.
zLower = double(zLower(:)); zUpper = double(zUpper(:));
x = problem.predictors(:,1);
switch problem.modelId
    case "kirby2"
        powers = [x,x.^2];
    case "hahn1"
        powers = [x,x.^2,x.^3];
    otherwise
        error("Problem %s is not a supported rational model.",problem.id);
end
lower = ones(size(x)); upper = ones(size(x));
for k = 1:size(powers,2)
    coefficient = powers(:,k);
    candidates = [coefficient*zLower(k),coefficient*zUpper(k)];
    lower = lower+min(candidates,[],2);
    upper = upper+max(candidates,[],2);
end
end
