function [cMinus,cPlus] = nist_interpolation_error_coefficients( ...
    problem,vertices,zLower,zUpper,nPoints)
%NIST_INTERPOLATION_ERROR_COEFFICIENTS Error coefficients at linear vertices.
problem = validate_nist_problem(problem);
vertices = double(vertices);
zLower = double(zLower(:)); zUpper = double(zUpper(:));
if nPoints < 2 || nPoints ~= floor(nPoints)
    error("nPoints must be an integer at least two.");
end
m = problem.nResidual; q = size(vertices,2); p = problem.nNonlinear;
h = (zUpper-zLower)/(nPoints-1);
cMinus = zeros(m,q); cPlus = zeros(m,q);
for k = 1:p
    [basisLower,basisUpper] = nist_basis_second_derivative_intervals( ...
        problem,zLower,zUpper,k);
    [valueLower,valueUpper] = localIntervalTimesConstantMatrix( ...
        basisLower,basisUpper,vertices);
    cMinus = cMinus+(h(k)^2/8)*min(0,-valueUpper);
    cPlus = cPlus+(h(k)^2/8)*max(0,-valueLower);
end
end

function [lower,upper] = localIntervalTimesConstantMatrix(aLower,aUpper,B)
[m,d] = size(aLower); q = size(B,2);
lower = zeros(m,q); upper = zeros(m,q);
for j = 1:d
    coefficient = B(j,:);
    positive = coefficient>=0;
    lo = aLower(:,j)*coefficient;
    hi = aUpper(:,j)*coefficient;
    if any(~positive)
        lo(:,~positive) = aUpper(:,j)*coefficient(~positive);
        hi(:,~positive) = aLower(:,j)*coefficient(~positive);
    end
    lower = lower+lo;
    upper = upper+hi;
end
end
