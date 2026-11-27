function [lower,upper] = nist_basis_second_derivative_intervals( ...
    problem,zLower,zUpper,coordinate)
%NIST_BASIS_SECOND_DERIVATIVE_INTERVALS Certified pure-second-derivative bounds.
%
% Returns componentwise bounds for d^2 Phi(z) / dz_coordinate^2 over
% zLower <= z <= zUpper. Mixed derivatives are not required by the
% multilinear interpolation enclosure used in the manuscript.
problem = validate_nist_problem(problem);
zLower = double(zLower(:)); zUpper = double(zUpper(:));
if numel(zLower) ~= problem.nNonlinear || any(zLower>zUpper)
    error("Invalid nonlinear box.");
end
if coordinate < 1 || coordinate > problem.nNonlinear
    error("Invalid nonlinear coordinate.");
end
m = problem.nResidual; d = problem.nLinear;
lower = zeros(m,d); upper = zeros(m,d);
X = problem.predictors; x = X(:,1);

switch problem.modelId
    case "boxbod"
        E = ivExp(ivScale([-zUpper(1),-zLower(1)],x));
        [lower(:,1),upper(:,1)] = ivScalePair(E.lo,E.hi,-x.^2);

    case "mgh17"
        column = coordinate+1;
        E = ivExp(ivScale([-zUpper(coordinate),-zLower(coordinate)],x));
        lower(:,column) = x.^2.*E.lo;
        upper(:,column) = x.^2.*E.hi;

    case "lanczos2"
        column = coordinate;
        E = ivExp(ivScale([-zUpper(coordinate),-zLower(coordinate)],x));
        lower(:,column) = x.^2.*E.lo;
        upper(:,column) = x.^2.*E.hi;

    case "nelson"
        x1 = X(:,1); x2 = X(:,2);
        E = ivExp(ivScale([-zUpper(1),-zLower(1)],x2));
        [lower(:,2),upper(:,2)] = ivScalePair(E.lo,E.hi,-x1.*x2.^2);

    case {"kirby2","hahn1"}
        [Dlo,Dhi] = nist_denominator_interval(problem,zLower,zUpper);
        if any(Dlo<=0)
            error("Rational denominator is not strictly positive on the box.");
        end
        denominatorLo = Dlo.^3;
        denominatorHi = Dhi.^3;
        q = x.^coordinate;
        if problem.modelId == "kirby2"
            powers = [ones(size(x)),x,x.^2];
        else
            powers = [ones(size(x)),x,x.^2,x.^3];
        end
        for j = 1:d
            numerator = 2*powers(:,j).*q.^2;
            candidates = [numerator./denominatorLo,numerator./denominatorHi];
            lower(:,j) = min(candidates,[],2);
            upper(:,j) = max(candidates,[],2);
        end

    case "mgh10"
        a = iv(zLower(1),zUpper(1));
        den = iv(x+zLower(2),x+zUpper(2));
        if any(den.lo<=0)
            error("MGH10 denominator crosses zero.");
        end
        phi = ivExp(ivDiv(a,den));
        if coordinate == 1
            term = ivDiv(iv(1,1),ivPowPositive(den,2));
        else
            term1 = ivDiv(ivSquare(a),ivPowPositive(den,4));
            term2 = ivScaleInterval(ivDiv(a,ivPowPositive(den,3)),2);
            term = ivAdd(term1,term2);
        end
        value = ivMul(phi,term);
        lower(:,1) = value.lo;
        upper(:,1) = value.hi;

    case "gauss1"
        if coordinate == 1
            % The common residual scaling is applied below, together with
            % all other models and Gauss1 coordinates. Do not return here.
            E = ivExp(ivScale([-zUpper(1),-zLower(1)],x));
            lower(:,1) = x.^2.*E.lo;
            upper(:,1) = x.^2.*E.hi;
        else
            if coordinate <= 3
                centerIndex = 2; widthIndex = 3; column = 2;
            else
                centerIndex = 4; widthIndex = 5; column = 3;
            end
            t = iv(x-zUpper(centerIndex),x-zLower(centerIndex));
            w = iv(zLower(widthIndex),zUpper(widthIndex));
            if any(w.lo<=0)
                error("Gaussian width must remain positive.");
            end
            t2 = ivSquare(t);
            w2 = ivPowPositive(w,2);
            a = ivDiv(t2,w2);
            phi = ivExp(ivNeg(a));
            if coordinate == centerIndex
                term = ivSub( ...
                    ivScaleInterval(ivDiv(t2,ivPowPositive(w,4)),4), ...
                    ivScaleInterval(ivDiv(iv(1,1),w2),2));
            elseif coordinate == widthIndex
                t4 = ivSquare(t2);
                term = ivSub( ...
                    ivScaleInterval(ivDiv(t4,ivPowPositive(w,6)),4), ...
                    ivScaleInterval(ivDiv(t2,ivPowPositive(w,4)),6));
            else
                error("Internal Gaussian coordinate mapping error.");
            end
            value = ivMul(phi,term);
            lower(:,column) = value.lo;
            upper(:,column) = value.hi;
        end

    otherwise
        error("Second-derivative bounds are not implemented for %s.",problem.id);
end

if isfield(problem,"residualScale")
    lower = lower/double(problem.residualScale);
    upper = upper/double(problem.residualScale);
end
% Protect against roundoff-induced reversed endpoints.
lo = min(lower,upper); hi = max(lower,upper);
lower = lo; upper = hi;
end

function result = iv(lo,hi)
result.lo = double(lo); result.hi = double(hi);
lo2 = min(result.lo,result.hi); hi2 = max(result.lo,result.hi);
result.lo = lo2; result.hi = hi2;
end

function result = ivAdd(a,b)
result = iv(a.lo+b.lo,a.hi+b.hi);
end

function result = ivSub(a,b)
result = iv(a.lo-b.hi,a.hi-b.lo);
end

function result = ivNeg(a)
result = iv(-a.hi,-a.lo);
end

function result = ivMul(a,b)
values = cat(ndims(a.lo)+1,a.lo.*b.lo,a.lo.*b.hi,a.hi.*b.lo,a.hi.*b.hi);
result = iv(min(values,[],ndims(values)),max(values,[],ndims(values)));
end

function result = ivDiv(a,b)
if any(b.lo<=0 & b.hi>=0,'all')
    error("Interval division by an interval containing zero.");
end
reciprocal = iv(1./b.hi,1./b.lo);
result = ivMul(a,reciprocal);
end

function result = ivExp(a)
result = iv(exp(a.lo),exp(a.hi));
end

function result = ivSquare(a)
lo = zeros(size(a.lo));
hi = max(a.lo.^2,a.hi.^2);
positive = a.lo>=0;
negative = a.hi<=0;
lo(positive) = a.lo(positive).^2;
lo(negative) = a.hi(negative).^2;
result = iv(lo,hi);
end

function result = ivPowPositive(a,power)
if any(a.lo<0,'all')
    error("ivPowPositive requires a nonnegative interval.");
end
result = iv(a.lo.^power,a.hi.^power);
end

function result = ivScaleInterval(a,scale)
if scale>=0
    result = iv(scale*a.lo,scale*a.hi);
else
    result = iv(scale*a.hi,scale*a.lo);
end
end

function result = ivScale(endpointPair,constant)
% endpointPair is [lower upper] for a scalar interval; constant is an array.
lo = endpointPair(1).*constant;
hi = endpointPair(2).*constant;
negative = constant<0;
tmp = lo(negative);
lo(negative) = hi(negative);
hi(negative) = tmp;
result = iv(lo,hi);
end

function [lo,hi] = ivScalePair(aLo,aHi,constant)
positive = constant>=0;
lo = constant.*aLo; hi = constant.*aHi;
lo(~positive) = constant(~positive).*aHi(~positive);
hi(~positive) = constant(~positive).*aLo(~positive);
end
