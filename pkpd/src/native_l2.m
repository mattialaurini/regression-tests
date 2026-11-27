function [val, v, ok, cert, decided] = native_l2(G, Clo, Chi, opts)
%NATIVE_L2 Structured solver for  min_{v in simplex, Clo v <= e <= Chi v} ||G v + e||_2.
%   Eliminating e row by row gives the smooth convex problem
%     f(v) = sum_k max(0, P_k v, -Q_k v)^2,  P = G + Clo, Q = G + Chi,
%   on the unit simplex (v only: no auxiliary variables, no constraints but the
%   simplex). Solved by accelerated projected gradient (FISTA with restarts),
%   started at the best vertex; stopped when the certified gap
%     sqrt(f(v)) - cert <= rtol * sqrt(f(v))
%   is reached, where cert is the weak-duality lower bound of NORM_MIN_LP computed
%   from the residual direction (valid for ANY v). val = sqrt(f(v)) is the exact
%   objective at a feasible point; NORM_MIN_LP then polishes v (active set) and
%   recomputes the certificate.
%   opts.tau (optional): stop as soon as the relaxation value is known to be
%   > tau (cert > tau) or <= tau (sqrt(f) <= tau); decided = true then.
if nargin < 4, opts = struct(); end
if ~isfield(opts, 'rtol'), opts.rtol = 1e-6; end
if ~isfield(opts, 'maxit'), opts.maxit = 3000; end
if ~isfield(opts, 'tau'), opts.tau = NaN; end
decided = false;
[m, r] = size(G);
if isempty(Clo), P = G; Q = G; Cl = zeros(m, r); Ch = Cl; else, P = G + Clo; Q = G + Chi; Cl = Clo; Ch = Chi; end
fval = @(v) sum(max(0, max(P * v, -(Q * v))).^2);
% start: best vertex (all columns at once)
fv = sum(max(0, max(P, -Q)).^2, 1);
[~, j] = min(fv); v = zeros(r, 1); v(j) = 1;
% step: Lipschitz constant of the gradient <= 2 (||P||^2 + ||Q||^2) (power iteration)
Lp = specnorm2(P) + specnorm2(Q); step = 1 / (2 * Lp);
y = v; t = 1; fbest = fval(v); vbest = v; cert = 0;
for it = 1:opts.maxit
    pv = P * y; qv = Q * y;
    g = 2 * (P' * max(pv, 0) + Q' * min(qv, 0));      % gradient (rows with P y > 0 or Q y < 0)
    vn = simplex_proj(y - step * g);
    fn = fval(vn);
    if fn > fval(v), t = 1; y = v; continue; end      % restart (monotone variant)
    tn = (1 + sqrt(1 + 4 * t^2)) / 2;
    y = vn + ((t - 1) / tn) * (vn - v); v = vn; t = tn;
    if fn < fbest, fbest = fn; vbest = vn; end
    if mod(it, 10) == 0 || it == opts.maxit
        cert = max(cert, certify(G, Cl, Ch, P, Q, vbest));
        if sqrt(fbest) - cert <= opts.rtol * max(sqrt(fbest), eps), break; end
        % decision w.r.t. the pruning threshold: the relaxation value lies in
        % [cert, sqrt(fbest)]; stop once it is known to be > tau or <= tau
        if ~isnan(opts.tau) && (cert > opts.tau || sqrt(fbest) <= opts.tau), decided = true; break; end
    end
end
v = vbest; val = sqrt(fbest); ok = true;
cert = max(cert, certify(G, Cl, Ch, P, Q, v));
if ~isnan(opts.tau) && (cert > opts.tau || val <= opts.tau), decided = true; end
end

function c = certify(G, Cl, Ch, P, Q, v)
% weak duality: lambda = w/||w||, w = signed residual with the best error vector
w = max(P * v, 0) + min(Q * v, 0);
if norm(w) == 0, c = 0; return; end
l = w / norm(w);
c = max(0, min(l' * G + sum(min(l .* Cl, l .* Ch), 1)));
end

function s = specnorm2(A)
% ||A||_2^2 by a few power iterations (upper estimate scaled by 1.05)
x = ones(size(A, 2), 1) / sqrt(size(A, 2)); s = 0;
for i = 1:30
    y = A' * (A * x); s = norm(y); if s == 0, break; end
    x = y / s;
end
s = 1.05 * s;
end
