function [lo, hi, ex1] = ml_d2_bounds(u, zmin, zmax, K)
%ML_D2_BOUNDS Interval enclosures, over the box [zmin, zmax] in
%   z = (gamma, ell, w), of the pure second derivatives of psi(k; z)
%   (columns: d2/dgamma2, d2/dell2, d2/dw2), for every sample k.
%
%   psi = sigma(t), t = gamma*L, L = ell + log x(k; w), sigma = logistic.
%     d2psi/dgamma2 = sigma''(t) L^2
%     d2psi/dell2   = sigma''(t) gamma^2
%     d2psi/dw2     = sigma''(t) gamma^2 R1^2 + sigma'(t) gamma Q,
%   R1 = d log x/dw, Q = d2 log x/dw2.
%   With u >= 0, x = sum_j u_j h_j(w), h_j = j e^{-2w}(1-e^{-w})^{j-1} > 0:
%     x in [sum u_j min h_j, sum u_j max h_j]   (exact range of each term),
%     R1 = weighted mean of a_j = -2 + (j-1)/(e^w - 1),
%     Q  = weighted mean of a_j' + weighted variance of a_j, the variance
%          being at most (range of a_j at fixed w)^2/4 (Popoviciu).
%   Rounding is not directed.
if nargin < 4 || isempty(K), K = 4; end   % (gamma, L) subdivision, ML_CHAIN_BOUNDS
u = u(:); n = numel(u);
if any(u < 0), error('ml_d2_bounds:input', 'Requires a nonnegative input.'); end
g1 = zmin(1); g2 = zmax(1); l1 = zmin(2); l2 = zmax(2); w1 = zmin(3); w2 = zmax(3);
s1 = 1/(exp(w1) - 1); s2 = 1/(exp(w2) - 1);           % s1 >= s2
phi1 = exp(w1)*s1^2; phi2 = exp(w2)*s2^2;            % e^w/(e^w-1)^2, decreasing
G = [g1, g2]; G2 = [g1^2, g2^2];
[sp1, sp2] = sigma_ranges_init();
lo = zeros(n, 3); hi = zeros(n, 3);
ex1.psi = zeros(n, 2); ex1.d1 = zeros(n, 3, 2);   % psi and first derivatives (zero where x = 0)
La = nan(n, 2); R1a = nan(n, 2); Qa = R1a;
for kk = 1:n
    k = kk - 1;                      % time index
    j = (1:k-1)';                    % x(k) uses u(k-1-j), j >= 1
    if isempty(j), continue; end
    ui = u(k - j);                   % u(k-1-j) -> MATLAB index k-j
    act = ui > 0;
    if ~any(act), continue; end
    j = j(act); ua = ui(act);
    h = @(w) j .* exp(-2*w) .* (1 - exp(-w)).^(j - 1);
    h1 = h(w1); h2 = h(w2);
    hmin = min(h1, h2); hmax = max(h1, h2);
    ws = log((j + 1) / 2);           % stationary point of h_j
    in = ws > w1 & ws < w2;
    if any(in)
        hs = h(ws);
        hmax(in) = max(hmax(in), hs(in));
    end
    xlo = ua' * hmin; xhi = ua' * hmax;
    R1 = [-2 + (min(j) - 1)*s2, -2 + (max(j) - 1)*s1];
    Q = [-(max(j) - 1)*phi1, -(min(j) - 1)*phi2 + ((max(j) - min(j))*s1)^2/4];
    La(kk, :) = [l1 + log(xlo), l2 + log(xhi)];
    R1a(kk, :) = R1; Qa(kk, :) = Q;
end
act = find(~isnan(La(:, 1)));
if isempty(act), return; end
[d2g, d2l, d2w, ex] = ml_chain_bounds(G, La(act, :), {R1a(act, :)}, {Qa(act, :)}, K);
ex1.psi(act, :) = ex.psi;
D = {ex.d1g, ex.d1l, ex.d1w{1}};
for j = 1:3, ex1.d1(act, j, 1) = D{j}(:, 1); ex1.d1(act, j, 2) = D{j}(:, 2); end
lo(act, :) = [d2g(:,1), d2l(:,1), d2w{1}(:,1)];
hi(act, :) = [d2g(:,2), d2l(:,2), d2w{1}(:,2)];
end

function [sp1, sp2] = sigma_ranges_init()
sig = @(t) 1 ./ (1 + exp(-t));
f1 = @(t) sig(t) .* (1 - sig(t));                    % sigma'
f2 = @(t) f1(t) .* (1 - 2*sig(t));                   % sigma''
tc = log(2 + sqrt(3));                               % extrema of sigma''
sp1 = @(t) range_with(f1, t, 0);
sp2 = @(t) range_with(f2, t, [-tc, tc]);
end

function r = range_with(f, t, crit)
c = crit(crit > t(1) & crit < t(2));
v = f([t(1), t(2), c]);
r = [min(v), max(v)];
end

function r = imul(a, b)
v = [a(1)*b(1), a(1)*b(2), a(2)*b(1), a(2)*b(2)];
r = [min(v), max(v)];
end

function r = iadd(a, b)
r = [a(1) + b(1), a(2) + b(2)];
end

function r = isq(a)
v = [a(1)^2, a(2)^2];
if a(1) <= 0 && a(2) >= 0, r = [0, max(v)]; else, r = [min(v), max(v)]; end
end
