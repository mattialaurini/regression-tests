function [lo, hi, ex1] = ml_d2_bounds_oc(u, zmin, zmax, K)
%ML_D2_BOUNDS_OC Interval enclosures over the box [zmin,zmax] in
%   z = (gamma, ell, w1, w2) of the pure second derivatives of psi(k; z) for
%   the one-compartment + effect-site model (columns: gamma, ell, w1, w2).
%
%   x(k) = sum_{i,l} u(k-2-i-l) f_i(w1) f_l(w2),  f_i(w) = e^{-w}(1-e^{-w})^i > 0.
%   * range of x: exact range of each f_i, combined by convolution;
%   * d log x/dw1 = -1 + s(w1) E[i],  s(w) = 1/(e^w - 1),
%     d2 log x/dw1^2 = -phi(w1) E[i] + s(w1)^2 Var(i),  phi(w) = e^w/(e^w-1)^2,
%     E, Var over the positive weights u(k-2-n) p1^i p2^l, n = i + l.
%     Given n, i is truncated geometric on {0..n} with ratio rho = p1/p2:
%     its mean m_n(rho) increases in rho and n, its variance v_n(rho) increases
%     in n and is maximal at rho = 1 (symmetric in log rho). Hence, over the box,
%       E[i]   in [m_nmin(rho_lo), m_nmax(rho_hi)]
%       Var(i) <= v_nmax(rho*) + (m_nmax(rho_hi) - m_nmin(rho_lo))^2/4
%     (rho* = point of [rho_lo, rho_hi] closest to 1; total variance = mean of
%     the conditional variances + variance of the conditional means).
%     nmin, nmax = smallest/largest lag with u(k-2-n) > 0. Same for w2 with
%     ratio p2/p1.
%   * psi = sigma(gamma(ell + log x)): chain rule as in ML_D2_BOUNDS.
%   Requires u >= 0. Rounding is not directed.
if nargin < 4 || isempty(K), K = 4; end   % (gamma, L) subdivision, ML_CHAIN_BOUNDS
u = u(:); n = numel(u);
if any(u < 0), error('ml_d2_bounds_oc:input', 'Requires a nonnegative input.'); end
g1 = zmin(1); g2 = zmax(1); l1 = zmin(2); l2 = zmax(2);
L = max(n - 2, 1);                          % maximal lag i+l
i = (0:L)';
[F1lo, F1hi] = frange(i, zmin(3), zmax(3));
[F2lo, F2hi] = frange(i, zmin(4), zmax(4));
glo = conv(F1lo, F2lo); ghi = conv(F1hi, F2hi);  % range of sum_i f_i(w1) f_{m-i}(w2)
glo = glo(1:L+1); ghi = ghi(1:L+1);
cvlo = conv(u, glo); cvhi = conv(u, ghi);
qfirst = find(u > 0, 1) - 1;                % first sample with positive input
s = @(w) 1 ./ (exp(w) - 1); phi = @(w) exp(w) ./ (exp(w) - 1).^2;
pp = @(w) 1 - exp(-w);
rl = pp(zmin(3)) / pp(zmax(4)); rh = pp(zmax(3)) / pp(zmin(4));   % rho = p1/p2
[m1lo, ~] = tgeom(rl, L); [m1hi, ~] = tgeom(rh, L); [~, v1s] = tgeom(closest1(rl, rh), L);
[m2lo, ~] = tgeom(1/rh, L); [m2hi, ~] = tgeom(1/rl, L); [~, v2s] = tgeom(closest1(1/rh, 1/rl), L);
sA = [s(zmax(3)), s(zmin(3))]; fA = [phi(zmax(3)), phi(zmin(3))];
sB = [s(zmax(4)), s(zmin(4))]; fB = [phi(zmax(4)), phi(zmin(4))];
G = [g1, g2]; G2 = [g1^2, g2^2];
lo = zeros(n, 4); hi = zeros(n, 4);
ex1.psi = zeros(n, 2); ex1.d1 = zeros(n, 4, 2);   % psi and first derivatives (zero where x = 0)
if isempty(qfirst), return; end
Lr = nan(n, 2); R1a = nan(n, 2); Q1a = R1a; R2a = R1a; Q2a = R1a;
pos = find(u > 0) - 1;                      % samples with positive input
for kk = 3:n                                % k = kk-1 >= 2
    k = kk - 1;
    q = pos(pos <= k - 2);
    if isempty(q), continue; end
    nmax = k - 2 - min(q); nmin = k - 2 - max(q);
    xlo = cvlo(kk - 2); xhi = cvhi(kk - 2);
    if ~(xlo > 0), continue; end
    [R1, Q1] = rq(m1lo(nmin+1), m1hi(nmax+1), v1s(nmax+1), sA, fA);
    [R2, Q2] = rq(m2lo(nmin+1), m2hi(nmax+1), v2s(nmax+1), sB, fB);
    Lr(kk, :) = [l1 + log(xlo), l2 + log(xhi)];
    R1a(kk, :) = R1; Q1a(kk, :) = Q1; R2a(kk, :) = R2; Q2a(kk, :) = Q2;
end
act = find(~isnan(Lr(:, 1)));
if isempty(act), return; end
[d2g, d2l, d2w, ex] = ml_chain_bounds(G, Lr(act, :), {R1a(act, :), R2a(act, :)}, ...
    {Q1a(act, :), Q2a(act, :)}, K);
ex1.psi(act, :) = ex.psi;
D = {ex.d1g, ex.d1l, ex.d1w{1}, ex.d1w{2}};
for j = 1:4, ex1.d1(act, j, 1) = D{j}(:, 1); ex1.d1(act, j, 2) = D{j}(:, 2); end
lo(act, :) = [d2g(:,1), d2l(:,1), d2w{1}(:,1), d2w{2}(:,1)];
hi(act, :) = [d2g(:,2), d2l(:,2), d2w{1}(:,2), d2w{2}(:,2)];
end

function [R, Q] = rq(Elo, Ehi, Vc, sr, fr)
% R = -1 + s E, Q = -phi E + s^2 Var, with E in [Elo,Ehi] (>= 0),
% Var in [0, Vc + (Ehi-Elo)^2/4], s in sr, phi in fr (positive intervals)
Vhi = Vc + (Ehi - Elo)^2 / 4;
R = [-1 + sr(1) * Elo, -1 + sr(2) * Ehi];
Q = [-fr(2) * Ehi, -fr(1) * Elo + sr(2)^2 * Vhi];
end

function [m, v] = tgeom(rho, L)
% mean and variance of the truncated geometric law P(i) ~ rho^i, i = 0..n,
% for n = 0..L (vectors of length L+1)
i = (0:L)';
lw = i * log(rho); lw = lw - max(lw);        % scaled weights for stability
w = exp(lw);
S0 = cumsum(w); S1 = cumsum(i .* w); S2 = cumsum(i.^2 .* w);
m = S1 ./ S0; v = max(S2 ./ S0 - m.^2, 0);
end

function r = closest1(a, b)
if a <= 1 && b >= 1, r = 1; elseif b < 1, r = b; else, r = a; end
end

function [flo, fhi] = frange(i, wa, wb)
% exact range of f_i(w) = e^{-w}(1-e^{-w})^i on [wa,wb]; maximum at w = log(i+1)
f = @(w) exp(-w) .* (1 - exp(-w)).^i;
fa = f(wa); fb = f(wb);
flo = min(fa, fb); fhi = max(fa, fb);
ws = log(i + 1); in = ws > wa & ws < wb;
if any(in)
    fs = exp(-ws) .* (1 - exp(-ws)).^i;
    fhi(in) = max(fhi(in), fs(in));
end
end

function r = srange(t, order)
sig = @(x) 1 ./ (1 + exp(-x));
if order == 1
    f = @(x) sig(x) .* (1 - sig(x)); crit = 0;
else
    f = @(x) sig(x) .* (1 - sig(x)) .* (1 - 2*sig(x)); tc = log(2 + sqrt(3)); crit = [-tc, tc];
end
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
