function [ggl, ggh, eel, eeh] = hill_inv_d2_bounds(s, g1, g2, E1, E2, K)
%HILL_INV_D2_BOUNDS Interval enclosures of the pure second derivatives of
%   c(gamma, Emax) = r^(1/gamma),  r = s/(Emax - s),  s = E0 - b >= 0,
%   over the box gamma in [g1,g2], Emax in [E1,E2] (requires E1 > max s):
%     d2c/dgamma^2 = c * L * (2*gamma + L) / gamma^4,   L = log r,
%     d2c/dEmax^2  = eta*(eta+1) * c / (Emax - s)^2,    eta = 1/gamma.
%   Returns per-sample lower/upper bounds (ggl,ggh) and (eel,eeh).
%   Samples with s == 0 have c == 0 identically and get zero bounds.
%
%   Method: natural interval extension, with the gamma and L ranges split
%   into K subintervals and the K^2 enclosures united (standard subdivision
%   to reduce the dependency effect). c = exp(L/gamma) is monotone in L and,
%   for fixed L, monotone in gamma, so its range on each sub-box is attained
%   at the corners. L*(2*gamma+L) is written as (L+gamma)^2 - gamma^2.
%   Floating-point rounding is not directed.
if nargin < 6 || isempty(K), K = 4; end
s = s(:);
n = numel(s);
ggl = zeros(n,1); ggh = zeros(n,1); eel = zeros(n,1); eeh = zeros(n,1);
p = s > 0;
if ~any(p), return; end
sp = s(p);
if any(E1 <= sp)
    error('hill_inv_d2_bounds:range', 'Requires Emax_min > E0 - min b.');
end
rlo = sp ./ (E2 - sp);
rhi = sp ./ (E1 - sp);
Llo = log(rlo);
Lhi = log(rhi);

% d2c/dEmax^2: product of positive, monotone factors
clo = min(rlo.^(1/g1), rlo.^(1/g2));
chi = max(rhi.^(1/g1), rhi.^(1/g2));
eel(p) = (1/g2)*(1/g2 + 1) * clo ./ (E2 - sp).^2;
eeh(p) = (1/g1)*(1/g1 + 1) * chi ./ (E1 - sp).^2;

% d2c/dgamma^2 with subdivision
gk = linspace(g1, g2, K+1);
lo = inf(size(sp)); hi = -inf(size(sp));
for i = 1:K
    a = gk(i); b = gk(i+1);
    for j = 1:K
        l1 = Llo + (Lhi - Llo)*(j-1)/K;
        l2 = Llo + (Lhi - Llo)*j/K;
        cc = [exp(l1/a), exp(l1/b), exp(l2/a), exp(l2/b)];
        cl = min(cc, [], 2); ch = max(cc, [], 2);
        t1 = l1 + a; t2 = l2 + b;                 % range of L + gamma
        sqh = max(t1.^2, t2.^2);
        sql = min(t1.^2, t2.^2);
        sql(t1 <= 0 & t2 >= 0) = 0;
        ql = sql - b^2; qh = sqh - a^2;           % range of L(2 gamma + L)
        [pl, ph] = imul(cl, ch, ql, qh);
        [pl, ph] = imul(pl, ph, 1/b^4, 1/a^4);
        lo = min(lo, pl); hi = max(hi, ph);
    end
end
ggl(p) = lo;
ggh(p) = hi;
end

function [lo, hi] = imul(alo, ahi, blo, bhi)
% interval product [alo,ahi]*[blo,bhi] (elementwise)
P = [alo.*blo, alo.*bhi, ahi.*blo, ahi.*bhi];
lo = min(P, [], 2);
hi = max(P, [], 2);
end
