function [d2g, d2l, d2w, ex] = ml_chain_bounds(G, Lr, R, Q, K)
%ML_CHAIN_BOUNDS Interval enclosures (vectorised over samples) of the second
%   derivatives of psi = sigma(t), t = gamma*L, L = ell + log x:
%     d2psi/dgamma2 = sigma''(t) L^2,   d2psi/dell2 = sigma''(t) gamma^2,
%     d2psi/dw_i^2  = sigma''(t) gamma^2 R_i^2 + sigma'(t) gamma Q_i.
%   G: 1x2 interval of gamma; Lr: n x 2 intervals of L; R, Q: cell arrays
%   (one n x 2 matrix per pole) of d log x/dw_i and d2 log x/dw_i^2.
%   The (gamma, L) box of each sample is split into K x K sub-boxes and the
%   enclosures are united (subdivision keeps the coupling between sigma(t)
%   and gamma, L). Valid for any K >= 1. Outputs: n x 2 (d2w: cell).
%   Optional ex: enclosures of psi (ex.psi) and of its first derivatives
%   ex.d1g = sigma'(t) L, ex.d1l = sigma'(t) gamma, ex.d1w{i} = sigma'(t) gamma R_i.
if nargin < 5 || isempty(K), K = 1; end
n = size(Lr, 1); m = numel(R);
d2g = [inf(n,1), -inf(n,1)]; d2l = d2g; d2w = repmat({d2g}, 1, m);
wantex = nargout > 3;
if wantex, ex.psi = d2g; ex.d1g = d2g; ex.d1l = d2g; ex.d1w = repmat({d2g}, 1, m); end
for a = 1:K
    Ga = G(1) + (G(2) - G(1)) * [(a-1)/K, a/K];
    G2 = [Ga(1)^2, Ga(2)^2];
    for b = 1:K
        Lb = [Lr(:,1) + (Lr(:,2) - Lr(:,1)) * (b-1)/K, Lr(:,1) + (Lr(:,2) - Lr(:,1)) * b/K];
        t = imul(repmat(Ga, n, 1), Lb);
        S1 = srange(t, 1); S2 = srange(t, 2);
        d2g = hull(d2g, imul(S2, isq(Lb)));
        d2l = hull(d2l, imul(S2, repmat(G2, n, 1)));
        if wantex
            sg = @(x) 1 ./ (1 + exp(-x));
            ex.psi = hull(ex.psi, [sg(t(:,1)), sg(t(:,2))]);   % sigma is increasing
            ex.d1g = hull(ex.d1g, imul(S1, Lb));
            ex.d1l = hull(ex.d1l, imul(S1, repmat(Ga, n, 1)));
            for i = 1:m
                ex.d1w{i} = hull(ex.d1w{i}, imul(imul(S1, repmat(Ga, n, 1)), R{i}));
            end
        end
        for i = 1:m
            v = iadd(imul(imul(S2, repmat(G2, n, 1)), isq(R{i})), imul(imul(S1, repmat(Ga, n, 1)), Q{i}));
            d2w{i} = hull(d2w{i}, v);
        end
    end
end
end

function h = hull(a, b)
h = [min(a(:,1), b(:,1)), max(a(:,2), b(:,2))];
end

function r = srange(t, order)
% range of sigma' (order 1) or sigma'' (order 2) on each interval t(i,:)
sig = @(x) 1 ./ (1 + exp(-x));
if order == 1
    f = @(x) sig(x) .* (1 - sig(x)); crit = 0;
else
    f = @(x) sig(x) .* (1 - sig(x)) .* (1 - 2*sig(x)); tc = log(2 + sqrt(3)); crit = [-tc, tc];
end
f1 = f(t(:,1)); f2 = f(t(:,2));
lo = min(f1, f2); hi = max(f1, f2);
for c = crit
    in = t(:,1) < c & t(:,2) > c;
    fc = f(c);
    lo(in) = min(lo(in), fc); hi(in) = max(hi(in), fc);
end
r = [lo, hi];
end

function r = imul(a, b)
P = [a(:,1).*b(:,1), a(:,1).*b(:,2), a(:,2).*b(:,1), a(:,2).*b(:,2)];
r = [min(P, [], 2), max(P, [], 2)];
end

function r = iadd(a, b)
r = a + b;
end

function r = isq(a)
s = [a(:,1).^2, a(:,2).^2];
r = [min(s, [], 2), max(s, [], 2)];
r(a(:,1) <= 0 & a(:,2) >= 0, 1) = 0;
end
