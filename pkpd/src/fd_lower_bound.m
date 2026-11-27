function [LB, info] = fd_lower_bound(zmin, zmax, ge, V, Q, N, M, opts)
%FD_LOWER_BOUND m2(R,N,M) (m1 if M == N) for r(theta,z) = w(z).*(c(z) - Q theta),
%   theta in C = conv(V). For vertex v_a, with the fixed vector m_a = Q v_a,
%     r_dd = w_dd (c - m_a) + 2 w_d c_d + w c_dd   (d = gamma, Emax),
%   enclosed with LW_ENCLOSURES; error coefficients as in Proposition 1.
if ~isfield(opts, 'Ksub'), opts.Ksub = 4; end
zmin = zmin(:); zmax = zmax(:);
s = ge.E0 - ge.y(:); n = numel(s); p = size(V, 2);
W1 = [ones(1, p); V];
gv = linspace(zmin(1), zmax(1), N); Ev = linspace(zmin(2), zmax(2), N);
h = (zmax - zmin) / (N - 1); T = ((N - 1) / (M - 1))^2 * h.^2 / 8;
G = zeros(n, p * N^2); col = 0;
for iE = 1:N
    for ig = 1:N
        G(:, col + (1:p)) = fd_matrix(ge, [gv(ig); Ev(iE)], Q) * W1; col = col + p;
    end
end
X = lw_enclosures(s, zmin(1), zmax(1), zmin(2), zmax(2), opts.Ksub);
MV = Q * V;
Clo1 = zeros(n, p); Chi1 = zeros(n, p);
for a = 1:p
    D = X.C - [MV(:, a), MV(:, a)];
    rgg = imul(X.Wgg, D) + 2 * imul(X.Wg, X.Cg) + imul(X.W, X.Cgg);
    rEE = imul(X.WEE, D) + 2 * imul(X.WE, X.CE) + imul(X.W, X.CEE);
    Clo1(:, a) = T(1) * min(0, -rgg(:, 2)) + T(2) * min(0, -rEE(:, 2));
    Chi1(:, a) = T(1) * max(0, -rgg(:, 1)) + T(2) * max(0, -rEE(:, 1));
end
Clo = repmat(Clo1, 1, N^2); Chi = repmat(Chi1, 1, N^2);
if isfield(opts, 'no_solve') && opts.no_solve
    ok = false; val = NaN;
else
    [~, ~, ok, val] = norm_min_lp(G, Clo, Chi, '2');
end
if ok, LB = max(0, val - opts.lb_safety * (1 + abs(val))); else, LB = -inf; end
info.ncols = size(G, 2);
if isfield(opts, 'return_matrices') && opts.return_matrices
    info.G = G; info.Clo = Clo; info.Chi = Chi; info.gv = gv; info.Ev = Ev;
end
end

function r = imul(a, b)
P = [a(:,1).*b(:,1), a(:,1).*b(:,2), a(:,2).*b(:,1), a(:,2).*b(:,2)];
r = [min(P, [], 2), max(P, [], 2)];
end
