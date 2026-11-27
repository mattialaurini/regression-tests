function D3 = ml_d3_bounds(u, zmin, zmax)
%ML_D3_BOUNDS Bounds on |d3 psi / dz_d^2 dz_j| over the box [zmin, zmax], for
%   every sample k (output n x d x d, entry (k,d,j)).
%   psi = sigma(t), t = gamma (ell + log x):
%     psi_ddj = sigma''' t_d^2 t_j + sigma'' (2 t_dj t_d + t_dd t_j) + sigma' t_ddj,
%   sigma''' = sigma'(1 - 6 sigma + 6 sigma^2) (extrema at t = 0, +-2.2924).
%   t derivatives: t_g = L, t_l = gamma, t_wi = gamma R_i; t_gl = 1, t_gwi = R_i,
%   t_wiwi = gamma Q_ii, t_w1w2 = gamma Q_12; t_wiwi g = Q_ii, t_wiwi wj = gamma S_iij
%   (all other second/third derivatives of t vanish).
%   R, Q, S: derivatives of log x, i.e. cumulants of the lag weights:
%     d log x = E[a], d2 = E[a'] + Var(a), d3 = E[a''] + 3 Cov(a,a') + k3(a),
%     mixed: Q_12 = Cov(a1,a2), S_112 = Cov(a1',a2) + E[(a1-m1)^2 (a2-m2)].
%   Third central moments: |k3| <= min(range*Var, range^3/(6 sqrt 3)).
%   d = 3: double pole (lag j, a_j = -2 + (j-1) s(w)); d = 4: one compartment +
%   effect site (lags i, l; truncated-geometric bounds as in ML_D2_BOUNDS_OC).
%   Rounding is not directed.
u = u(:); n = numel(u); d = numel(zmin);
g1 = zmin(1); g2 = zmax(1); G = [g1, g2]; l1 = zmin(2); l2 = zmax(2);
s = @(w) 1 ./ (exp(w) - 1); phi = @(w) exp(w) ./ (exp(w) - 1).^2;
dphi = @(w) exp(w) .* (exp(w) + 1) ./ (exp(w) - 1).^3;          % |phi'|, decreasing
c6 = 1 / (6 * sqrt(3));
D3 = zeros(n, d, d);
if d == 3
    wa = zmin(3); wb = zmax(3);
    sr = [s(wb), s(wa)]; fr = [phi(wb), phi(wa)]; dfr = [dphi(wb), dphi(wa)];
    LA = nan(n, 2); RA = LA; QA = LA; SA = LA;
    for kk = 1:n
        k = kk - 1; j = (1:k-1)';
        if isempty(j), continue; end
        ui = u(k - j); act = ui > 0;
        if ~any(act), continue; end
        j = j(act); ua = ui(act);
        h = @(w) j .* exp(-2*w) .* (1 - exp(-w)).^(j - 1);
        ha = h(wa); hb = h(wb); hmin = min(ha, hb); hmax = max(ha, hb);
        ws = log((j + 1) / 2); in = ws > wa & ws < wb;
        if any(in), hs = h(ws); hmax(in) = max(hmax(in), hs(in)); end
        xlo = ua' * hmin; xhi = ua' * hmax;
        J0 = min(j) - 1; J1 = max(j) - 1; Dl = J1 - J0;
        R = [-2 + J0 * sr(1), -2 + J1 * sr(2)];
        Var = Dl^2 / 4;                                    % Var of J (Popoviciu)
        Q = [-J1 * fr(2), -J0 * fr(1) + sr(2)^2 * Var];
        K3 = min(Dl * Var, Dl^3 * c6);
        S = [dfr(1) * J0 - 3 * fr(2) * sr(2) * Var - sr(2)^3 * K3, dfr(2) * J1 + sr(2)^3 * K3];
        LA(kk, :) = [l1 + log(xlo), l2 + log(xhi)]; RA(kk, :) = R; QA(kk, :) = Q; SA(kk, :) = S;
    end
    act = ~isnan(LA(:, 1));
    if any(act)
        D3(act, :, :) = assemble_v(G, LA(act, :), {RA(act, :)}, {QA(act, :)}, [], {SA(act, :)}, {}, 3);
    end
else
    wa1 = zmin(3); wb1 = zmax(3); wa2 = zmin(4); wb2 = zmax(4);
    Lg = max(n - 2, 1); i = (0:Lg)';
    [F1lo, F1hi] = frange(i, wa1, wb1); [F2lo, F2hi] = frange(i, wa2, wb2);
    glo = conv(F1lo, F2lo); ghi = conv(F1hi, F2hi); glo = glo(1:Lg+1); ghi = ghi(1:Lg+1);
    cvlo = conv(u, glo); cvhi = conv(u, ghi);
    pp = @(w) 1 - exp(-w);
    rl = pp(wa1) / pp(wb2); rh = pp(wb1) / pp(wa2);
    [m1lo, ~] = tgeom(rl, Lg); [m1hi, ~] = tgeom(rh, Lg); [~, v1s] = tgeom(closest1(rl, rh), Lg);
    [m2lo, ~] = tgeom(1/rh, Lg); [m2hi, ~] = tgeom(1/rl, Lg); [~, v2s] = tgeom(closest1(1/rh, 1/rl), Lg);
    sA = [s(wb1), s(wa1)]; fA = [phi(wb1), phi(wa1)]; dA = [dphi(wb1), dphi(wa1)];
    sB = [s(wb2), s(wa2)]; fB = [phi(wb2), phi(wa2)]; dB = [dphi(wb2), dphi(wa2)];
    pos = find(u > 0) - 1;
    LA = nan(n, 2); R1A = LA; R2A = LA; Q11A = LA; Q22A = LA; Q12A = LA; S1A = LA; S2A = LA; S12A = LA; S21A = LA;
    for kk = 3:n
        k = kk - 1; q = pos(pos <= k - 2);
        if isempty(q), continue; end
        nmax = k - 2 - min(q); nmin = k - 2 - max(q);
        xlo = cvlo(kk - 2); xhi = cvhi(kk - 2);
        if ~(xlo > 0), continue; end
        [R1, Q11, S111, Vi, Ei] = rqs(m1lo(nmin+1), m1hi(nmax+1), v1s(nmax+1), sA, fA, dA, nmax, c6);
        [R2, Q22, S222, Vl, El] = rqs(m2lo(nmin+1), m2hi(nmax+1), v2s(nmax+1), sB, fB, dB, nmax, c6);
        cv = sqrt(Vi * Vl);                                % |Cov(i,l)|
        Q12 = sA(2) * sB(2) * cv * [-1, 1];
        S112 = (fA(2) * sB(2) * cv + sA(2)^2 * sB(2) * nmax * Vi) * [-1, 1];
        S221 = (fB(2) * sA(2) * cv + sB(2)^2 * sA(2) * nmax * Vl) * [-1, 1];
        LA(kk, :) = [l1 + log(xlo), l2 + log(xhi)];
        R1A(kk, :) = R1; R2A(kk, :) = R2; Q11A(kk, :) = Q11; Q22A(kk, :) = Q22; Q12A(kk, :) = Q12;
        S1A(kk, :) = S111; S2A(kk, :) = S222; S12A(kk, :) = S112; S21A(kk, :) = S221;
    end
    act = ~isnan(LA(:, 1));
    if any(act)
        D3(act, :, :) = assemble_v(G, LA(act, :), {R1A(act, :), R2A(act, :)}, {Q11A(act, :), Q22A(act, :)}, ...
            Q12A(act, :), {S1A(act, :), S2A(act, :)}, {S12A(act, :), S21A(act, :)}, 4);
    end
end
end

function M = assemble_v(G, L, R, Q, Q12, S, Sx, d)
% vectorised over rows: all intervals are n x 2
n = size(L, 1); Gv = repmat(G, n, 1); one = ones(n, 2); Z = zeros(n, 2);
t = vmul(Gv, L);
sg = @(x) 1 ./ (1 + exp(-x));
f1 = @(x) sg(x) .* (1 - sg(x)); f2 = @(x) f1(x) .* (1 - 2*sg(x));
f3 = @(x) f1(x) .* (1 - 6*sg(x) + 6*sg(x).^2);
S1 = vrng(f1, t, 0); S2 = vrng(f2, t, log(2 + sqrt(3)) * [-1, 1]); S3 = vrng(f3, t, [0, -2.2924, 2.2924]);
m = numel(R);
td = cell(1, d); td{1} = L; td{2} = Gv;
for i = 1:m, td{2+i} = vmul(Gv, R{i}); end
tdj = repmat({Z}, d, d);
tdj{1,2} = one; tdj{2,1} = one;
for i = 1:m
    tdj{1,2+i} = R{i}; tdj{2+i,1} = R{i}; tdj{2+i,2+i} = vmul(Gv, Q{i});
end
if m == 2, tdj{3,4} = vmul(Gv, Q12); tdj{4,3} = tdj{3,4}; end
tddj = repmat({Z}, d, d);
for i = 1:m
    tddj{2+i,1} = Q{i}; tddj{2+i,2+i} = vmul(Gv, S{i});
    if m == 2, tddj{2+i,2+(3-i)} = vmul(Gv, Sx{i}); end
end
M = zeros(n, d, d);
for dd = 1:d
    for jj = 1:d
        v = vmul(S3, vmul(vsq(td{dd}), td{jj}));
        v = v + vmul(S2, 2 * vmul(tdj{dd,jj}, td{dd}) + vmul(tdj{dd,dd}, td{jj}));
        v = v + vmul(S1, tddj{dd,jj});
        M(:, dd, jj) = max(abs(v), [], 2);
    end
end
end

function r = vmul(a, b)
P = [a(:,1).*b(:,1), a(:,1).*b(:,2), a(:,2).*b(:,1), a(:,2).*b(:,2)];
r = [min(P, [], 2), max(P, [], 2)];
end
function r = vsq(a)
q = [a(:,1).^2, a(:,2).^2]; r = [min(q, [], 2), max(q, [], 2)];
r(a(:,1) <= 0 & a(:,2) >= 0, 1) = 0;
end
function r = vrng(f, t, crit)
f1 = f(t(:,1)); f2 = f(t(:,2)); lo = min(f1, f2); hi = max(f1, f2);
for c = crit
    in = t(:,1) < c & t(:,2) > c; fc = f(c);
    lo(in) = min(lo(in), fc); hi(in) = max(hi(in), fc);
end
r = [lo, hi];
end

function M = assemble(G, L, R, Q, Q12, S, Sx, d)
% max |psi_ddj| from interval data (one sample)
t = imul(G, L);
sg = @(x) 1 ./ (1 + exp(-x));
f1 = @(x) sg(x) .* (1 - sg(x)); f2 = @(x) f1(x) .* (1 - 2*sg(x));
f3 = @(x) f1(x) .* (1 - 6*sg(x) + 6*sg(x).^2);
S1 = rng_(f1, t, 0); S2 = rng_(f2, t, log(2 + sqrt(3)) * [-1, 1]); S3 = rng_(f3, t, [0, -2.2924, 2.2924]);
m = numel(R);
% first derivatives of t
td = cell(1, d); td{1} = L; td{2} = G;
for i = 1:m, td{2+i} = imul(G, R{i}); end
% second derivatives of t (symmetric table)
Z = [0, 0]; tdj = repmat({Z}, d, d);
tdj{1,2} = [1, 1]; tdj{2,1} = [1, 1];
for i = 1:m
    tdj{1,2+i} = R{i}; tdj{2+i,1} = R{i};
    tdj{2+i,2+i} = imul(G, Q{i});
end
if m == 2, tdj{3,4} = imul(G, Q12); tdj{4,3} = tdj{3,4}; end
% third derivatives t_ddj
tddj = repmat({Z}, d, d);
for i = 1:m
    tddj{2+i,1} = Q{i};
    tddj{2+i,2+i} = imul(G, S{i});
    if m == 2, tddj{2+i,2+(3-i)} = imul(G, Sx{i}); end
end
M = zeros(1, d, d);
for dd = 1:d
    for jj = 1:d
        v = imul(S3, imul(isq(td{dd}), td{jj}));
        v = iadd(v, imul(S2, iadd(imul([2, 2], imul(tdj{dd,jj}, td{dd})), imul(tdj{dd,dd}, td{jj}))));
        v = iadd(v, imul(S1, tddj{dd,jj}));
        M(1, dd, jj) = max(abs(v));
    end
end
end

function [R, Q, S, Vhi, E] = rqs(Elo, Ehi, Vc, sr, fr, dfr, nmax, c6)
% derivatives of log x along one pole: a = -1 + s*i, a' = -phi*i, a'' = |phi'|*i
Vhi = Vc + (Ehi - Elo)^2 / 4; E = [Elo, Ehi];
R = [-1 + sr(1) * Elo, -1 + sr(2) * Ehi];
Q = [-fr(2) * Ehi, -fr(1) * Elo + sr(2)^2 * Vhi];
K3 = min(nmax * Vhi, nmax^3 * c6);
S = [dfr(1) * Elo - 3 * sr(2) * fr(2) * Vhi - sr(2)^3 * K3, dfr(2) * Ehi + sr(2)^3 * K3];
end

function r = rng_(f, t, crit)
c = crit(crit > t(1) & crit < t(2)); v = f([t(1), t(2), c]); r = [min(v), max(v)];
end
function r = imul(a, b)
v = [a(1)*b(1), a(1)*b(2), a(2)*b(1), a(2)*b(2)]; r = [min(v), max(v)];
end
function r = iadd(a, b)
r = [a(1) + b(1), a(2) + b(2)];
end
function r = isq(a)
v = [a(1)^2, a(2)^2];
if a(1) <= 0 && a(2) >= 0, r = [0, max(v)]; else, r = [min(v), max(v)]; end
end
function [flo, fhi] = frange(i, wa, wb)
f = @(w) exp(-w) .* (1 - exp(-w)).^i;
fa = f(wa); fb = f(wb); flo = min(fa, fb); fhi = max(fa, fb);
ws = log(i + 1); in = ws > wa & ws < wb;
if any(in), fs = exp(-ws) .* (1 - exp(-ws)).^i; fhi(in) = max(fhi(in), fs(in)); end
end
function [m, v] = tgeom(rho, L)
i = (0:L)'; lw = i * log(rho); lw = lw - max(lw); w = exp(lw);
S0 = cumsum(w); S1 = cumsum(i .* w); S2 = cumsum(i.^2 .* w);
m = S1 ./ S0; v = max(S2 ./ S0 - m.^2, 0);
end
function r = closest1(a, b)
if a <= 1 && b >= 1, r = 1; elseif b < 1, r = b; else, r = a; end
end
