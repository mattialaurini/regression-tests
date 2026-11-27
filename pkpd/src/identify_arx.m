function [gamma, emax, theta, stat] = identify_arx(ge, demo, varargin)
% Reproduction note: this objective is a residual norm, unlike normalized RSS
% in the NIST suite. The adaptive linear-domain expansion is part of the fit.
%   model = 'fixed' (default): denominator fixed to the full 4th-order Schnider
%   model, numerator free (FD_BASIS).
%   model = 'linearized': ARX model of order n_a = n_b = 'order' (1..4) linearized
%   (Gauss-Newton) around the output-error reduction of that order (A_S, B_S) of
%   the Schnider model (order 4: the Schnider model itself):
%     c_hat = (B/A_S) u - sum_i dalpha_i (q^-i / A_S) c_S,  c_S = (B_S/A_S) u,
%   linear in (beta, dalpha). The pole shift dalpha is confined to the prior
%   neighbourhood [dlo, dhi] = 2nd-98th percentiles of the reduced denominators
%   of nprior Schnider models perturbed with log-normal omega_prior; beta is
%   represented in an orthonormal basis (theta_beta) with a box of radius
%   rho*||theta_beta0|| around the preliminary fit. C = product of the two boxes
%   (lin_box = 'prior'), or, default lin_box = 'orth', a box of radius rho*||theta0||
%   around the preliminary fit in an orthonormal basis of [Q_beta, Phi] (neighbourhood
%   in response space; stat.dalpha_in_prior reports whether the certified pole shift
%   lies in the prior range).
%IDENTIFY_ARX Certified identification with the BIS-linearized likelihood and
%   the denominator fixed to the Schnider model of the demographics demo:
%     min_{z in B, theta in C} || w(z) .* (c(z) - Q theta) ||_2,
%   c = g_z^{-1}(b), w = slope of the Hill map at c, Q = orthonormal basis of
%   the fixed-denominator model (FD_BASIS). The criterion is affine in theta,
%   and well conditioned, so the paper's separable relaxation (FD_LOWER_BOUND)
%   certifies it with BNB_SCREENING.
%   Linear-parameter domain C: a preliminary multi-start local fit (theta by
%   unconstrained least squares) gives (z0, theta0); C = theta0 + [-R, R]^nb,
%   R = rho * ||theta0||, through its 2^nb vertices (Remark after Prop. 4).
%   The certificate holds for the problem over B x C; stat.theta_interior
%   reports whether the certified theta lies inside C (if not, enlarge rho).
%   A box from the worst corner of B would be enormous (c explodes for
%   gamma -> 1, Emax -> E0 - min b) and makes the relaxation useless.
%   ge: y (BIS, clipped to <= E0), u, E0.
%   Options: model ('fixed'), order (2, linearized model), lin_box ('orth'), rank_tol
%   (1e-9, rank-revealing basis), omega_prior,
%   nprior, nb (4, fixed model), rho (0.5), prelim_starts (20), Reltol (0.05), Abstol, i_gamma, i_emax, Ngrid,
%   q0, Nhat, screening, Nnoscreen, branching, rel_branch, branch_alpha,
%   localsearch, ls_starts, ls_maxevals, timeout, progress, verbose, emax_margin.
d = struct('rank_tol', 1e-9, 'order', 2, 'lin_box', 'orth', 'model', 'fixed', 'omega_prior', 0.4, 'nprior', 50, 'nb', 4, 'rho', 0.1, 'rho_growth', 2.5, 'rho_max', 1.25, 'prelim_starts', 20, 'Reltol', 0.05, 'Abstol', 1e-7, 'i_gamma', [1 7], ...
    'i_emax', [50 160], 'Ngrid', 3, 'q0', 2, 'Nhat', 3, 'screening', false, 'Nnoscreen', 3, ...
    'branching', 'pseudo', 'rel_branch', 4, 'branch_alpha', 0.25, 'pc_score', 'product', 'split', 'mid', 'split_beta', 0.25, 'localsearch', true, ...
    'ls_starts', 5, 'ls_maxevals', 300, 'timeout', inf, 'progress', 0, 'verbose', false, ...
    'emax_margin', 0.01, 'lb_safety', 1e-9);
o = parse_opts(d, varargin{:});
ge.y = ge.y(:); ge.u = ge.u(:);
s = ge.E0 - ge.y;
if any(s < 0), error('identify_arx:data', 'BIS samples must satisfy b <= E0.'); end
Emin = max(o.i_emax(1), (1 + o.emax_margin) * max(s));
z0min = [o.i_gamma(1); Emin]; z0max = [o.i_gamma(2); o.i_emax(2)];
if strcmp(o.model, 'linearized')
    [Q, lo, hi, info_lin] = lin_basis(demo, ge.u, o);
else
    Q = fd_basis(demo, ge.u, o.nb); lo = -inf(o.nb, 1); hi = inf(o.nb, 1); info_lin = [];
end
nv = size(Q, 2);
% preliminary fit (uncertified): least squares, with the pole shift bounded
prof = @(zz) fd_profile(zz, ge, Q, z0min, z0max, lo, hi);
H = halton_points(o.prelim_starts, 2); best = inf; z0 = (z0min + z0max) / 2;
tp = tic; tlast = 0;
for i = 1:o.prelim_starts
    zs = z0min + (z0max - z0min) .* H(i, :)';
    [zl, fl] = local_search(prof, zs, z0min, z0max, o.ls_maxevals);
    if fl < best, best = fl; z0 = zl; end
    if o.progress > 0 && toc(tp) - tlast >= o.progress
        tlast = toc(tp);
        fprintf('  preliminary fit: start %d/%d, best ||r|| = %.6g at gamma = %.4f, Emax = %.3f (%.0fs)\n', ...
            i, o.prelim_starts, best, z0(1), z0(2), toc(tp));
    end
end
if o.progress > 0
    fprintf('  preliminary fit done: ||r|| = %.6g at gamma = %.4f, Emax = %.3f (%.0fs)\n', best, z0(1), z0(2), toc(tp));
end
[~, theta0] = prof(z0);
Mraw = Q;
if strcmp(o.model, 'linearized') && strcmp(o.lin_box, 'orth')
    % neighbourhood in response space: orthonormalise [Q_beta, Phi] and put the box
    % around the preliminary fit in those coordinates (the prior box on dalpha
    % maps to responses far larger than the concentration and makes the relaxation weak)
    % rank-revealing basis: only directions with singular value > rank_tol*s_max
    % (with poles close to 1 the columns of [Q_beta, Phi] can be almost dependent)
    [Uq, Sq] = svd(Mraw, 0); sv = diag(Sq);
    rk = sum(sv > o.rank_tol * sv(1));
    Q = Uq(:, 1:rk); theta0 = Q' * (Mraw * theta0);
    nv = rk; lo = -inf(nv, 1); hi = inf(nv, 1);
    if o.progress > 0 && rk < size(Mraw, 2)
        fprintf('  linearized basis: rank %d of %d (dependent directions removed)\n', rk, size(Mraw, 2));
    end
end
free = isinf(lo);                                     % coordinates with a neighbourhood around the fit
rho = o.rho; tot_iter = 0; tot_time = 0; attempts = 0; zstart = z0;
while true
attempts = attempts + 1;
R = rho * norm(theta0(free));
if all(free)
    % linear parameters enter through the n+1 vertices of a simplex (as in the
    % paper): regular simplex centred at theta0 with circumradius n*R, which
    % contains the ball of radius R around the preliminary fit
    V = theta0 + nv * R * regular_simplex(nv);
else
    % prior box on dalpha (lin_box = 'prior'): box, 2^nv vertices
    L = lo; Hh = hi; L(free) = theta0(free) - R; Hh(free) = theta0(free) + R;
    B01 = (dec2bin(0:2^nv - 1) - '0')';
    V = L + (Hh - L) .* B01;
end
lbopt = struct('lb_safety', o.lb_safety);
if size(V, 2) == size(V, 1) + 1, Vinv = inv([V; ones(1, size(V, 2))]); else, Vinv = []; end
o.z_init = zstart;
ubfun = @(zz) fd_ub_fast(zz, ge, V, Q, Vinv);
o.ub_lower = @(zz) fd_ub_lower(zz, ge, Q);
lbfun = @(zmin, zmax, N, M, yd) fd_lower_bound(zmin, zmax, ge, V, Q, N, M, lbopt);
ob = o; ob.timeout = o.timeout - tot_time;
[z, theta, stat] = bnb_screening(ubfun, lbfun, z0min, z0max, ob);
if isempty(theta), theta = theta0; z = z0; end      % no iteration (e.g. timeout 0)
tot_iter = tot_iter + stat.iterations; tot_time = tot_time + stat.time;
if all(free)
    lamv = [V; ones(1, size(V, 2))] \ [theta; 1];      % barycentric coordinates in the simplex
    interior = all(lamv > 1e-6);
else
    interior = all(abs(theta(free) - theta0(free)) < R * (1 - 1e-6));
end
% adaptive size of C: if the optimum lies on the boundary of C, the minimum over
% C is limited by C itself; enlarge and repeat (starting from the best point)
if interior || rho * o.rho_growth > o.rho_max * (1 + 1e-12) || tot_time >= o.timeout, break; end
if o.progress > 0
    fprintf('  optimum on the boundary of C (rho = %.3g): enlarging to rho = %.3g\n', rho, rho * o.rho_growth);
end
rho = rho * o.rho_growth; zstart = z;
end
gamma = z(1); emax = z(2);
stat.iterations = tot_iter; stat.time = tot_time; stat.rho = rho; stat.rho_attempts = attempts;
stat.theta_interior = interior;
stat.Q = Q; stat.R = R; stat.theta0 = theta0; stat.z0 = z0; stat.emax_box = [Emin, o.i_emax(2)];
stat.dalpha_at_bound = any(abs(theta(~free) - lo(~free)) < 1e-9 | abs(theta(~free) - hi(~free)) < 1e-9);
stat.lin = info_lin; stat.V = V;
if ~isempty(info_lin)
    % is there a coefficient vector in [Q_beta, Phi] reproducing the certified
    % response with the pole shift inside the prior range?
    target = Q * theta;
    [th_raw, res] = bounded_ls(Mraw, -target, [-inf(o.order, 1); info_lin.dlo], [inf(o.order, 1); info_lin.dhi]);
    stat.dalpha = th_raw(o.order+1:end);
    stat.dalpha_in_prior = res <= 1e-6 * max(norm(target), eps);
end
end

function [M, lo, hi, info] = lin_basis(demo, u, o)
% linearized ARX basis of order o.order around the Schnider model, prior box on dalpha
u = u(:); n = numel(u); na = o.order;
[AS, BS] = schnider_ahat(demo, na, u);
cS = filter(BS, AS, u);
Psi = zeros(n, na); Phi = zeros(n, na);
for i = 1:na
    Psi(:, i) = filter(1, AS, [zeros(i, 1); u(1:n-i)]);
    Phi(:, i) = -filter(1, AS, [zeros(i, 1); cS(1:n-i)]);
end
[Qb, ~] = qr(Psi, 0);
M = [Qb, Phi];
% prior neighbourhood of the pole shift (perturbed Schnider models)
row = [1, demo(1), demo(2), demo(3), demo(4), 1, 1, 100, 100];
[k10, k12, k13, k21, k31, k1e, ke0, V1] = propofol_patient_parameters(1, row);
pnom = [k10, k12, k13, k21, k31, ke0, V1];
persistent prior_cache
if isempty(prior_cache), prior_cache = containers.Map(); end
key = sprintf('%g_', demo, na, o.omega_prior, o.nprior, numel(u), sum(u .* (1:numel(u))'));
if isKey(prior_cache, key)
    A = prior_cache(key);
    if o.progress > 0, fprintf('  prior of the pole shift: cached (%d models, order %d)\n', o.nprior, na); end
else
    st = rng; rng(7, 'twister'); A = zeros(o.nprior, na);
    tp = tic; tlast = 0;
    if o.progress > 0, fprintf('  prior of the pole shift: %d perturbed Schnider models, order %d\n', o.nprior, na); end
    for i = 1:o.nprior
        Ai = schnider_ahat(demo, na, u, 1, pnom .* exp(o.omega_prior * randn(1, 7)));
        A(i, :) = Ai(2:na+1);
        if o.progress > 0 && toc(tp) - tlast >= o.progress
            tlast = toc(tp); fprintf('  prior: %d/%d models (%.0fs)\n', i, o.nprior, toc(tp));
        end
    end
    rng(st);
    if o.progress > 0, fprintf('  prior done (%.0fs)\n', toc(tp)); end
    prior_cache(key) = A;
end
dlo = prct(A, 2)' - AS(2:na+1)'; dhi = prct(A, 98)' - AS(2:na+1)';
lo = [-inf(na, 1); dlo]; hi = [inf(na, 1); dhi];
info = struct('AS', AS, 'BS', BS, 'dlo', dlo, 'dhi', dhi);
end

function S = regular_simplex(n)
% n+1 unit vectors in R^n (columns) forming a regular simplex centred at 0
E = eye(n + 1) - 1 / (n + 1);
[U, ~] = qr(E(:, 1:n), 0);                            % orthonormal basis of the hyperplane
S = U' * E;
S = S ./ sqrt(sum(S.^2, 1));
end

function v = prct(X, pc)
X = sort(X, 1); n = size(X, 1); r = 1 + (n - 1) * pc / 100; f = floor(r);
v = X(f, :) + (r - f) * (X(min(n, f + 1), :) - X(f, :));
end

function U = fd_ub_lower(z, ge, Q)
% unconstrained least squares over all theta: a lower bound of min over C
F = fd_matrix(ge, z, Q);
U = norm(F(:, 1) - F(:, 2:end) * (pinv(F(:, 2:end)) * F(:, 1)));
end

function [U, th] = fd_profile(z, ge, Q, zlo, zhi, lo, hi)
% min over theta of ||F(z)[1; theta]|| with lo <= theta <= hi (BOUNDED_LS)
if any(z(:) < zlo) || any(z(:) > zhi), U = 1e12; th = []; return; end
F = fd_matrix(ge, z, Q);
[th, U] = bounded_ls(F(:, 2:end), F(:, 1), lo, hi);
end
