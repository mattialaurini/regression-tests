function [LB, info] = ml_lower_bound(zmin, zmax, ge, Erange, N, M, opts)
%ML_LOWER_BOUND m2(R,N,M) (Eq. eqn_bound_m2) for the Gaussian ML criterion;
%   m1(R,N) when M == N. Separable affine form (Remark after Prop. 4):
%     f(y, z) = (E0 - b) - Emax * psi(z),
%   with y = Emax in [E-, E+] (E0 known, ge.E0), or y = (E0, Emax) in the rectangle
%   ge.E0range x [E-, E+] (E0 estimated; its four vertices). E0 does not depend on z,
%   so the interpolation error involves only Emax * psi: the error coefficients and
%   the third-order refinement scale with the Emax component of each vertex.
%   z = (gamma, ell, w) [double pole] or (gamma, ell, w1, w2) [one compartment +
%   effect site] on the uniform N^d grid of [zmin, zmax].
%   Error coefficients for vertex v_a (> 0), clipped as in Proposition 1:
%     c-(k,a) = sum_d h_d^2/8 min{0,  v_a * min_B d2psi/dz_d^2}
%     c+(k,a) = sum_d h_d^2/8 max{0,  v_a * max_B d2psi/dz_d^2}
%   The relaxation is a convex QP (2-norm). opts: lb_safety.
zmin = zmin(:); zmax = zmax(:);
yb = ge.y(:); n = numel(yb);
V = Erange(:)';                                   % Emax component of the vertices
if isfield(ge, 'E0range') && ~isempty(ge.E0range)
    e0 = ge.E0range(:)';
    Ve0 = kron(e0, ones(1, numel(V)));            % (E0, Emax) vertices of the rectangle
    V = repmat(V, 1, numel(e0));
else
    Ve0 = ge.E0 * ones(1, numel(V));
end
p = numel(V);
d = numel(zmin);
gr = cell(d, 1);
for i = 1:d, gr{i} = linspace(zmin(i), zmax(i), N); end
h = (zmax - zmin) / (N - 1);

% grid points in column-major order (first coordinate fastest), so the grid
% point with indices (i1,...,id) is column block sub2ind(N*ones(1,d), i1,...,id)
mdl = ''; if isfield(ge, 'mlmodel'), mdl = ge.mlmodel; end
G = zeros(n, p*N^d);
idx = cell(1, d);
[idx{:}] = ndgrid(1:N);
for b = 1:N^d
    z = zeros(d, 1);
    for dd = 1:d, z(dd) = gr{dd}(idx{dd}(b)); end
    psi = ml_psi(z, ge.u, mdl);
    G(:, (b-1)*p + (1:p)) = (ones(n, 1) * Ve0 - yb) - psi * V;
end

if strcmp(mdl, 'firstorder')
    [lo, hi] = ml_d2_bounds_fo(ge.u, zmin, zmax);
elseif d == 3
    [lo, hi] = ml_d2_bounds(ge.u, zmin, zmax);
else
    [lo, hi] = ml_d2_bounds_oc(ge.u, zmin, zmax);
end
Clo1 = zeros(n, p); Chi1 = zeros(n, p);
for a = 1:p
    Clo1(:, a) = min(0, V(a) * lo) * (h.^2 / 8);
    Chi1(:, a) = max(0, V(a) * hi) * (h.^2 / 8);
end
sc = ((N - 1) / (M - 1))^2;
Clo = repmat(sc * Clo1, 1, N^d);
Chi = repmat(sc * Chi1, 1, N^d);

if isfield(opts, 'no_solve') && opts.no_solve
    val = NaN; v = []; ok = false;
elseif isfield(opts, 'errmodel') && strcmp(opts.errmodel, 'structured')
    % structured error model: e = sum_d om_d m_d + rem, 0 <= om_d <= T_d y,
    % |rem_k| <= y sum_d T_d hw_dk (m, hw: midpoint / half-width of the enclosure
    % of d2psi/dz_d^2; the factor theta_d = (z_d-a_d)(b_d-z_d)/2 of the
    % interpolation remainder is common to all rows), intersected with the
    % per-row bounds Clo v <= e <= Chi v.
    T = sc * h(:).^2 / 8;
    [val, v, ok] = structured_bound(G, Clo, Chi, (lo + hi) / 2, ((hi - lo) / 2) * T, T, V, N^d);
elseif isfield(opts, 'errmodel') && strcmp(opts.errmodel, 'third')
    % third-order error model, expanded around the centre of each grid cell:
    % for z in cell c, e = y sum_d theta_d (psi_dd(c_c) + delta), theta_d in [0, T_d],
    % |delta_k| <= sum_j (h_j/2) max_R |psi_ddj| (third-derivative bounds over R).
    % Relaxation: e = sum_{c,d} om_cd M_cd + rem, 0 <= om_cd <= T_d * (weight on
    % the corners of c), sum_c om_cd <= T_d y, |rem| <= rho y; plus Clo v <= e <= Chi v.
    T = sc * h(:).^2 / 8;
    D3 = ml_d3_bounds(ge.u, zmin, zmax);
    rho = zeros(n, 1);
    for dd = 1:d
        rho = rho + T(dd) * (reshape(D3(:, dd, :), n, d) * (h / 2));
    end
    cidx = cell(1, d); [cidx{:}] = ndgrid(1:N-1); nc = (N-1)^d;
    Mc = zeros(n, nc * d); Bc = zeros(nc * d, size(G, 2)); Ag = zeros(d, nc * d);
    Vc = repmat(V(:), N^d, 1)';
    for c = 1:nc
        lowc = zeros(d, 1); for dd = 1:d, lowc(dd) = cidx{dd}(c); end
        zc = zeros(d, 1); for dd = 1:d, zc(dd) = (gr{dd}(lowc(dd)) + gr{dd}(lowc(dd) + 1)) / 2; end
        Mpt = ml_point_d2(zc, ge.u);
        incell = false(1, N^d);
        for b = 1:N^d
            ok_ = true;
            for dd = 1:d, ok_ = ok_ && any(idx{dd}(b) == lowc(dd) + [0, 1]); end
            incell(b) = ok_;
        end
        colmask = kron(incell, true(1, p));
        for dd = 1:d
            k = (c - 1) * d + dd;
            Mc(:, k) = Mpt(:, dd);
            Bc(k, :) = T(dd) * (Vc .* colmask);
            Ag(dd, k) = 1;
        end
        rho = rho + 1e-9 * abs(Mpt) * T;            % rounding margin on the point values
    end
    if isfield(opts, 'rho_scale'), rho = rho * opts.rho_scale; end   % EXPERIMENT ONLY (invalid if < 1)
    [val, v, ok] = structured_bound2(G, Clo, Chi, Mc, Bc, Ag, T, rho, V, N^d);
else
    [~, v, ok, val] = norm_min_lp(G, Clo, Chi, '2');
end
if ok
    LB = max(0, val - opts.lb_safety * (1 + abs(val)));
else
    LB = -inf;
end
info.ncols = size(G, 2); info.ok = ok;
if nargout > 1 && isfield(opts, 'return_matrices') && opts.return_matrices
    info.G = G; info.Clo = Clo; info.Chi = Chi; info.v = v; info.grid = gr;
end
end

function [cert, v, ok] = structured_bound(G, Clo, Chi, Mm, rho, T, V, nb)
% Certified lower bound for the structured error model
%   e = sum_d om_d Mm(:,d) + rem, 0 <= om_d <= T_d y, |rem| <= rho y,
% intersected with Clo v <= e <= Chi v.
[m, r] = size(G); d = numel(T);
Vc = repmat(V(:), nb, 1)';
nv = r + d + 2*m; iw = r + d + m + (1:m);
A = [sparse([-G, -Mm, -speye(m), speye(m)]); ...
     sparse([Clo, -Mm, -speye(m), sparse(m, m)]); ...
     sparse([-Chi, Mm, speye(m), sparse(m, m)]); ...
     sparse([-T(:) * Vc, eye(d), sparse(d, 2*m)]); ...
     sparse([-rho * Vc, sparse(m, d), speye(m), sparse(m, m)]); ...
     sparse([-rho * Vc, sparse(m, d), -speye(m), sparse(m, m)]); ...
     sparse([ones(1, r), zeros(1, d + 2*m)])];
model.A = A;
model.rhs = [zeros(3*m + d + 2*m, 1); 1];
model.sense = [repmat('=', m, 1); repmat('<', 2*m + d + 2*m, 1); '='];
model.obj = zeros(nv, 1);
model.Q = blkdiag(sparse(r + d + m, r + d + m), speye(m));
model.lb = [zeros(r + d, 1); -inf(2*m, 1)]; model.ub = inf(nv, 1);
model.modelsense = 'min';
params.OutputFlag = 0; params.FeasibilityTol = 1e-8; params.OptimalityTol = 1e-8;
params.BarConvTol = 1e-9; params.TimeLimit = 60;
if gurobi_threads() > 0, params.Threads = gurobi_threads(); end
params.Presolve = gurobi_presolve();
cert = -inf; v = []; ok = false;
try
    res = solve_model(model, params);
catch
    return
end
if ~isfield(res, 'x') || isempty(res.x), return; end
v = max(res.x(1:r), 0); v = v / max(sum(v), eps); ok = true;
w = res.x(iw);
if norm(w) == 0, cert = 0; return; end
lam = w / norm(w);
% (i) weak norm duality with one direction (each error description separately)
c1 = min(lam' * G + sum(min(lam .* Clo, lam .* Chi), 1));
c2 = min(lam' * G + Vc .* (sum(T(:)' .* min(0, lam' * Mm)) - sum(abs(lam) .* rho)));
cert = max([0, c1, c2]);
% (ii) safe Lagrangian bound from the solver's multipliers (valid for ANY
% multipliers y, >= 0 on the '<' rows; solver error can only weaken it):
%   g(y) = min_{x in X0} ||w||^2 + y'(A x - b),  X0: v in [0,1], om in [0, T E+],
%   |rem| <= rho E+, w free (all implied by the constraints),
%        = -||c_w||^2/4 + sum min(0,c_v) + sum min(0,c_om) T E+ - sum |c_rem| rho E+ - y'b,
%   with c = A'y split by blocks. min ||w|| >= sqrt(max(g, 0)).
if isfield(res, 'pi') && numel(res.pi) == size(A, 1)
    isle = (model.sense == '<'); Ep = max(V);
    for sgn = [-1, 1]
        y = sgn * res.pi(:);
        y(isle) = max(y(isle), 0);
        c = A' * y;
        cv = c(1:r); co = c(r+1:r+d); cr = c(r+d+1:r+d+m); cw = c(r+d+m+1:end);
        g = -sum(cw.^2) / 4 + sum(min(0, cv)) + sum(min(0, co) .* T(:) * Ep) ...
            - sum(abs(cr) .* rho * Ep) - y' * model.rhs;
        cert = max(cert, sqrt(max(g, 0)));
    end
end
end

function [cert, v, ok] = structured_bound2(G, Clo, Chi, Mc, Bc, Ag, T, rho, V, nb)
% Certified lower bound for e = Mc*om + rem with 0 <= om <= Bc v, Ag om <= T y,
% |rem| <= rho y, Clo v <= e <= Chi v (y = Vc v), by the safe Lagrangian bound.
[m, r] = size(G); K = size(Mc, 2); d = numel(T);
Vc = repmat(V(:), nb, 1)';
nv = r + K + 2*m; iw = r + K + m + (1:m);
A = [sparse([-G, -Mc, -speye(m), speye(m)]); ...
     sparse([Clo, -Mc, -speye(m), sparse(m, m)]); ...
     sparse([-Chi, Mc, speye(m), sparse(m, m)]); ...
     sparse([-Bc, eye(K), sparse(K, 2*m)]); ...
     sparse([-T(:) * Vc, Ag, sparse(d, 2*m)]); ...
     sparse([-rho * Vc, sparse(m, K), speye(m), sparse(m, m)]); ...
     sparse([-rho * Vc, sparse(m, K), -speye(m), sparse(m, m)]); ...
     sparse([ones(1, r), zeros(1, K + 2*m)])];
model.A = A;
model.rhs = [zeros(3*m + K + d + 2*m, 1); 1];
model.sense = [repmat('=', m, 1); repmat('<', 2*m + K + d + 2*m, 1); '='];
model.obj = zeros(nv, 1);
model.Q = blkdiag(sparse(r + K + m, r + K + m), speye(m));
model.lb = [zeros(r + K, 1); -inf(2*m, 1)]; model.ub = inf(nv, 1);
model.modelsense = 'min';
params.OutputFlag = 0; params.FeasibilityTol = 1e-8; params.OptimalityTol = 1e-8;
params.BarConvTol = 1e-9; params.TimeLimit = 60;
if gurobi_threads() > 0, params.Threads = gurobi_threads(); end
params.Presolve = gurobi_presolve();
cert = -inf; v = []; ok = false;
try, res = solve_model(model, params); catch, return; end
if ~isfield(res, 'x') || isempty(res.x), return; end
v = max(res.x(1:r), 0); v = v / max(sum(v), eps); ok = true;
cert = 0;
if isfield(res, 'pi') && numel(res.pi) == size(A, 1)
    isle = (model.sense == '<'); Ep = max(V);
    omax = min(max(Bc, [], 2), (Ep * T(:)' * Ag)');   % implied: om_k <= max_j Bc(k,j) (sum v = 1)
    for sgn = [-1, 1]
        y = sgn * res.pi(:); y(isle) = max(y(isle), 0);
        c = A' * y;
        cv = c(1:r); co = c(r+1:r+K); cr = c(r+K+1:r+K+m); cw = c(r+K+m+1:end);
        g = -sum(cw.^2) / 4 + sum(min(0, cv)) + sum(min(0, co) .* omax) ...
            - sum(abs(cr) .* rho * Ep) - y' * model.rhs;
        cert = max(cert, sqrt(max(g, 0)));
    end
end
end
