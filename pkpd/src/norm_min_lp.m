function [val, v, ok, cert, status] = norm_min_lp(G, Clo, Chi, norma)
% Reproduction note: hybrid mode first brackets the QP optimum numerically;
% the active pruning threshold determines whether a Gurobi solve is needed.
%NORM_MIN_LP  min ||G v + e||  s.t.  v in the unit simplex,  Clo v <= e <= Chi v.
%   With Clo = Chi = [] (e = 0) this is the upper-bound subproblem
%   min_{y in C} ||f(y,z)||; with the error coefficients it is the
%   relaxation (eqn_bnd_m1)/(eqn_bound_m2). Both bounds use this routine,
%   hence the same linear-parameter domain.  norma: '1', 'inf' or '2'.
%
%   The result does not rely on the solver status (Gurobi may stop with
%   SUBOPTIMAL, especially the barrier method on the QP):
%   val  = exact objective at a FEASIBLE point: v is projected onto the
%          simplex and e is the best error vector for that v (an upper
%          estimate of the minimum; used for UB),
%   cert = RIGOROUS lower bound on the minimum, by weak duality: for any
%          lambda with ||lambda||_dual <= 1 and every feasible (v,e),
%            ||Gv+e|| >= lambda'(Gv+e) >= min_j [lambda'G_j + sum_i min(lambda_i Clo_ij, lambda_i Chi_ij)].
%          lambda is taken from the returned solution (residual direction)
%          and, when available, from the solver duals; the best is kept.
%          At an optimal solution cert equals the minimum.
%   ok   = true when the solver returned a point; status = solver status.
[m, r] = size(G);
Gs = sparse(G);
if isempty(Clo)
    Clo = sparse(m, r); Chi = sparse(m, r);
else
    Clo = sparse(Clo); Chi = sparse(Chi);
end
if strcmp(norma, '2') && strcmp(qp_solver(), 'hybrid')
    % hybrid: native solver for upper-bound QPs (no error terms: any feasible point
    % gives a valid value) and as a first pass for lower bounds; Gurobi only when
    % the native bracket [cert, val] does not decide pruning w.r.t. the threshold
    tau = hybrid_tau();
    if isempty(Clo)
        [val, v, ok, cert] = native_l2(full(G), [], []);
        status = 'NATIVE'; return
    end
    if ~isnan(tau)
        [val, v, ok, cert, decided] = native_l2(full(G), full(Clo), full(Chi), struct('tau', tau));
        if decided, status = 'NATIVE'; return; end
    end
    % undecided: exact solve with Gurobi (generic model below)
end
if strcmp(norma, '2') && strcmp(qp_solver(), 'native')
    % structured solver: v only, on the simplex (NATIVE_L2), then active-set polish
    if isempty(Clo), Clo = []; Chi = []; end
    [val, v, ok, cert] = native_l2(full(G), full(Clo), full(Chi));
    status = 'NATIVE';
    if isempty(Clo), Cl = zeros(m, r); Ch = Cl; else, Cl = full(Clo); Ch = full(Chi); end
    [vp, ok2] = polish2(full(G), Cl, Ch, v);
    if ok2
        ap = full(G) * vp; wp = ap + min(max(-ap, Cl * vp), Ch * vp);
        if norm(wp) <= val, v = vp; val = norm(wp); end
    end
    a = full(G) * v; w = a + min(max(-a, Cl * v), Ch * v);
    if norm(w) > 0
        l = w / norm(w); cert = max(cert, max(0, min(l' * full(G) + sum(min(l .* Cl, l .* Ch), 1))));
    end
    return
end
model = struct();
switch norma
    case '1'
        model.A = [Gs + Clo, -speye(m); -(Gs + Chi), -speye(m); sparse(ones(1, r)), sparse(1, m)];
        model.rhs = [zeros(2*m, 1); 1];
        model.sense = [repmat('<', 2*m, 1); '='];
        model.obj = [zeros(r, 1); ones(m, 1)];
        model.lb = zeros(r + m, 1); model.ub = inf(r + m, 1);
    case 'inf'
        model.A = sparse([Gs + Clo, -ones(m, 1); -(Gs + Chi), -ones(m, 1); ones(1, r), 0]);
        model.rhs = [zeros(2*m, 1); 1];
        model.sense = [repmat('<', 2*m, 1); '='];
        model.obj = [zeros(r, 1); 1];
        model.lb = zeros(r + 1, 1); model.ub = inf(r + 1, 1);
    case '2'
        % variables [v; e; w]: w = G v + e, minimise w'w
        model.A = [-Gs, -speye(m), speye(m); -Clo, speye(m), sparse(m, m); ...
                   Chi, -speye(m), sparse(m, m); sparse(ones(1, r)), sparse(1, 2*m)];
        model.rhs = [zeros(3*m, 1); 1];
        model.sense = [repmat('=', m, 1); repmat('>', 2*m, 1); '='];
        model.obj = zeros(r + 2*m, 1);
        model.Q = blkdiag(sparse(r + m, r + m), speye(m));
        % implied finite bounds on e and w (redundant, they help numerically)
        elo = full(min(Clo, [], 2)); ehi = full(max(Chi, [], 2));
        wlo = full(min(Gs + Clo, [], 2)); whi = full(max(Gs + Chi, [], 2));
        model.lb = [zeros(r, 1); elo; wlo];
        model.ub = [inf(r, 1); ehi; whi];
    otherwise
        error('norm_min_lp:norm', 'Unknown norm "%s".', norma);
end
model.modelsense = 'min';
params.OutputFlag = 0;
if gurobi_threads() > 0, params.Threads = gurobi_threads(); end
params.Presolve = gurobi_presolve();
params.FeasibilityTol = 1e-8;
params.OptimalityTol = 1e-8;
params.BarConvTol = 1e-9;
params.TimeLimit = 60;
val = inf; v = []; ok = false; cert = -inf; status = 'NO_SOLUTION';
try
    result = solve_model(model, params);
    status = result.status;
    ok = isfield(result, 'x') && ~isempty(result.x) && all(isfinite(result.x));
catch err
    warning('norm_min_lp:solver', 'Solver failure: %s', err.message);
end
if ~ok, return; end

% --- feasible point and its exact objective
v = max(full(result.x(1:r)), 0);
if sum(v) <= 0, v = ones(r, 1) / r; else, v = v / sum(v); end
a = full(G * v); lo = full(Clo * v); hi = full(Chi * v);
e = min(max(-a, lo), hi);                 % best e for this v (row-wise)
w = a + e;
val = pnorm(w, norma);

% --- 2-norm: polish the (possibly inexact) solution before certifying.
% The certificate is sensitive to the direction of w because the columns of
% G are much larger than the residual; an exact solve on the support of v
% with the current clipping pattern of e removes the solver error.
if strcmp(norma, '2')
    [vp, ok2] = polish2(full(G), full(Clo), full(Chi), v);
    if ok2
        ap = full(G * vp); ep = min(max(-ap, full(Clo * vp)), full(Chi * vp)); wp = ap + ep;
        if norm(wp) <= val, v = vp; a = ap; e = ep; w = wp; val = norm(wp); end
    end
end

% --- dual certificates
lams = {};
switch norma
    case '2'
        if norm(w) > 0, lams{end+1} = w / norm(w); end
    case '1'
        lams{end+1} = sign(w);
    case 'inf'
        [~, i] = max(abs(w)); l = zeros(m, 1); l(i) = sign(w(i)); lams{end+1} = l;
end
if any(strcmp(norma, {'1', 'inf'})) && isfield(result, 'pi') && numel(result.pi) >= 2*m
    l = full(result.pi(m+1:2*m) - result.pi(1:m));   % sign convention irrelevant:
    if strcmp(norma, '1'), l = min(max(l, -1), 1);  % both signs are tried
    else, l = l / max(1, sum(abs(l))); end
    lams{end+1} = l; lams{end+1} = -l;
end
cert = 0;                                  % norms are nonnegative
Gf = full(G); Cl = full(Clo); Ch = full(Chi);
for t = 1:numel(lams)
    l = lams{t};
    cj = l' * Gf + sum(min(l .* Cl, l .* Ch), 1);  % one value per column j
    cert = max(cert, min(cj));
end
end

function [v, ok] = polish2(G, Cl, Ch, v)
% Active-set polishing of min ||Gv+e||_2, v in simplex, Cl v <= e <= Ch v:
% for a fixed support S of v and a fixed clipping pattern of e, the problem is
% min ||A_S v_S|| s.t. sum v_S = 1, solved exactly from its KKT system;
% negative weights are dropped and the pattern is updated (a few rounds).
ok = false;
r = numel(v);
for round = 1:8
    S = find(v > 1e-4 * max(v));                % columns carrying real weight
    a = G * v; lo = Cl * v; hi = Ch * v;
    rlo = -a < lo; rhi = -a > hi;              % rows with e at a bound
    A = [G(rlo, S) + Cl(rlo, S); G(rhi, S) + Ch(rhi, S)];
    for inner = 1:numel(S)
        H = A' * A; k = numel(S);
        K = [2*H, ones(k, 1); ones(1, k), 0];
        sol = pinv(K) * [zeros(k, 1); 1];
        vs = sol(1:k);
        if all(vs >= -1e-14), break; end
        [~, i] = min(vs); S(i) = []; A(:, i) = [];
        if isempty(S), return; end
    end
    vn = zeros(r, 1); vn(S) = max(vs, 0); 
    if sum(vn) <= 0, return; end
    vn = vn / sum(vn);
    an = G * vn;
    same = isequal(-an < Cl * vn, rlo) && isequal(-an > Ch * vn, rhi);
    v = vn; ok = true;
    if same, break; end
end
end

function n = pnorm(w, norma)
switch norma
    case '1', n = sum(abs(w));
    case 'inf', n = max(abs(w));
    otherwise, n = norm(w);
end
end
