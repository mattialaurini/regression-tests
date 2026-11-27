function [zmin2, zmax2, Erange2, status] = ml_obbt(zmin, zmax, ge, Erange, N, tau, opts)
%ML_OBBT Optimality-based bound tightening of a region for the ML criterion.
%   For every coordinate z_d (and for y = Emax), computes
%     min / max  sum_j v_j c_j   over the relaxation of ML_LOWER_BOUND(R, N, N)
%     restricted to ||G v + e||_2 <= tau,
%   where c_j is the coordinate of the grid point (or vertex) of column j.
%   Validity: for any (y,z) in R with ||f(y,z)|| <= tau, the decomposition
%   v_ab = lambda_a mu_b(z) is feasible and sum_ab v_ab z_b = z (multilinear
%   weights reproduce linear functions), sum_ab v_ab E_a = y. Hence every point
%   of R with objective <= tau lies in the returned box; the removed parts have
%   objective > tau and can be discarded exactly as a pruned region.
%   status: 'ok', 'infeasible' (no point of R has objective <= tau: prune R),
%   or 'failed' (solver trouble: the input box is returned unchanged).
%   opts.coords: 'all' (default) or 'y' (tighten only the linear parameter).
%   opts.margin (default 1e-6): each new bound is relaxed by margin*(width)
%   to absorb solver tolerances.
if nargin < 7, opts = struct(); end
if ~isfield(opts, 'margin'), opts.margin = 1e-6; end
if ~isfield(opts, 'coords'), opts.coords = 'all'; end   % 'all' or 'y' (Emax only)
zmin = zmin(:); zmax = zmax(:); d = numel(zmin);
zmin2 = zmin; zmax2 = zmax; Erange2 = Erange(:)'; status = 'ok';
[~, info] = ml_lower_bound(zmin, zmax, ge, Erange, N, N, ...
    struct('lb_safety', 0, 'return_matrices', true, 'no_solve', true));
G = info.G; Clo = info.Clo; Chi = info.Chi;
[m, r] = size(G); p = numel(Erange);
% coordinates of the columns: column (b-1)*p + a  <->  grid point b, vertex a
idx = cell(1, d); [idx{:}] = ndgrid(1:N);
C = zeros(d + 1, r);
for b = 1:N^d
    for a = 1:p
        j = (b-1)*p + a;
        for dd = 1:d, C(dd, j) = info.grid{dd}(idx{dd}(b)); end
        C(d+1, j) = Erange(a);
    end
end
Gs = sparse(G); Cl = sparse(Clo); Ch = sparse(Chi);
model.A = [-Gs, -speye(m), speye(m); -Cl, speye(m), sparse(m, m); ...
           Ch, -speye(m), sparse(m, m); sparse(ones(1, r)), sparse(1, 2*m)];
model.rhs = [zeros(3*m, 1); 1];
model.sense = [repmat('=', m, 1); repmat('>', 2*m, 1); '='];
elo = full(min(Clo, [], 2)); ehi = full(max(Chi, [], 2));
model.lb = [zeros(r, 1); elo; -tau * ones(m, 1)];
model.ub = [inf(r, 1); ehi; tau * ones(m, 1)];
model.quadcon(1).Qc = blkdiag(sparse(r + m, r + m), speye(m));
model.quadcon(1).q = zeros(r + 2*m, 1);
model.quadcon(1).rhs = tau^2;
model.modelsense = 'min';
params.OutputFlag = 0; params.FeasibilityTol = 1e-8; params.OptimalityTol = 1e-8;
params.BarQCPConvTol = 1e-9; params.TimeLimit = 60;
if gurobi_threads() > 0, params.Threads = gurobi_threads(); end
params.Presolve = gurobi_presolve();
lo = [zmin; Erange(1)]; hi = [zmax; Erange(end)];
if strcmp(opts.coords, 'y'), ks = d + 1; else, ks = 1:d + 1; end
for k = ks
    for sgn = [1, -1]
        model.obj = [sgn * C(k, :)'; zeros(2*m, 1)];
        try
            res = gurobi(model, params);
        catch
            status = 'failed'; return
        end
        if strcmp(res.status, 'INFEASIBLE') || strcmp(res.status, 'INF_OR_UNBD')
            status = 'infeasible'; return
        end
        if ~strcmp(res.status, 'OPTIMAL'), continue; end   % keep the old bound
        val = sgn * res.objval; mg = opts.margin * (hi(k) - lo(k));
        if sgn == 1, lo(k) = max(lo(k), val - mg); else, hi(k) = min(hi(k), val + mg); end
    end
end
zmin2 = lo(1:d); zmax2 = hi(1:d); Erange2 = [lo(d+1), hi(d+1)];
end
