function [zstar, ystar, stat] = bnb_screening(ubfun, lbfun, z0min, z0max, o)
% Reproduction note: screening controls queue priorities and whether to compute
% m1. Only valid m1 bounds (including inherited bounds) enter the certificate.
%BNB_SCREENING Algorithm 1 (branch-and-bound with screening) on a box of
%   nonlinear parameters; the linear-parameter domain C is fixed and handled
%   inside ubfun/lbfun.
%     [U, y] = ubfun(z)            min over y in C of ||f(y,z)|| (routine r(R))
%     LB = lbfun(zmin, zmax, N, M, ydom) m2(R,N,M); equals the valid bound m1 when
%                                  M == N; ydom = range of the linear parameter
%   o fields: Reltol, Abstol, Ngrid, q0, Nhat, screening, Nnoscreen,
%             localsearch, ls_maxevals, timeout, verbose, progress,
%             ls_starts (optional: extra local searches at the root from
%             ls_starts Halton points of the initial box; 0 = none),
%             tighten (optional handle, bound tightening):
%                 [zmin2, zmax2, ydom2, status] = tighten(zmin, zmax, ydom, tau)
%               must return a box containing every point of the region with
%               objective <= tau ('infeasible' if there is none). It is called
%               when the valid bound m1 of a region (if computed) exceeds
%               obbt_trigger*tau (default 0.9); obbt_rounds (default 2) rounds,
%               after which m1 is recomputed on the tightened region,
%             refine (optional handle, tighter but costlier valid bound):
%               L2 = refine(zmin, zmax, ydom); called when a computed valid bound
%               m1 satisfies refine_trigger*tau <= m1 <= tau (default 0.95);
%               the larger of m1 and L2 is used,
%             z_init (optional): starting point of an extra local search at the
%               root (e.g. a preliminary fit), in addition to ls_starts,
%             ls_fun (optional handle, same signature as ubfun): surrogate used inside
%               the local searches, then evaluated with ubfun. Use only if its
%               minimisers are close to those of ubfun (identify_arx no longer uses it),
%             ub_lower (optional handle): a LOWER bound of ubfun; ubfun is not
%               evaluated at a box centre when ub_lower there is >= the incumbent
%               (it cannot improve it). A surrogate that is not a lower bound must
%               not be used here: it could hide improving points and stall the search,
%             ydom0 (optional): initial range of the linear parameter,
%             branching (optional): 'longest' (default, normalised longest edge),
%               'strong' (evaluate the children of every eligible direction and
%               take the direction maximising the smaller child estimate), or
%               'pseudo' (reliability branching: strong evaluation until a
%               direction has rel_branch (default 4) observations, then a
%               pseudo-cost prediction of the gain, scaled by the squared
%               normalised width). Eligible directions: normalised width >=
%               branch_alpha (default 0.25) times the largest one, so every
%               direction keeps shrinking (box diameters -> 0).
%             pc_score (optional): 'min' (score = smaller child bound) or 'product'
%               (score = product of the bound increases of the two children, with
%               separate pseudo-costs for the lower and upper child; floor 1e-6 |L|).
%               'product' is the default of IDENTIFY_ARX and IDENTIFY_ML: -6% bound
%               evaluations overall, -17% with three nonlinear parameters (README.md).
%             split (optional): 'mid' (default) or 'incumbent' (split at the
%               incumbent's coordinate when it lies in the box, at least split_beta of
%               the width from either end); not beneficial in our tests.
%               Child estimates: m2(child,N,M0) with screening (its queue
%               priority), the valid m1(child) without screening (kept; a child
%               with m1 > tau is pruned at once).
%             discard (optional handle: discard(zmin,zmax) true when every point
%             of the region has a symmetric copy with the same objective
%             elsewhere in the initial box; such regions are dropped).
%   progress > 0 prints a status line every 'progress' seconds (and whenever
%   the incumbent improves, marked '*'); 0 = silent.
%
%   Certified global lower bound (printed as LB, returned in stat.lb):
%   every region carries the largest valid bound m1 computed on it or on an
%   ancestor (0 at the root, since the objective is a norm); pruned regions
%   contribute their m1. LB = min over open regions and pruned regions, so
%   f* >= LB at all times. The screening values m2 are NOT valid bounds and
%   are used only for the queue order.
%   On termination with an empty queue, UB - f* <= Reltol*UB + Abstol.
z0min = z0min(:); z0max = z0max(:); d = numel(z0min); w0 = z0max - z0min;
N = o.Ngrid; M0 = 1 + o.q0 * (N - 1);
if o.screening && (o.q0 < 1 || N < 2 || o.Nhat < N)
    error('bnb_screening:grid', 'Requires q0 >= 1, Ngrid >= 2, Nhat >= Ngrid.');
end
if ~isfield(o, 'progress'), o.progress = 0; end
if ~isfield(o, 'verbose'), o.verbose = false; end
% queue rows: [priority, Mbar, cached m2 or NaN, valid LB, ylo, yhi, zmin', zmax']
q = PriorityQueue(1);
q.insert([-inf, M0, NaN, 0, NaN, NaN, z0min', z0max']);
first = true;
UB = inf; zstar = (z0min + z0max) / 2; ystar = [];
pruned_min = inf;
cnt = struct('iterations', 0, 'm1', 0, 'm2', 0, 'reinserted', 0, 'branched', 0, 'pruned', 0, 'symmetric', 0);
hasdiscard = isfield(o, 'discard') && ~isempty(o.discard);
hastight = isfield(o, 'tighten') && ~isempty(o.tighten);
if ~isfield(o, 'obbt_trigger'), o.obbt_trigger = 0.9; end
if ~isfield(o, 'obbt_rounds'), o.obbt_rounds = 2; end
if isfield(o, 'ydom0') && ~isempty(o.ydom0), yd0 = o.ydom0(:)'; else, yd0 = [NaN, NaN]; end
cnt.tightened = 0; cnt.tight_pruned = 0; cnt.branch_evals = 0; cnt.refined = 0; cnt.refine_pruned = 0;
hasrefine = isfield(o, 'refine') && ~isempty(o.refine);
if ~isfield(o, 'refine_trigger'), o.refine_trigger = 0.95; end
if ~isfield(o, 'branching') || isempty(o.branching), o.branching = 'longest'; end
if ~isfield(o, 'rel_branch'), o.rel_branch = 4; end
if ~isfield(o, 'branch_alpha'), o.branch_alpha = 0.25; end
if ~isfield(o, 'pc_score'), o.pc_score = 'min'; end        % 'min' | 'product'
if ~isfield(o, 'split'), o.split = 'mid'; end              % 'mid' | 'incumbent'
if ~isfield(o, 'split_beta'), o.split_beta = 0.25; end
pc_sum = zeros(d, 1); pc_cnt = zeros(d, 1);   % pseudo-costs: gain / (normalised width)^2
pc_lo = zeros(d, 1); pc_hi = zeros(d, 1);     % per child (lower / upper half), for 'product'
t0 = tic; converged = true; tlast = -inf; lastprinted = -1;
if o.progress > 0
    fprintf(['Columns: open%% = volume of the initial box not yet pruned; LB = certified\n', ...
             'lower bound on f* (stays low until every open region has an m1 bound);\n', ...
             'min m2 = smallest screening value in the queue (estimate, not certified).\n']);
    fprintf('%7s %8s %7s %8s %13s %13s %8s %11s %6s %7s %7s %7s\n', 'iter', 'time[s]', 'queue', ...
        'open%', 'UB', 'LB', 'gap', 'min m2', 'm1', 'm2', 'reins', 'pruned');
end
hybrid_tau(NaN);
while q.size() > 0
    if toc(t0) > o.timeout, converged = false; break; end
    node = q.remove();
    Mbar = node(2); m2c = node(3); vlb = node(4); yd = node(5:6);
    if first, yd = yd0; first = false; end
    zmin = node(7:6+d)'; zmax = node(7+d:6+2*d)';
    if hasdiscard && o.discard(zmin, zmax)
        cnt.symmetric = cnt.symmetric + 1;
        continue
    end
    cnt.iterations = cnt.iterations + 1;
    m1 = [];

    % --- feasible point r(R)
    zc = (zmin + zmax) / 2;
    if isfield(o, 'ub_lower') && ~isempty(o.ub_lower) && cnt.iterations > 1 && o.ub_lower(zc) >= UB
        U = inf; yc = [];            % a LOWER bound of ubfun(zc) is >= UB: no improvement possible
    else
        [U, yc] = ubfun(zc);
    end
    zh = zc;
    if o.localsearch && U < UB
        [zl, Ul, yl] = ls_run(ubfun, o, zc, zmin, zmax);
        if Ul < U, U = Ul; yc = yl; zh = zl; end
    end
    if cnt.iterations == 1 && o.localsearch
        starts = zeros(d, 0);
        if isfield(o, 'z_init') && ~isempty(o.z_init), starts = o.z_init(:); end
        if isfield(o, 'ls_starts') && o.ls_starts > 0   % multi-start at the root
            starts = [starts, zmin + (zmax - zmin) .* halton_points(o.ls_starts, d)'];
        end
        for i = 1:size(starts, 2)
            [zl, Ul, yl] = ls_run(ubfun, o, starts(:, i), zmin, zmax);
            if Ul < U, U = Ul; yc = yl; zh = zl; end
        end
    end
    improved = U < UB;
    if improved
        UB = U; zstar = zh; ystar = yc;
        if o.verbose, fprintf('UB = %.6g at z = %s\n', UB, mat2str(zstar', 6)); end
    end
    tau = (1 - o.Reltol) * UB - o.Abstol;
    hybrid_tau(tau);                 % threshold for the 'hybrid' solver mode

    % --- screening / bounding
    branch = false;
    if o.screening
        if isnan(m2c), m2 = lbfun(zmin, zmax, N, Mbar, yd); cnt.m2 = cnt.m2 + 1; else, m2 = m2c; end
        if m2 <= tau
            branch = true;
        elseif Mbar <= o.Nhat
            m1 = lbfun(zmin, zmax, Mbar, Mbar, yd); cnt.m1 = cnt.m1 + 1;
            [m1, cnt] = do_refine(m1, tau, zmin, zmax, yd, o, hasrefine, cnt);
            vlb = max(vlb, m1);
            if m1 <= tau
                branch = true;
            else
                cnt.pruned = cnt.pruned + 1; pruned_min = min(pruned_min, m1);
            end
        else
            qn = max(1, floor((Mbar - 1) / (2 * (N - 1))));
            q.insert([m2, 1 + qn * (N - 1), NaN, vlb, yd, zmin', zmax']);
            cnt.reinserted = cnt.reinserted + 1;
        end
    else
        if ~isnan(m2c)
            m1 = m2c;                        % valid m1 computed when the node was created
        else
            m1 = lbfun(zmin, zmax, o.Nnoscreen, o.Nnoscreen, yd); cnt.m1 = cnt.m1 + 1;
        end
        [m1, cnt] = do_refine(m1, tau, zmin, zmax, yd, o, hasrefine, cnt);
        vlb = max(vlb, m1);
        if m1 <= tau
            branch = true;
        else
            cnt.pruned = cnt.pruned + 1; pruned_min = min(pruned_min, m1);
        end
    end

    % --- optional bound tightening of a nearly prunable region
    if branch && hastight
        if ~isempty(m1), est = m1; else, est = -inf; end   % trigger on the valid bound only
        if est >= o.obbt_trigger * tau
            cnt.tightened = cnt.tightened + 1;
            st = 'ok';
            for rr = 1:o.obbt_rounds
                [zmin, zmax, yd, st] = o.tighten(zmin, zmax, yd, tau);
                if ~strcmp(st, 'ok'), break; end
            end
            if strcmp(st, 'infeasible')
                branch = false; cnt.tight_pruned = cnt.tight_pruned + 1;
                pruned_min = min(pruned_min, tau);   % every removed point has f > tau
            elseif strcmp(st, 'ok')
                pruned_min = min(pruned_min, tau);   % parts cut off have f > tau
                m1 = lbfun(zmin, zmax, o.Nnoscreen, o.Nnoscreen, yd); cnt.m1 = cnt.m1 + 1;
                vlb = max(vlb, m1);
                if m1 > tau
                    branch = false; cnt.tight_pruned = cnt.tight_pruned + 1;
                    pruned_min = min(pruned_min, m1);
                end
            end
        end
    end

    % --- branching (C stays fixed)
    if branch
        cnt.branched = cnt.branched + 1;
        wn = (zmax - zmin) ./ w0;                     % normalised widths
        if o.screening, Lpar = m2; else, Lpar = m1; end
        evalchild = @(a, b) child_estimate(lbfun, a, b, yd, o, N, M0);
        cand = find(wn >= o.branch_alpha * max(wn));  % eligible directions
        % split point in each direction: midpoint, or the incumbent's coordinate if
        % it lies in the box (at least split_beta of the width from either end)
        tsplit = 0.5 * ones(d, 1);
        if strcmp(o.split, 'incumbent') && ~isempty(zstar) && all(zstar >= zmin) && all(zstar <= zmax)
            tsplit = min(max((zstar - zmin) ./ max(zmax - zmin, realmin), o.split_beta), 1 - o.split_beta);
        end
        kids = []; k = [];
        if strcmp(o.branching, 'longest')
            [~, k] = max(wn);                         % ties: smallest index
        else
            best = -inf;
            for kk = cand(:)'
                if strcmp(o.branching, 'pseudo') && pc_cnt(kk) >= o.rel_branch
                    if strcmp(o.pc_score, 'product')
                        eps_ = 1e-6 * max(abs(Lpar), 1e-12);
                        sc = max(pc_lo(kk) / pc_cnt(kk) * wn(kk)^2, eps_) * max(pc_hi(kk) / pc_cnt(kk) * wn(kk)^2, eps_);
                    else
                        sc = Lpar + pc_sum(kk) / pc_cnt(kk) * wn(kk)^2;   % predicted
                    end
                    Lk = [];
                else
                    [a1, b1, a2, b2] = halves(zmin, zmax, kk, tsplit(kk));
                    Lk = [evalchild(a1, b1), evalchild(a2, b2)];      % evaluated
                    cnt.branch_evals = cnt.branch_evals + 2;
                    if isfinite(Lpar) && all(isfinite(Lk))
                        pc_sum(kk) = pc_sum(kk) + max(min(Lk) - Lpar, 0) / wn(kk)^2;
                        pc_lo(kk) = pc_lo(kk) + max(Lk(1) - Lpar, 0) / wn(kk)^2;
                        pc_hi(kk) = pc_hi(kk) + max(Lk(2) - Lpar, 0) / wn(kk)^2;
                        pc_cnt(kk) = pc_cnt(kk) + 1;
                    end
                    if strcmp(o.pc_score, 'product') && isfinite(Lpar)
                        eps_ = 1e-6 * max(abs(Lpar), 1e-12);
                        sc = max(Lk(1) - Lpar, eps_) * max(Lk(2) - Lpar, eps_);
                    else
                        sc = min(Lk);
                    end
                end
                sc = sc + 1e-12 * wn(kk);             % ties: wider direction
                if sc > best, best = sc; k = kk; kids = Lk; end
            end
            if isempty(k), [~, k] = max(wn); end
        end
        [a1, b1, a2, b2] = halves(zmin, zmax, k, tsplit(k));
        kvalid = ~isempty(kids);
        if ~kvalid && (o.screening || ~strcmp(o.branching, 'longest'))
            kids = [evalchild(a1, b1), evalchild(a2, b2)];
            cnt.branch_evals = cnt.branch_evals + 2; kvalid = true;
            if strcmp(o.branching, 'pseudo') && isfinite(Lpar) && all(isfinite(kids))
                pc_sum(k) = pc_sum(k) + max(min(kids) - Lpar, 0) / wn(k)^2;
                pc_lo(k) = pc_lo(k) + max(kids(1) - Lpar, 0) / wn(k)^2;
                pc_hi(k) = pc_hi(k) + max(kids(2) - Lpar, 0) / wn(k)^2;
                pc_cnt(k) = pc_cnt(k) + 1;
            end
        end
        boxes = {a1, b1; a2, b2};
        for c = 1:2
            a = boxes{c, 1}; b = boxes{c, 2};
            if hasdiscard && o.discard(a, b)
                cnt.symmetric = cnt.symmetric + 1; continue
            end
            if o.screening
                Lc = kids(c);
                q.insert([Lc, M0, Lc, vlb, yd, a', b']);
            elseif kvalid
                Lc = kids(c);                         % valid m1 of the child
                if Lc > tau
                    cnt.pruned = cnt.pruned + 1; pruned_min = min(pruned_min, Lc);
                else
                    q.insert([Lc, o.Nnoscreen, Lc, max(vlb, Lc), yd, a', b']);
                end
            else
                q.insert([m1, o.Nnoscreen, NaN, vlb, yd, a', b']);
            end
        end
    end

    if o.progress > 0 && (improved || toc(t0) - tlast >= o.progress)
        tlast = toc(t0); lastprinted = cnt.iterations;
        print_line(cnt, toc(t0), q, UB, global_lb(q, pruned_min, UB, d), improved, d, w0);
    end
end
LB = global_lb(q, pruned_min, UB, d);
if o.progress > 0
    if lastprinted ~= cnt.iterations, print_line(cnt, toc(t0), q, UB, LB, false, d, w0); end
    if converged
        fprintf('Converged: UB - f* <= %g*UB + %g (certified LB = %.6g)\n', o.Reltol, o.Abstol, LB);
    else
        fprintf('Stopped by timeout: f* in [%.6g, %.6g]\n', LB, UB);
    end
end
stat.fobj = UB; stat.lb = LB; stat.time = toc(t0); stat.converged = converged;
stat.iterations = cnt.iterations; stat.counts = cnt;
end

function [zl, Ul, yl] = ls_run(ubfun, o, z0, zmin, zmax)
% local search from z0, on the surrogate ls_fun if given; the result is always
% evaluated with ubfun (feasible incumbent)
if isfield(o, 'ls_fun') && ~isempty(o.ls_fun)
    zl = local_search(o.ls_fun, z0, zmin, zmax, o.ls_maxevals);
    [Ul, yl] = ubfun(zl);
else
    [zl, Ul, yl] = local_search(ubfun, z0, zmin, zmax, o.ls_maxevals);
end
end

function [m1, cnt] = do_refine(m1, tau, zmin, zmax, yd, o, hasrefine, cnt)
% tighter valid bound for nearly prunable regions (the larger value is kept)
if hasrefine && m1 <= tau && m1 >= o.refine_trigger * tau
    cnt.refined = cnt.refined + 1;
    m1 = max(m1, o.refine(zmin, zmax, yd));
    if m1 > tau, cnt.refine_pruned = cnt.refine_pruned + 1; end
end
end

function [a1, b1, a2, b2] = halves(zmin, zmax, k, t)
if nargin < 4, t = 0.5; end
mid = zmin(k) + t * (zmax(k) - zmin(k));
a1 = zmin; b1 = zmax; b1(k) = mid;
a2 = zmin; a2(k) = mid; b2 = zmax;
end

function L = child_estimate(lbfun, a, b, yd, o, N, M0)
% estimate of the lower bound of a child region: m2(child, N, M0) with
% screening (its queue priority), the valid m1(child) without screening
if o.screening
    L = lbfun(a, b, N, M0, yd);
else
    L = lbfun(a, b, o.Nnoscreen, o.Nnoscreen, yd);
end
end

function LB = global_lb(q, pruned_min, UB, d)
% certified lower bound on f*: min over open regions (inherited valid bound)
% and pruned regions (their m1); capped by UB
LB = min(pruned_min, UB);
n = q.size();
if n > 0
    els = q.elements();
    for i = 1:n
        LB = min(LB, els{i}(4));   % column 4: inherited valid bound
    end
end
end

function print_line(cnt, t, q, UB, LB, improved, d, w0)
if UB > 0 && isfinite(UB), gap = 100 * (UB - LB) / UB; else, gap = 0; end
vol = 0; mp = inf; n = q.size();
if n > 0
    els = q.elements();
    for i = 1:n
        r = els{i};
        vol = vol + prod((r(7+d:6+2*d) - r(7:6+d)) ./ w0');
        mp = min(mp, r(1));
    end
end
mark = ' '; if improved, mark = '*'; end
fprintf('%7d %8.1f %7d %7.3f%% %13.6g %13.6g %7.2f%% %11.6g %6d %7d %7d %7d %s\n', ...
    cnt.iterations, t, n, 100 * vol, UB, LB, gap, mp, cnt.m1, cnt.m2, ...
    cnt.reinserted, cnt.pruned, mark);
end
