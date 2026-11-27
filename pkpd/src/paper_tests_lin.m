% PAPER_TESTS_LIN  Certified identification with the BIS-linearized criterion and
% the ARX model of order 1..4 linearized around the Schnider model of each patient
% (IDENTIFY_ARX, model 'linearized'). The TRUE PK model is perturbed (log-normal
% omega on every PK parameter), so the true ARX parameters differ from the
% Schnider nominal ones. Same data as PAPER_TESTS_ML (seed = patient id).
%
% Phase 1 (tuning of N, M, N^): every setting is run on the tuning patients; per
%   order, the best setting without screening and the best with screening are
%   selected by the PAR2 score (time if converged, 2*timeout otherwise).
%   Default grid: N = M = N^ in cfg.tune_ns (2,3,4,5,6,8,10) without screening and
%   (N, M, N^) in the rows of cfg.tune_s with screening (M up to 11); M must satisfy
%   M = 1 + q(N-1), q integer (grid nesting). If cfg.tune_extend (default true) and,
%   for an order, the best setting lies at the upper edge of the tested range
%   (largest N without screening, largest M with screening), larger values are
%   added and tested (N + 2; q + 2 with N^ = M if N^ was M, N^ + 2 otherwise) until
%   the best is interior or cfg.tune_N_max (20) is reached: the time may keep
%   decreasing with N and M. Symmetrically, at the lower edge N - 1 (N >= 2) and
%   q - 1 (q >= 1) are added.
%   Parallel mode (cfg.parallel): the tuning runs are grouped by (order, patient);
%   each worker runs the settings of its pair in sequence, with capping. Gurobi uses
%   one thread per solve (GUROBI_THREADS), avoiding oversubscription.
%   Racing (sequential mode): after each tuning patient, the settings whose
%   cumulative PAR2 exceeds cfg.race (2) times the best are dropped (cfg.race = Inf
%   disables it); tuning runs use cfg.tune_timeout (600 s).
%   Capping (sequential mode): a tuning run is stopped at cfg.cap_factor (2) times the
%   best time of the same order and patient (at least cfg.cap_min = 30 s): a setting
%   slower than that cannot be the best (cfg.cap_factor = Inf disables it).
%   Reuse (cfg.reuse, default true): converged tuning runs with the selected settings
%   are copied into the main phase instead of being repeated.
%   Defaults chosen to limit the total time: 3 tuning patients, 3 + 3 starting
%   settings around the coarse grids selected in all the batches so far (the first
%   setting is the capping reference; the extension explores larger grids when the
%   best lies at the edge).
%   Execution (default): runs one after another, each Gurobi solve multi-threaded
%   (cfg.gurobi_threads = 0: Gurobi's automatic choice; k > 0: k threads).
%   cfg.parallel = true runs them with parfor instead, with one Gurobi thread each.
% Phase 2: the two selected settings are run on all patients, noise levels and
%   orders.
% Results: results/tests_lin.mat (T1: tuning, T2: main, best, settings, cfg).
% Override any field of cfg by defining it before running the script.
% Tolerances: relative cfg.Reltol (0.05) and absolute cfg.Abstol (1e-3). The absolute one
%   only matters for noise-free data (sigma = 0), where the optimal residual is tiny and a
%   relative gap alone would require extremely fine boxes; 1e-3 is a root-mean-square
%   residual below 1e-4 BIS units over 150 samples.
% cfg.resume_redo (with cfg.resume): noise levels whose runs are redone (e.g. 0).
% cfg.resume_unconverged (with cfg.resume): redo the runs that did not converge; fields
%   set explicitly (e.g. a longer cfg.timeout) take precedence over those of the file.
% Noise levels (default): sigma in {0, 0.01, 0.1, 1}.
% cfg.resume = 'file.mat': reuse the tuning (selected settings) and the runs of an
%   earlier results file; only the missing (patient, order, sigma, setting) runs are
%   executed, e.g. cfg = struct('resume', 'results/tests_lin.mat') adds the noise
%   levels of cfg.sigmas that the file does not contain. The tuning options, omega,
%   tolerances and branching rule are taken from the file.
% cfg.parallel = true runs the independent runs of each phase with PARFOR (needs
% the Parallel Computing Toolbox; set cfg.progress = 0 to avoid interleaved output);
% results are then saved at the end of each phase.
% Status: a header line per run (position in the batch, elapsed time, estimated
% time remaining) and, inside each run, a status line every cfg.progress seconds
% (default 15; 0 = silent).
if ~exist('cfg', 'var'), cfg = struct(); end
user_progress = isfield(cfg, 'progress'); user_fields = fieldnames(cfg);
cfg = paper_cfg(cfg, struct('resume', '', 'ids', 1:13, 'orders', 1:4, 'sigmas', [0, 0.01, 0.1, 1], 'omega', 0.2, ...
    'tune_ids', [1, 6, 9], 'tune_sigma', 0.1, 'timeout', 1800, 'tune_timeout', 600, 'race', 2, 'Reltol', 0.05, 'Abstol', 1e-3, 'resume_redo', [], 'resume_unconverged', false, ...
    'tune_ns', [2 3 4], ...                          % N = M = N^ without screening (first = cap reference)
    'tune_s', [2 3 3; 2 5 5; 3 5 5], ...               % (N, M, N^) with screening
    'tune_extend', true, 'tune_N_max', 20, 'cap_factor', 2, 'cap_min', 30, 'reuse', true, ...
    'pc_score', 'product', 'progress', 15, 'parallel', false, 'gurobi_threads', 0, 'gurobi_presolve', -1, 'qp_solver', 'hybrid', 'out', fullfile('results', 'tests_lin.mat')));
if ischar(cfg.parallel)                                % 'auto': parallel if the toolbox is there
    try
        cfg.parallel = license('test', 'Distrib_Computing_Toolbox') && exist('parpool', 'file') > 0;
    catch
        cfg.parallel = false;
    end
end
if cfg.parallel && ~user_progress, cfg.progress = 0; end   % avoid interleaved output
% sequential runs (default): Gurobi multi-threaded (cfg.gurobi_threads, 0 = automatic);
% parallel runs: one Gurobi thread per worker (set automatically by GUROBI_THREADS)
if cfg.parallel, cfg.gurobi_threads = 1; else, gurobi_threads(cfg.gurobi_threads); end
gurobi_presolve(cfg.gurobi_presolve);                  % -1 = automatic (default), 0 = presolve off
qp_solver(cfg.qp_solver);                              % 'hybrid' (default) or 'gurobi'
% ---- resume: reuse the tuning and the runs of an earlier results file
resumed = false; T2prev = struct();
if ~isempty(cfg.resume)
    R = load(cfg.resume);
    T1 = R.T1; best = R.best; settings = R.settings;
    if isfield(R, 'T2'), T2prev = R.T2; end
    if ~isempty(cfg.resume_redo) && isfield(T2prev, 'sigma')       % redo these noise levels
        keepr = ~ismember(T2prev.sigma, cfg.resume_redo); f = fieldnames(T2prev);
        for ii = 1:numel(f), T2prev.(f{ii}) = T2prev.(f{ii})(keepr); end
        fprintf('Resume: the runs with sigma in %s will be redone\n', mat2str(cfg.resume_redo));
    end
    keep = {'tune_ids', 'tune_sigma', 'tune_timeout', 'race', 'tune_ns', 'tune_s', 'tune_extend', ...
            'tune_N_max', 'cap_factor', 'cap_min', 'reuse', 'pc_score', 'omega', 'Reltol', 'timeout'};
    for i = 1:numel(keep)        % from the file, unless set explicitly now
        if isfield(R.cfg, keep{i}) && ~ismember(keep{i}, user_fields), cfg.(keep{i}) = R.cfg.(keep{i}); end
    end
    if cfg.resume_unconverged && isfield(T2prev, 'converged')     % redo the unconverged runs
        keepr = T2prev.converged == 1; f = fieldnames(T2prev);
        fprintf('Resume: %d unconverged runs will be redone\n', sum(~keepr));
        for ii = 1:numel(f), T2prev.(f{ii}) = T2prev.(f{ii})(keepr); end
    end
    cfg.sigmas = unique([R.cfg.sigmas(:)', cfg.sigmas(:)']);
    if (isfield(R.cfg, 'parallel') && ~isequal(R.cfg.parallel, cfg.parallel)) || ...
            (isfield(R.cfg, 'gurobi_threads') && ~isequal(R.cfg.gurobi_threads, cfg.gurobi_threads)) || ...
            ~isequal(getfield_or(R.cfg, 'gurobi_presolve', -1), cfg.gurobi_presolve) || ...
            ~isequal(getfield_or(R.cfg, 'qp_solver', 'gurobi'), cfg.qp_solver)
        warning('paper_tests_lin:resume', ['The execution mode (parallel / Gurobi threads / presolve / QP solver) differs from %s: ', ...
            'the computation times of old and new runs are not comparable.'], cfg.resume);
    end
    cfg.resumed_from = cfg.resume; resumed = true;
    fprintf('Resuming from %s: tuning and %d runs reused; noise levels %s\n', cfg.resume, ...
        paper_nrows(T2prev), mat2str(cfg.sigmas));
end
if ~exist('settings', 'var')
    settings = cell(0, 2);
    for N = cfg.tune_ns(:)', settings(end+1, :) = paper_setting(N, N, N); end
    for r = 1:size(cfg.tune_s, 1), settings(end+1, :) = paper_setting(cfg.tune_s(r, 1), cfg.tune_s(r, 2), cfg.tune_s(r, 3)); end
end
screened = cellfun(@(c) any(cellfun(@(x) isequal(x, true), c(2:2:end))), settings(:, 2));
if ~exist('results', 'dir'), mkdir('results'); end
common = {'model', 'linearized', 'omega', cfg.omega, 'Reltol', cfg.Reltol, 'Abstol', cfg.Abstol, 'pc_score', cfg.pc_score, ...
          'timeout', cfg.timeout, 'progress', cfg.progress};
run1 = @(id, k, sg, s) fd_test_bnb('id', id, 'noise', sg, 'order', k, common{:}, settings{s, 2}{:});
run1t = @(id, k, sg, s) fd_test_bnb('id', id, 'noise', sg, 'order', k, common{:}, settings{s, 2}{:}, 'timeout', cfg.tune_timeout);
row_of = @(ph, id, k, sg, s, o) struct('phase', ph, 'id', id, 'order', k, 'sigma', sg, 'omega', cfg.omega, ...
    'setting', s, 'converged', o.converged, 'time', o.total_time, 'iterations', o.iterations, ...
    'bounds', o.counts.m1 + o.counts.m2 + o.counts.branch_evals, 'fobj', o.fobj, 'lb', o.lb, ...
    'err_gamma', o.err_gamma, 'err_Emax', o.err_Emax, 'nrmse_c', o.nrmse_c, 'rmse_b', o.rmse_b, ...
    'interior', o.theta_interior, 'dalpha_in_prior', o.dalpha_in_prior);

% ---- phase 1: tuning, in rounds (the grid is extended where the best lies at its edge)
batch = struct('t0', tic, 'done', 0, 'n', NaN);
if ~resumed
T1 = struct(); best = struct('order', zeros(0, 1), 'ns', zeros(0, 1), 's', zeros(0, 1));
todo = cell(1, max(cfg.orders));
for k = cfg.orders, todo{k} = 1:size(settings, 1); end
round_ = 0;
while true
    round_ = round_ + 1;
    jobs = zeros(0, 4);
    for k = cfg.orders, for s = todo{k}(:)', for id = cfg.tune_ids
        jobs(end+1, :) = [id, k, cfg.tune_sigma, s];
    end, end, end
    if isempty(jobs), break; end
    batch.n = batch.done + size(jobs, 1);
    fprintf('\n==== tuning round %d: %d runs\n', round_, size(jobs, 1));
    grp = cellfun(@(c) any(cellfun(@(x) isequal(x, true), c(2:2:end))), settings(:, 2));   % true = screening
    % (re)create the run function: the settings list may have grown
    run1t = @(id, k, sg, s) fd_test_bnb('id', id, 'noise', sg, 'order', k, common{:}, settings{s, 2}{:}, 'timeout', cfg.tune_timeout);
    if cfg.parallel
        % parallel over (order, patient) pairs; each worker runs the settings of its
        % pair in sequence, with capping (it needs only the results of that pair)
        G = unique(jobs(:, 1:2), 'rows');                % [id, order]
        grows = cell(size(G, 1), 1); T1c = T1;
        parfor g = 1:size(G, 1)
            id = G(g, 1); k = G(g, 2); rr = {}; tbest = [inf inf];   % best time per group (NS, S)
            if isfield(T1c, 'order')                     % earlier rounds
                for gg = 0:1
                    bsel = T1c.order == k & T1c.id == id & T1c.converged == 1 & grp(T1c.setting) == gg;
                    if any(bsel), tbest(gg + 1) = min(T1c.time(bsel)); end
                end
            end
            for s = todo{k}(:)'
                gi = grp(s) + 1;
                tl = min(cfg.tune_timeout, max(cfg.cap_min, cfg.cap_factor * tbest(gi)));
                o = fd_test_bnb('id', id, 'noise', cfg.tune_sigma, 'order', k, common{:}, settings{s, 2}{:}, 'timeout', tl);
                rr{end+1} = row_of(1, id, k, cfg.tune_sigma, s, o);
                if o.converged, tbest(gi) = min(tbest(gi), o.total_time); end
            end
            grows{g} = rr;
            fprintf('[tuning] done: order %d, patient %d (%d settings)\n', k, id, numel(rr));
        end
        for g = 1:numel(grows), for j = 1:numel(grows{g}), T1 = paper_append(T1, grows{g}{j}); end, end
        batch.done = batch.done + size(jobs, 1);
        save(cfg.out, 'T1', 'settings', 'cfg');
    else
        for k = cfg.orders
            S = todo{k}(:); if isempty(S), continue; end
            alive = true(numel(S), 1); cum = zeros(numel(S), 1);
            for ip = 1:numel(cfg.tune_ids)
                id = cfg.tune_ids(ip);
                for a = find(alive)'
                    s = S(a);
                    tl = cfg.tune_timeout;                 % capping: best time on this order/patient
                    if isfield(T1, 'order')
                        bsel = T1.order == k & T1.id == id & T1.converged == 1 & grp(T1.setting) == grp(s);
                        if any(bsel), tl = min(tl, max(cfg.cap_min, cfg.cap_factor * min(T1.time(bsel)))); end
                    end
                    paper_status(batch, sprintf('[tuning] order %d, %s, patient %d (limit %.0f s)', k, settings{s, 1}, id, tl));
                    o = fd_test_bnb('id', id, 'noise', cfg.tune_sigma, 'order', k, common{:}, settings{s, 2}{:}, 'timeout', tl);
                    batch.done = batch.done + 1;
                    T1 = paper_append(T1, row_of(1, id, k, cfg.tune_sigma, s, o));
                    save(cfg.out, 'T1', 'settings', 'cfg');
                    if o.converged, cum(a) = cum(a) + o.total_time; else, cum(a) = cum(a) + 2 * cfg.tune_timeout; end
                end
                % racing within each group (with / without screening), also against the
                % settings of the same group completed in earlier rounds
                drop = false(size(alive));
                for g = [false, true]
                    ing = grp(S) == g;
                    if ~any(alive & ing), continue; end
                    ref = min([cum(alive & ing); paper_prev_cum(T1, k, cfg.tune_ids(1:ip), S, ...
                        numel(cfg.tune_ids), cfg.tune_timeout, find(grp == g))]);
                    drop = drop | (alive & ing & cum > cfg.race * ref);
                end
                if any(drop)
                    fprintf('  racing: dropping %s for order %d\n', strjoin(settings(S(drop), 1)', ', '), k);
                    batch.n = batch.n - sum(drop) * (numel(cfg.tune_ids) - ip);
                end
                alive = alive & ~drop;
            end
        end
    end
    % selection
    screened = cellfun(@(c) any(cellfun(@(x) isequal(x, true), c(2:2:end))), settings(:, 2));
    NMN = paper_setting_nmn(settings);
    best = struct('order', zeros(0, 1), 'ns', zeros(0, 1), 's', zeros(0, 1));
    todo = cell(1, max(cfg.orders));
    for k = cfg.orders
        screened = cellfun(@(c) any(cellfun(@(x) isequal(x, true), c(2:2:end))), settings(:, 2));
        NMN = paper_setting_nmn(settings);              % (the list may have grown)
        score = inf(size(settings, 1), 1); tested = false(size(settings, 1), 1);
        for s = 1:size(settings, 1)
            sel = T1.order == k & T1.setting == s;
            tested(s) = any(sel);
            par2 = T1.time(sel); par2(~T1.converged(sel)) = 2 * cfg.tune_timeout;
            if sum(sel) == numel(cfg.tune_ids), score(s) = mean(par2); end   % else dropped by racing
        end
        sn = score; sn(screened) = inf; [~, bns] = min(sn);
        ss = score; ss(~screened) = inf; [~, bs] = min(ss);
        if ~any(screened), bs = bns; end
        if all(screened), bns = bs; end
        best.order(end+1, 1) = k; best.ns(end+1, 1) = bns; best.s(end+1, 1) = bs;
        fprintf('order %d: best without screening %s (PAR2 %.0f s), with screening %s (PAR2 %.0f s)\n', ...
            k, settings{bns, 1}, score(bns), settings{bs, 1}, score(bs));
        if ~cfg.tune_extend, continue; end
        % extension: best at the upper edge of the tested range
        if ~screened(bns) && isfinite(score(bns))
            Nb = NMN(bns, 1); Nt = NMN(tested & ~screened, 1);
            if Nb == max(Nt) && Nb + 2 <= cfg.tune_N_max
                [settings, idx] = paper_add_setting(settings, Nb + 2, Nb + 2, Nb + 2);
                todo{k}(end+1) = idx;
                fprintf('  order %d: best N = %d at the upper edge of the grid: adding N = %d\n', k, Nb, Nb + 2);
            elseif Nb == min(Nt) && Nb > 2                % lower edge: try N - 1
                [settings, idx] = paper_add_setting(settings, Nb - 1, Nb - 1, Nb - 1);
                todo{k}(end+1) = idx;
                fprintf('  order %d: best N = %d at the lower edge of the grid: adding N = %d\n', k, Nb, Nb - 1);
            end
        end
        % screening: for every N, extend M if the best setting with that N lies at the
        % edge of its M range and is not hopeless (within cfg.race of the best)
        if any(screened & isfinite(score))
            for Nv = unique(NMN(screened & tested, 1))'
                cand = find(screened & NMN(:, 1) == Nv & isfinite(score));
                if isempty(cand), continue; end
                [sc, b] = min(score(cand)); b = cand(b); t = NMN(b, :);
                Mt = NMN(screened & tested & NMN(:, 1) == Nv, 2);
                if sc > cfg.race * score(bs), continue; end
                q0 = (t(2) - 1) / (t(1) - 1);
                if t(2) == max(Mt)                            % upper edge: q + 2
                    q = q0 + 2;
                elseif t(2) == min(Mt) && q0 > 1             % lower edge: q - 1
                    q = q0 - 1;
                else
                    continue
                end
                Mn = 1 + q * (t(1) - 1);
                if t(3) == t(2), Nhn = Mn; else, Nhn = max(t(1), min(Mn, t(3) + sign(q - q0) * 2)); end
                if Mn == t(1), continue; end                  % (N, N, N) is the setting without screening
                if Mn <= cfg.tune_N_max
                    [settings, idx] = paper_add_setting(settings, t(1), Mn, Nhn);
                    todo{k}(end+1) = idx;
                    fprintf('  order %d: best M = %d for N = %d at the edge of the grid: adding (N,M,N^) = (%d,%d,%d)\n', k, t(2), Nv, t(1), Mn, Nhn);
                end
            end
        end
        todo{k} = setdiff(todo{k}, find(tested));       % never rerun a tested setting
    end
    save(cfg.out, 'T1', 'best', 'settings', 'cfg');
    if ~cfg.tune_extend, break; end
end
save(cfg.out, 'T1', 'best', 'settings', 'cfg');
end   % ~resumed

% ---- phase 2: all patients, noise levels and orders, best settings
T2 = T2prev;
run1 = @(id, k, sg, s) fd_test_bnb('id', id, 'noise', sg, 'order', k, common{:}, settings{s, 2}{:});
jobs = zeros(0, 4);
nreused = 0;
for sg = cfg.sigmas, for k = cfg.orders
    for which = unique([best.ns(best.order == k), best.s(best.order == k)])
        for id = cfg.ids
            if isfield(T2, 'order') && any(T2.order == k & T2.id == id & T2.sigma == sg & T2.setting == which)
                continue                        % already in the resumed file
            end
            r = [];
            if cfg.reuse && isfield(T1, 'order')
                r = find(T1.order == k & T1.id == id & T1.sigma == sg & T1.setting == which & T1.converged == 1, 1);
            end
            if ~isempty(r)                      % identical run already done in the tuning
                row = struct(); f = fieldnames(T1);
                for ff = 1:numel(f), row.(f{ff}) = T1.(f{ff})(r); end
                row.phase = 2; T2 = paper_append(T2, row); nreused = nreused + 1;
            else
                jobs(end+1, :) = [id, k, sg, which];
            end
        end
    end
end, end
fprintf('\n==== main phase: %d runs (%d reused from the tuning, %d from the resumed file)\n', ...
    size(jobs, 1), nreused, paper_nrows(T2prev));
batch.n = batch.done + size(jobs, 1);
if cfg.parallel
    rows = cell(size(jobs, 1), 1);
    parfor j = 1:size(jobs, 1)
        o = run1(jobs(j, 1), jobs(j, 2), jobs(j, 3), jobs(j, 4));
        rows{j} = row_of(2, jobs(j, 1), jobs(j, 2), jobs(j, 3), jobs(j, 4), o);
        fprintf('[main] done: sigma %g, order %d, setting %d, patient %d (%.0f s)\n', jobs(j, 3), jobs(j, 2), jobs(j, 4), jobs(j, 1), o.total_time);
    end
    for j = 1:numel(rows), T2 = paper_append(T2, rows{j}); end
    save(cfg.out, 'T1', 'T2', 'best', 'settings', 'cfg');
else
    for j = 1:size(jobs, 1)
        id = jobs(j, 1); k = jobs(j, 2); sg = jobs(j, 3); which = jobs(j, 4);
        paper_status(batch, sprintf('[main] sigma %g, order %d, %s, patient %d', sg, k, settings{which, 1}, id));
        o = run1(id, k, sg, which); batch.done = batch.done + 1;
        T2 = paper_append(T2, row_of(2, id, k, sg, which, o));
        save(cfg.out, 'T1', 'T2', 'best', 'settings', 'cfg');
    end
end
fprintf('Saved %s\n', cfg.out);
