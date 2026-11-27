% PAPER_TESTS_ML  Certified Gaussian maximum-likelihood identification
% (IDENTIFY_ML) with one pole parameter ('doublepole': two equal poles) and
% two pole parameters ('onecomp': one compartment + effect site, distinct poles),
% with and without screening, on the same data as PAPER_TESTS_LIN (true PK
% perturbed with log-normal omega, seed = patient id).
% Results: results/tests_ml.mat (T, settings, cfg). Override any field of cfg
% by defining it before running the script. Tolerances: cfg.Reltol (0.05),
% cfg.Abstol (1e-3, only relevant for sigma = 0). cfg.resume_redo: noise levels to redo
% when resuming; cfg.resume_unconverged: redo the unconverged runs (explicit fields, e.g.
% cfg.timeout, take precedence over the file). Noise levels (default): sigma in
% {0, 0.01, 0.1, 1}. cfg.resume = 'file.mat' reuses the runs of an earlier results
% file (same models and settings) and executes only the missing ones.
% cfg.estimate_e0 (default true): E0 estimated with Emax (second linear parameter);
% false: the measured baseline is used. Status: a header line per run
% (position, elapsed time, estimated time remaining) and a status line every
% cfg.progress seconds inside each run (default 15; 0 = silent).
if ~exist('cfg', 'var'), cfg = struct(); end
user_progress = isfield(cfg, 'progress'); user_fields = fieldnames(cfg);
cfg = paper_cfg(cfg, struct('ids', 1:13, 'models', {{'doublepole', 'onecomp'}}, ...
    'resume', '', 'estimate_e0', true, 'sigmas', [0, 0.01, 0.1, 1], 'omega', 0.2, 'timeout', 3600, 'Reltol', 0.05, 'Abstol', 1e-3, 'resume_redo', [], 'resume_unconverged', false, 'pc_score', 'product', 'progress', 15, 'parallel', false, 'gurobi_threads', 0, 'gurobi_presolve', -1, 'qp_solver', 'hybrid', ...
    'out', fullfile('results', 'tests_ml.mat')));
if ischar(cfg.parallel)                                % 'auto': parallel if the toolbox is there
    try
        cfg.parallel = license('test', 'Distrib_Computing_Toolbox') && exist('parpool', 'file') > 0;
    catch
        cfg.parallel = false;
    end
end
if cfg.parallel && ~user_progress, cfg.progress = 0; end
% sequential runs (default): Gurobi multi-threaded (cfg.gurobi_threads, 0 = automatic);
% parallel runs: one Gurobi thread per worker (set automatically by GUROBI_THREADS)
if cfg.parallel, cfg.gurobi_threads = 1; else, gurobi_threads(cfg.gurobi_threads); end
gurobi_presolve(cfg.gurobi_presolve);                  % -1 = automatic (default), 0 = presolve off
qp_solver(cfg.qp_solver);                              % 'hybrid' (default) or 'gurobi'
if ~exist('settings_ml', 'var')
    settings_ml = { ...
        'NS N=3',         {'screening', false, 'Nnoscreen', 3}; ...
        'S N=2 M=3 Nh=3', {'screening', true, 'Ngrid', 2, 'q0', 2, 'Nhat', 3}};
end
if ~exist('results', 'dir'), mkdir('results'); end
T = struct();
if ~isempty(cfg.resume)                  % reuse the runs of an earlier results file
    R = load(cfg.resume); T = R.T;
    if ~isempty(cfg.resume_redo) && isfield(T, 'sigma')              % redo these noise levels
        keepr = ~ismember(T.sigma, cfg.resume_redo); f = fieldnames(T);
        for ii = 1:numel(f), T.(f{ii}) = T.(f{ii})(keepr); end
    end
    prev = false; if isfield(R.cfg, 'estimate_e0'), prev = R.cfg.estimate_e0; end
    if ~isequal(logical(prev), logical(cfg.estimate_e0))
        error('paper_tests_ml:resume', ['%s was produced with estimate_e0 = %d: runs with E0 known and ', ...
            'estimated are different methods and cannot be merged. Start a new results file.'], cfg.resume, prev);
    end
    keep = {'models', 'omega', 'timeout', 'Reltol', 'pc_score'};
    for i = 1:numel(keep)        % from the file, unless set explicitly now
        if isfield(R.cfg, keep{i}) && ~ismember(keep{i}, user_fields), cfg.(keep{i}) = R.cfg.(keep{i}); end
    end
    if cfg.resume_unconverged && isfield(T, 'converged')          % redo the unconverged runs
        keepr = T.converged == 1; f = fieldnames(T);
        fprintf('Resume: %d unconverged runs will be redone\n', sum(~keepr));
        for ii = 1:numel(f), T.(f{ii}) = T.(f{ii})(keepr); end
    end
    if isfield(R, 'settings_ml'), settings_ml = R.settings_ml; end
    cfg.sigmas = unique([R.cfg.sigmas(:)', cfg.sigmas(:)']); cfg.resumed_from = cfg.resume;
    fprintf('Resuming from %s: %d runs reused; noise levels %s\n', cfg.resume, paper_nrows(T), mat2str(cfg.sigmas));
end
jobs = zeros(0, 4);
for sg = cfg.sigmas, for m = 1:numel(cfg.models), for s = 1:size(settings_ml, 1), for id = cfg.ids
    if isfield(T, 'order') && any(T.order == m & T.id == id & T.sigma == sg & T.setting == s), continue; end
    jobs(end+1, :) = [id, m, sg, s];
end, end, end, end
batch = struct('t0', tic, 'done', 0, 'n', size(jobs, 1));
runml = @(id, m, sg, s) ml_test_bnb('id', id, 'noise', sg, 'omega', cfg.omega, 'model', cfg.models{m}, 'estimate_e0', cfg.estimate_e0, ...
    'Reltol', cfg.Reltol, 'Abstol', cfg.Abstol, 'pc_score', cfg.pc_score, 'timeout', cfg.timeout, 'progress', cfg.progress, settings_ml{s, 2}{:});
rowml = @(id, m, sg, s, o) struct('id', id, 'order', m, 'sigma', sg, 'omega', cfg.omega, 'setting', s, ...
    'converged', o.converged, 'time', o.total_time, 'iterations', o.iterations, ...
    'bounds', o.counts.m1 + o.counts.m2 + o.counts.branch_evals, 'fobj', o.fobj, 'lb', o.lb, ...
    'err_gamma', o.err_gamma, 'err_Emax', o.err_Emax, 'nrmse_c', o.nrmse_c, 'rmse_b', o.rmse_b, ...
    'e0', o.e0, 'err_E0', o.err_E0);
if cfg.parallel
    rows = cell(size(jobs, 1), 1);
    parfor j = 1:size(jobs, 1)
        o = runml(jobs(j, 1), jobs(j, 2), jobs(j, 3), jobs(j, 4));
        rows{j} = rowml(jobs(j, 1), jobs(j, 2), jobs(j, 3), jobs(j, 4), o);
        fprintf('[ML] done: sigma %g, model %d, setting %d, patient %d (%.0f s)\n', jobs(j, 3), jobs(j, 2), jobs(j, 4), jobs(j, 1), o.total_time);
    end
    for j = 1:numel(rows), T = paper_append(T, rows{j}); end
    save(cfg.out, 'T', 'settings_ml', 'cfg');
else
    for j = 1:size(jobs, 1)
        id = jobs(j, 1); m = jobs(j, 2); sg = jobs(j, 3); s = jobs(j, 4);
        paper_status(batch, sprintf('[ML] sigma %g, %s, %s, patient %d', sg, cfg.models{m}, settings_ml{s, 1}, id));
        o = runml(id, m, sg, s); batch.done = batch.done + 1;
        T = paper_append(T, rowml(id, m, sg, s, o));
        save(cfg.out, 'T', 'settings_ml', 'cfg');
    end
end
fprintf('Saved %s\n', cfg.out);
