function paper_lin_outputs(matfile, outdir)
%PAPER_LIN_OUTPUTS Figures, tables and text numbers of the ARX experiment for the paper.
%   paper_lin_outputs                     reads results/tests_lin.mat, writes into '.'
%   paper_lin_outputs(matfile, outdir)    outdir = folder of the paper (the manuscript)
%
%   Writes (from the results of PAPER_TESTS_LIN):
%     outdir/Figures/lin_errors_order.pdf    relative errors in gamma and Emax vs order
%                                            (box plots; black circles = minimiser on the
%                                            boundary of the final domain C)
%     outdir/Figures/lin_time_screening.pdf  time with vs without screening
%     outdir/Tables/tab_lin_tuning.tex       table of the tuning of (N, M, N^)  (PAR2)
%     outdir/Tables/tab_lin_results.tex      table of the main results
%     outdir/Tables/lin_numbers.tex          \newcommand macros with every number quoted
%                                            in the text (medians, counts, times, ...)
%   The paper includes the tables with \input{Tables/...} and the macros with
%   \input{Tables/lin_numbers} in the preamble, so re-running the experiments and then
%   this command updates the paper without manual edits.
if nargin < 1 || isempty(matfile), matfile = fullfile('results', 'tests_lin.mat'); end
if nargin < 2 || isempty(outdir), outdir = '.'; end
S = load(matfile);
T = S.T2; T1 = S.T1; cfg = S.cfg; settings = S.settings; best = S.best;
if ~isfield(cfg, 'tune_timeout'), cfg.tune_timeout = cfg.timeout; end
figdir = fullfile(outdir, 'Figures'); tabdir = fullfile(outdir, 'Tables');
if ~exist(figdir, 'dir'), mkdir(figdir); end
if ~exist(tabdir, 'dir'), mkdir(tabdir); end
orders = unique(T.order)'; sig = unique(T.sigma)';
% reference noise levels for the text: 0.1 (or the smallest positive level) and the largest
if any(abs(sig - 0.1) < 1e-12), sl = 0.1; else, sl = min(sig(sig > 0)); end
if isempty(sl), sl = min(sig); end
sh = max(sig);
ns_of = @(k) best.ns(best.order == k); s_of = @(k) best.s(best.order == k);
selb = @(k, sg, w) T.order == k & T.sigma == sg & T.setting == w;
q = @(x, p) paper_prctile(x, p);
npat = numel(unique(T.id));

%% ---- (N, M, N^) of every setting
NMN = zeros(size(settings, 1), 3); scr = false(size(settings, 1), 1);
for s = 1:size(settings, 1)
    o = struct('screening', false, 'Nnoscreen', 3, 'Ngrid', 3, 'q0', 2, 'Nhat', 3);
    a = settings{s, 2};
    for i = 1:2:numel(a), o.(a{i}) = a{i+1}; end
    scr(s) = o.screening;
    if o.screening, NMN(s, :) = [o.Ngrid, 1 + o.q0 * (o.Ngrid - 1), o.Nhat];
    else, NMN(s, :) = o.Nnoscreen * [1 1 1]; end
end

%% ---- tuning table
P = nan(size(settings, 1), numel(orders)); full = true(size(P)); ran = false(size(P)); tst = false(size(P));
for s = 1:size(settings, 1)
    for j = 1:numel(orders)
        sel = T1.order == orders(j) & T1.setting == s;
        if ~any(sel), continue; end
        tst(s, j) = true;                                 % tested for this order
        ran(s, j) = sum(sel) == numel(cfg.tune_ids);      % completed (not dropped by racing)
        t = T1.time(sel); t(~T1.converged(sel)) = 2 * cfg.tune_timeout;
        P(s, j) = mean(t); full(s, j) = all(T1.converged(sel));
    end
end
L = {'\begin{table}[tb]', ...
    ['\caption{Tuning of the grid parameters for the approximate maximum-likelihood criterion: ', ...
     'mean PAR2 score [s] on patients ', listing(cfg.tune_ids), ' ($\sigma=', num(cfg.tune_sigma), ...
     '$, time limit ', num(cfg.tune_timeout), '~s). All runs converged except those marked with $^\dagger$; --: configuration not tested for that order or discarded by racing. ', ...
     'Bold: configurations selected for each order, without and with screening}'], ...
    '\label{tab:lin_tuning}', '\centering', ['\begin{tabular}{lc', repmat('r', 1, numel(orders)), '}'], '\toprule', ...
    [' & $(N,M,\widehat N)$', sprintf(' & $n_c=%d$', orders), ' \\'], '\midrule'};
groups = {find(~scr), 'without', 'screening'; find(scr), 'with', 'screening'};
for g = 1:2
    idx = groups{g, 1};
    for r = 1:numel(idx)
        s = idx(r);
        if r == 1, lab = groups{g, 2}; elseif r == 2, lab = groups{g, 3}; else, lab = ''; end
        row = sprintf('%s & $(%d,%d,%d)$', lab, NMN(s, :));
        for j = 1:numel(orders)
            if isnan(P(s, j)) || ~ran(s, j), c = '--';
            else
                c = sprintf('%.0f', P(s, j));
                if ~full(s, j), c = [c '$^\dagger$']; end
                if s == ns_of(orders(j)) || s == s_of(orders(j)), c = ['\textbf{' c '}']; end
            end
            row = [row ' & ' c];
        end
        L{end+1} = [row ' \\'];
    end
    if g == 1, L{end+1} = '\midrule'; end
end
L = [L, {'\bottomrule', '\end{tabular}', '\end{table}'}];
writelines_(fullfile(tabdir, 'tab_lin_tuning.tex'), L);

%% ---- results table
L = {'\begin{table}[tb]', ...
    ['\caption{Certified identification with the approximate maximum-likelihood criterion (', ...
     sprintf('%d', npat), ' patients, PK parameters perturbed with $\omega=', num(cfg.omega), '$): ', ...
     'median and 90th percentile of the computation time (best configuration without screening), ', ...
     'of the relative errors in $\gamma$ and $E_{\max}$, and of the normalized RMS error of the identified ', ...
     'concentration; last column: runs whose certified minimizer lies in the interior of the final $\C$}'], ...
    '\label{tab:lin_results}', '\centering', '\setlength{\tabcolsep}{4pt}', ...
    '\begin{tabular}{cc rr rr rr rr c}', '\toprule', ...
    ['$\sigma$ & $n_c$ & \multicolumn{2}{c}{time [s]} & \multicolumn{2}{c}{$e_\gamma$ [\%]} & ', ...
     '\multicolumn{2}{c}{$e_{E_{\max}}$ [\%]} & \multicolumn{2}{c}{NRMSE $c$ [\%]} & interior \\'], ...
    ' & & med & 90\% & med & 90\% & med & 90\% & med & 90\% & \\', '\midrule'};
for a = 1:numel(sig)
    for j = 1:numel(orders)
        k = orders(j); sel = selb(k, sig(a), ns_of(k));
        if j == 1, lab = num(sig(a)); else, lab = ''; end
        L{end+1} = sprintf('%s & %d & %.0f & %.0f & %.1f & %.1f & %.1f & %.1f & %.1f & %.1f & %d/%d \\\\', ...
            lab, k, q(T.time(sel), 50), q(T.time(sel), 90), 100*q(T.err_gamma(sel), 50), 100*q(T.err_gamma(sel), 90), ...
            100*q(T.err_Emax(sel), 50), 100*q(T.err_Emax(sel), 90), 100*q(T.nrmse_c(sel), 50), 100*q(T.nrmse_c(sel), 90), ...
            sum(T.interior(sel)), sum(sel));
    end
    if a < numel(sig), L{end+1} = '\midrule'; end
end
L = [L, {'\bottomrule', '\end{tabular}', '\end{table}'}];
writelines_(fullfile(tabdir, 'tab_lin_results.tex'), L);

%% ---- numbers quoted in the text
med = @(k, sg, f) 100 * q(T.(f)(selb(k, sg, ns_of(k))), 50);
hi_orders = orders(orders >= 3); M = {};
add = @(name, val) {sprintf('\\newcommand{\\%s}{%s}', name, val)};
M = [M, add('LinNpatients', sprintf('%d', npat)), add('LinOmega', num(cfg.omega)), ...
     add('LinSigmaLow', num(sl)), add('LinSigmaHigh', num(sh)), add('LinReltol', num(cfg.Reltol)), ...
     add('LinTimeout', num(cfg.timeout)), add('LinTuneTimeout', num(cfg.tune_timeout)), ...
     add('LinTuneSigma', num(cfg.tune_sigma)), add('LinTuneIds', listing(cfg.tune_ids)), ...
     add('LinRaceFactor', words(cfg.race))];
M = [M, add('LinTuneNSList', setlist(NMN(~scr, 1)', true)), add('LinTuneSList', tuplelist(NMN(scr, :)))];
M = [M, add('LinTuneNSCount', words(sum(~scr))), add('LinTuneSCount', words(sum(scr)))];
if isfield(cfg, 'tune_N_max'), M = [M, add('LinTuneNmax', sprintf('%d', cfg.tune_N_max))];
else, M = [M, add('LinTuneNmax', sprintf('%d', max(NMN(:))))]; end
nadd = 0;
if isfield(cfg, 'tune_ns') && isfield(cfg, 'tune_s')
    nadd = size(settings, 1) - numel(cfg.tune_ns) - size(cfg.tune_s, 1);
end
if ~isfield(cfg, 'tune_extend') || ~cfg.tune_extend, ext = 'the extension was not used in this run';
elseif nadd == 0, ext = 'no configuration had to be added';
elseif nadd == 1, ext = 'one configuration was added in this way';
else, ext = sprintf('%s configurations were added in this way', words(nadd)); end
M = [M, add('LinTuneExtOutcome', ext)];
% capping and reuse (sentences, empty when not used in the run)
if isfield(cfg, 'cap_factor') && isfinite(cfg.cap_factor)
    if cfg.cap_factor == 2, cf = 'twice'; else, cf = [words(cfg.cap_factor) ' times']; end
    capS = sprintf(['Moreover, each tuning run was stopped after %s the best time obtained ', ...
        'for the same order and patient by a configuration of the same group (at least %s~s): ', ...
        'such a configuration cannot be the best, and the run counts as unsuccessful in the ', ...
        'PAR2 score.'], cf, num(cfg.cap_min));
else, capS = ''; end
nre = 0;
if isfield(cfg, 'reuse') && cfg.reuse
    for r = 1:numel(T.id)
        nre = nre + any(T1.id == T.id(r) & T1.order == T.order(r) & T1.sigma == T.sigma(r) & ...
            T1.setting == T.setting(r) & T1.time == T.time(r));
    end
end
if nre > 0
    reS = sprintf(['The %d converged tuning runs that coincide with runs of the main experiment ', ...
        '(same patient, order, noise level and configuration) were reused.'], nre);
else, reS = ''; end
M = [M, add('LinCapSentence', capS), add('LinReuseSentence', reS)];
% branching score actually used in the run (files without the field: 'min')
sc = 'min'; if isfield(cfg, 'pc_score'), sc = cfg.pc_score; end
if strcmp(sc, 'product')
    br = ['both children are evaluated and the direction maximizing the product ', ...
          'of the bound increases of the two children is selected; afterwards, the two ', ...
          'increases are predicted from the average increases observed for the lower ', ...
          'and upper child, per squared normalized width'];
else
    br = ['both children are evaluated and the direction maximizing the smaller ', ...
          'child bound is selected; afterwards, the bound increase is predicted from ', ...
          'the average observed gain per squared normalized width'];
end
M = [M, add('LinBranchRule', br)];
dropped = sum(tst(:) & ~ran(:));                       % dropped by racing
if dropped == 0, rs = 'none was';
else, rs = sprintf('%s configuration-order pairs were', words(dropped)); end
if isfield(cfg, 'parallel') && isequal(cfg.parallel, true)
    rs = 'this rule applies only to sequential execution and was not used here (parallel execution)';
end
M = [M, add('LinRacingOutcome', rs)];
M = [M, add('LinNruns', sprintf('%d', numel(T.id))), add('LinNconverged', sprintf('%d', sum(T.converged))), ...
     add('LinTimeMainH', sprintf('%.1f', sum(T.time) / 3600)), add('LinTimeTuneH', sprintf('%.1f', sum(T1.time) / 3600))];
M = [M, add('LinGammaTwoLow', sprintf('%.1f', med(2, sl, 'err_gamma'))), add('LinEmaxTwoLow', sprintf('%.1f', med(2, sl, 'err_Emax'))), ...
     add('LinConcTwoLow', sprintf('%.1f', med(2, sl, 'nrmse_c'))), add('LinTimeTwoLow', sprintf('%.0f', q(T.time(selb(2, sl, ns_of(2))), 50))), ...
     add('LinGammaTwoHigh', sprintf('%.1f', med(2, sh, 'err_gamma'))), add('LinEmaxTwoHigh', sprintf('%.1f', med(2, sh, 'err_Emax')))];
if any(orders == 1)
    g1 = min(arrayfun(@(x) med(1, x, 'err_gamma'), sig));
    M = [M, add('LinGammaOneMin', sprintf('%.0f', floor(g1 / 5) * 5))];
else
    M = [M, add('LinGammaOneMin', 'n/a')];                 % order 1 not in this file
end
if ~isempty(hi_orders)
    gh = min(arrayfun(@(k) med(k, sh, 'err_gamma'), hi_orders)); ch = min(arrayfun(@(k) med(k, sh, 'nrmse_c'), hi_orders));
    M = [M, add('LinGammaHighOrdMin', sprintf('%.0f', floor(gh / 10) * 10)), add('LinConcHighOrdMin', sprintf('%.0f', floor(ch / 10) * 10))];
end
bl = selb(2, sl, ns_of(2)) & ~T.interior; bh = selb(2, sh, ns_of(2)) & ~T.interior;
% noise levels, order-2 errors and boundary cases at every level
cs = arrayfun(@(x) num(x), sig, 'UniformOutput', false);
M = [M, add('LinSigmaSet', strjoin(cs, ',\\,')), add('LinSigmaList', join_and(cs)), add('LinNsigma', words(numel(sig)))];
c2 = arrayfun(@(x) sprintf('$%.1f\\%%$ and $%.1f\\%%$ for $\\sigma=%s$', med(2, x, 'err_gamma'), med(2, x, 'err_Emax'), num(x)), sig, 'UniformOutput', false);
M = [M, add('LinOrderTwoBySigma', join_and(c2))];
cb = {};
for x = sig
    nb = sum(selb(2, x, ns_of(2)) & ~T.interior);
    if nb > 0, cb{end+1} = sprintf('%s with $\\sigma=%s$', count_phrase(nb, 'patient'), num(x)); end
end
if isempty(cb), cbs = 'no patient'; else, cbs = join_and(cb); end
M = [M, add('LinBoundTwoBySigma', cbs)];
% runs that did not converge: what they have in common (factual)
nu = find(~T.converged);
if isempty(nu), us = '';
else
    allscr = all(scr(T.setting(nu)));
    cs_u = arrayfun(@(x) num(x), unique(T.sigma(nu))', 'UniformOutput', false);
    if numel(nu) == 1, us = 'The run that reached the time limit'; else, us = sprintf('The %s runs that reached the time limit', words(numel(nu))); end
    if allscr, us = [us ' used screening']; else, us = [us ' include runs with and without screening']; end
    us = sprintf('%s, with $\\sigma\\in\\{%s\\}$.', us, strjoin(cs_u, ',\\,'));
end
M = [M, add('LinUnconvSentence', us)];
if any(sig == 0)
    nf = sprintf(['With noise-free data ($\\sigma=0$), where the errors are due only to the ', ...
        'mismatch between the true PK model and the population model, they are $%.1f\\%%$ and ', ...
        '$%.1f\\%%$, and the identified concentration differs from the true one by $%.1f\\%%$.'], ...
        med(2, 0, 'err_gamma'), med(2, 0, 'err_Emax'), med(2, 0, 'nrmse_c'));
else, nf = ''; end
M = [M, add('LinNoiseFreeSentence', nf)];
M = [M, add('LinBoundTwoLow', count_phrase(sum(bl), 'patient')), add('LinBoundTwoHigh', count_phrase(sum(bh), 'patient'))];
bsel = T.order == 2 & T.setting == ns_of(2) & ~T.interior;          % boundary cases of order 2
if any(bsel), tmax = max(T.time(bsel)); else, tmax = max(T.time(T.order == 2 & T.setting == ns_of(2))); end
M = [M, add('LinTimeMax', sprintf('%.0f', tmax))];
% median time of the other (interior) order-2 runs at the noise levels with boundary cases
ci = {};
for x = sig
    sb = T.order == 2 & T.setting == ns_of(2) & T.sigma == x;
    if any(sb & ~T.interior) && any(sb & T.interior)
        ci{end+1} = sprintf('%.0f~s for $\\sigma=%s$', q(T.time(sb & T.interior), 50), num(x));
    end
end
M = [M, add('LinTimeInteriorBySigma', join_and(ci))];
cand = find(T.order == 2 & ~T.interior & T.sigma == sl);
if isempty(cand), cand = find(T.order == 2 & ~T.interior); end
if ~isempty(cand)
    [~, w] = max(T.err_gamma(cand)); w = cand(w);
    M = [M, add('LinWorstPatient', sprintf('%d', T.id(w))), add('LinWorstGamma', sprintf('%.0f', 100 * T.err_gamma(w))), ...
         add('LinWorstSigma', num(T.sigma(w)))];
end
% screening: time ratios and iteration increase, pairs (patient, order, sigma)
r = []; inc = [];
for k = orders
    for sg = sig
        a = selb(k, sg, ns_of(k)); b = selb(k, sg, s_of(k));
        if ns_of(k) == s_of(k), continue; end
        [ia, ib] = pair(T.id(a), T.id(b)); ta = T.time(a); tb = T.time(b);
        r = [r; tb(ib) ./ ta(ia)];
        ita = T.iterations(a); itb = T.iterations(b);
        inc(end+1) = q(itb(ib), 50) / q(ita(ia), 50) - 1;
    end
end
if ~isempty(r)
    M = [M, add('LinScrRatio', sprintf('%.2f', q(r, 50))), add('LinScrFaster', sprintf('%d', sum(r < 1))), ...
         add('LinScrPairs', sprintf('%d', numel(r))), add('LinItIncMin', sprintf('%+.0f', 100 * min(inc))), ...
         add('LinItIncMax', sprintf('%+.0f', 100 * max(inc)))];
end
% verdict on screening and selected grid sizes (for the text)
if ~isempty(r)
    rm = q(r, 50);
    if rm < 0.9, vd = 'reduces'; elseif rm > 1.1, vd = 'increases'; else, vd = 'has no systematic effect on'; end
    M = [M, add('LinScrVerdict', vd)];
end
nsN = unique(arrayfun(@(k) NMN(ns_of(k), 1), orders));
c = arrayfun(@(x) sprintf('%d', x), nsN, 'UniformOutput', false);
if numel(c) == 1, gs = c{1}; else, gs = [strjoin(c(1:end-1), ', '), ' or ', c{end}]; end
M = [M, add('LinBestNSPoints', gs)];
ab = 1e-7; if isfield(cfg, 'Abstol'), ab = cfg.Abstol; end
M = [M, add('LinAbstol', tol_tex(ab))];
if ab >= 1e-4
    nsamp = 150; if isfield(cfg, 'horizon'), nsamp = cfg.horizon; end
    abs_s = sprintf(['The absolute tolerance only matters for noise-free data, where the optimal ', ...
        'residual is small and a relative gap alone would require extremely fine regions; ', ...
        '$\\epsilon_a=%s$ corresponds to a root-mean-square residual of %s BIS units over the ', ...
        '%d samples.'], tol_tex(ab), tol_tex_rms(ab / sqrt(nsamp)), nsamp);
else, abs_s = ''; end
M = [M, add('LinAbstolSentence', abs_s)];
M = [M, add('LinExecSentence', exec_sentence(cfg))];
M = [M, add('LinSolverSentence', solver_sentence(cfg))];
writelines_(fullfile(tabdir, 'lin_numbers.tex'), [{'% generated by paper_lin_outputs.m: do not edit'}, M]);

%% ---- figure: errors vs order
col = [0.84 0.15 0.16; 0.12 0.47 0.71];
[ev, zero_plot_level] = pkpd_error_plot_values([T.err_gamma; T.err_Emax]);
ev = ev(isfinite(ev) & ev > 0);
ylim_common = [10^floor(log10(min(ev))), 10^ceil(log10(max(ev)))];    % same scale in all panels
if ylim_common(1) == ylim_common(2), ylim_common(2) = 10 * ylim_common(1); end
if numel(sig) > 2, nr = ceil(numel(sig) / 2); nc = 2; else, nr = 1; nc = numel(sig); end
f = figure('Units', 'centimeters', 'Position', [2 2 18 7.5 * nr], 'Color', 'w');
for a = 1:numel(sig)
    ax = subplot(nr, nc, a); hold(ax, 'on');
    for k = orders
        sel = selb(k, sig(a), s_of(k));
        % Display only exact zeros below the smallest positive error, as in
        % the manuscript caption. Never floor positive errors or statistics.
        vals = {pkpd_error_plot_values(T.err_gamma(sel), zero_plot_level), ...
                pkpd_error_plot_values(T.err_Emax(sel), zero_plot_level)};
        inter = T.interior(sel) == 1;
        for v = 1:2
            x0 = k + (v - 1.5) * 0.34; y = vals{v};
            box_(ax, x0, y, col(v, :));
            jit = linspace(-0.07, 0.07, numel(y))';
            h(v) = plot(ax, x0 + jit(inter), y(inter), 'o', 'MarkerSize', 3.5, 'MarkerFaceColor', col(v, :), 'MarkerEdgeColor', col(v, :));
            if any(~inter)
                hb = plot(ax, x0 + jit(~inter), y(~inter), 'ko', 'MarkerSize', 6, 'LineWidth', 0.9);
            end
        end
    end
    set(ax, 'YScale', 'log', 'XTick', orders, 'XLim', [min(orders) - 0.5, max(orders) + 0.5], ...
        'YLim', ylim_common, 'FontSize', 9);
    grid(ax, 'on'); box(ax, 'on'); xlabel(ax, 'ARX order');
    title(ax, sprintf('\\sigma = %s', num(sig(a))));
    if a == 1
        ylabel(ax, 'relative error');
        if exist('hb', 'var'), legend(ax, [h(1), h(2), hb], {'\gamma', 'E_{max}', 'minimizer on \partial C'}, 'Location', 'southwest');
        else, legend(ax, h, {'\gamma', 'E_{max}'}, 'Location', 'southwest'); end
    end
end
save_pdf(f, fullfile(figdir, 'lin_errors_order.pdf'));

%% ---- figure: time with vs without screening
f = figure('Units', 'centimeters', 'Position', [2 2 9.5 9], 'Color', 'w'); ax = axes(f); hold(ax, 'on');
cols = lines(max(orders)); mk = {'s', 'd', 'o', '^', 'v', 'p'}; hh = []; lab = {}; allt = [];
for k = orders
    if ns_of(k) == s_of(k), continue; end
    for a = 1:numel(sig)
        A = selb(k, sig(a), ns_of(k)); B = selb(k, sig(a), s_of(k));
        [ia, ib] = pair(T.id(A), T.id(B)); ta = T.time(A); tb = T.time(B);
        hk = plot(ax, ta(ia), tb(ib), mk{min(a, numel(mk))}, 'MarkerSize', 4, 'MarkerFaceColor', cols(k, :), 'MarkerEdgeColor', cols(k, :));
        allt = [allt; ta(ia); tb(ib)];
        if a == 1, hh(end+1) = hk; lab{end+1} = sprintf('order %d', k); end
    end
end
if ~isempty(allt)
    lim = [10^floor(log10(min(allt))), 10^ceil(log10(max(allt)))];
    plot(ax, lim, lim, 'k--'); set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', lim, 'YLim', lim, 'FontSize', 9);
    grid(ax, 'on'); box(ax, 'on');
    xlabel(ax, 'time without screening [s]'); ylabel(ax, 'time with screening [s]');
    legend(ax, hh, lab, 'Location', 'northwest');
    tt = arrayfun(@(a) sprintf('%s: \\sigma = %s', mk{min(a, numel(mk))}, num(sig(a))), 1:numel(sig), 'UniformOutput', false);
    title(ax, strjoin(tt, ',  '), 'FontWeight', 'normal', 'FontSize', 8);
end
save_pdf(f, fullfile(figdir, 'lin_time_screening.pdf'));
fprintf('Written: %s, %s, %s (figures and tables of the ARX experiment)\n', figdir, ...
    fullfile(tabdir, 'tab_lin_*.tex'), fullfile(tabdir, 'lin_numbers.tex'));
end

%% ---------------------------------------------------------------- helpers
function box_(ax, x0, y, c)
% box plot: quartiles, median, whiskers at min/max
q1 = paper_prctile(y, 25); q2 = paper_prctile(y, 50); q3 = paper_prctile(y, 75); w = 0.13;
patch(ax, x0 + [-w w w -w], [q1 q1 q3 q3], c, 'FaceAlpha', 0.25, 'EdgeColor', c);
plot(ax, x0 + [-w w], [q2 q2], '-', 'Color', c, 'LineWidth', 1.2);
plot(ax, [x0 x0], [min(y) q1], '-', 'Color', c); plot(ax, [x0 x0], [q3 max(y)], '-', 'Color', c);
end

function [ia, ib] = pair(a, b)
[~, ia, ib] = intersect(a(:), b(:));
end

function save_pdf(f, file)
try
    exportgraphics(f, file, 'ContentType', 'vector');    % MATLAB R2020a+
catch
    set(f, 'PaperPositionMode', 'auto');
    print(f, file, '-dpdf');
end
close(f);
end

function s = num(x)
s = sprintf('%g', x);
end

function s = listing(v)
v = v(:)'; c = arrayfun(@(x) sprintf('%d', x), v, 'UniformOutput', false);
if numel(c) == 1, s = c{1}; else, s = [strjoin(c(1:end-1), ', '), ' and ', c{end}]; end
end

function s = setlist(v, ~)
s = ['N=M=\widehat N\in\{', strjoin(arrayfun(@(x) sprintf('%d', x), v, 'UniformOutput', false), ','), '\}'];
end

function s = tuplelist(A)
c = arrayfun(@(i) sprintf('(%d,%d,%d)', A(i, :)), 1:size(A, 1), 'UniformOutput', false);
s = ['(N,M,\widehat N)\in\{', strjoin(c, ','), '\}'];
end

function s = words(n)
w = {'zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten', 'eleven', 'twelve', 'thirteen'};
if n == round(n) && n >= 0 && n <= 13, s = w{n + 1}; else, s = sprintf('%g', n); end
end

function s = count_phrase(n, noun)
if n == 1, s = ['one ' noun]; else, s = [words(n) ' ' noun 's']; end
end

function writelines_(file, L)
fid = fopen(file, 'w');
fprintf(fid, '%s\n', L{:});
fclose(fid);
end

function s = exec_sentence(cfg)
% how the runs were executed (for the text of the paper)
s = exec_sentence0(cfg);
if isfield(cfg, 'gurobi_presolve') && cfg.gurobi_presolve == 0
    s = [s(1:end-1), ' and its presolve disabled.'];
end
end

function s = exec_sentence0(cfg)
if isfield(cfg, 'parallel') && isequal(cfg.parallel, true)
    s = 'The runs were executed in parallel, one per core, with one Gurobi thread each.';
else
    th = 0; if isfield(cfg, 'gurobi_threads'), th = cfg.gurobi_threads; end
    if th == 0
        s = 'The runs were executed one after another, with the automatic multi-threading of Gurobi.';
    elseif th == 1
        s = 'The runs were executed one after another, with one Gurobi thread.';
    else
        s = sprintf('The runs were executed one after another, with %d Gurobi threads.', th);
    end
end
end

function s = join_and(c)
% 'a', 'a and b', 'a, b and c'
if isempty(c), s = ''; elseif numel(c) == 1, s = c{1};
else, s = [strjoin(c(1:end-1), ', '), ' and ', c{end}]; end
end

function s = solver_sentence(cfg)
% how the bound subproblems were solved (for the text of the paper)
mode = 'gurobi'; if isfield(cfg, 'qp_solver'), mode = cfg.qp_solver; end
if strcmp(mode, 'hybrid')
    s = ['Each bound subproblem is first solved approximately by an accelerated ', ...
         'projected-gradient method applied to the problem ', ...
         '\[\min_{v\in\Delta^{r(N)-1}}\;\sum_i\max\{0,\,((G_N+C^-)v)_i,\,-((G_N+C^+)v)_i\}^2,\]', ...
         'which is equivalent to~\eqref{eqn_bnd_m1} with the 2-norm ($C^-$ and $C^+$ collect ', ...
         'the coefficients $c^-_{i,j}$ and $c^+_{i,j}$, and the error vector $e$ is eliminated ', ...
         'row by row); its weak-duality bound and its objective value bracket the optimal value ', ...
         'of the subproblem, and Gurobi is called only when this bracket does not decide the ', ...
         'pruning test. The upper-bound subproblems are solved by the projected-gradient method alone.'];
else
    s = 'All the bound subproblems were solved with Gurobi.';
end
end

function s = tol_tex(x)
% tolerance as TeX: 10^{k} for powers of ten, otherwise %g
k = round(log10(x));
if abs(x - 10^k) <= 1e-12 * x, s = sprintf('10^{%d}', k); else, s = sprintf('%g', x); end
end

function s = tol_tex_rms(x)
% value as TeX, rounded to one significant digit: $a\cdot 10^{k}$
k = floor(log10(x)); a = round(x / 10^k);
if a == 10, a = 1; k = k + 1; end
if a == 1, s = sprintf('$10^{%d}$', k); else, s = sprintf('$%d\\cdot 10^{%d}$', a, k); end
end
