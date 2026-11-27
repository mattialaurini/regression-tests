function [gamma, emax, z, stat] = identify_ml(ge, varargin)
% Reproduction note: raw BIS samples are used. E0 and Emax are linear
% parameters; gamma, log gain and transformed pole rates are nonlinear.
% Jointly estimating E0 accounts for uncertainty in a finite/noisy baseline
% record without adding a nonlinear search dimension. The inverse-Hill ARX
% criterion instead treats the measured baseline as known.
%IDENTIFY_ML Certified Gaussian maximum-likelihood identification of the
%   PK/PD Wiener model. Option 'model':
%   'onecomp' (default): one-compartment PK + effect-site filter (ML_SIM_X2),
%       z = (gamma, ell, w1, w2), p_i = 1 - exp(-w_i). With a free gain the
%       BIS response is symmetric in the two poles, so only the unordered pair
%       of rates is identifiable: both w use the interval i_w and regions with
%       w1 > w2 throughout are discarded (their mirror images are searched).
%       stat.rates = [fast, slow] rates in 1/s; which one is plasma (k10) and
%       which is effect site (ke0) cannot be told from BIS alone.
%   'doublepole': p1 = p2 (ML_SIM_X), z = (gamma, ell, w).
%   'firstorder': single pole (ML_SIM_X1), z = (gamma, ell, w).
%   Model of the BIS:
%     b(k) = E0 - Emax * psi(k; z) + noise,  noise ~ N(0, sigma^2) i.i.d.
%     psi = c^gamma/(1+c^gamma),  c(k) = exp(ell) * x(k; w),
%     x = double-pole filtered input, pole p = 1 - exp(-w) (ML_SIM_X).
%   For known or unknown sigma, ML is  min ||(E0 - b) - Emax psi(z)||_2.
%   'estimate_e0' (default true): E0 is estimated with Emax, as a second linear
%   parameter, in 'e0_range' (default max(b) + [-10, 5]); stat.e0 is the estimate.
%   With 'estimate_e0', false, ge.E0 is used (the measured baseline).
%   Branch on z=(gamma,ell,w) or z=(gamma,ell,w1,w2), depending on the model.
%   The default linear block is (E0,Emax); only Emax remains free when
%   estimate_e0=false. No clipping of b or condition on min b is needed.
%   Returns gamma, emax, the model-dependent z and stat (stat.fobj is the norm;
%   RSS = fobj^2, sigma_hat = fobj/sqrt(n)).
d = struct('Reltol', 0.01, 'Abstol', 1e-7, 'i_gamma', [1 7], 'i_ell', [-3 3], ...
    'estimate_e0', true, 'e0_range', [], 'i_w', [3 7.6], 'i_emax', [50 160], 'Ngrid', [], 'q0', [], 'Nhat', [], ...
    'screening', [], 'Nnoscreen', [], 'localsearch', true, 'ls_maxevals', 300, ...
    'lb_safety', 1e-9, 'timeout', inf, 'verbose', false, 'progress', 0, ...
    'model', 'onecomp', 'ls_starts', 20, 'Ts', 1, ...
    'obbt', false, 'obbt_trigger', 0.9, 'obbt_rounds', 2, 'obbt_coords', 'all', ...
    'branching', [], 'rel_branch', 4, 'branch_alpha', 0.25, 'pc_score', 'product', 'split', 'mid', 'split_beta', 0.25, 'errmodel', 'interval', 'refine', [], 'refine_trigger', 0.9);
o = parse_opts(d, varargin{:});
% model-dependent defaults for the bounding strategy (measured on patient 1,
% sigma = 0.5, epsilon_r = 0.05; see README.md):
%   onecomp    : screening (m2 on 2 points, m1 on 3) + pseudo-cost branching
%                + third-order refinement (820 s; without: 870 s; longest-edge: 1549 s)
%   doublepole : no screening, m1 on 3 points + pseudo-cost branching
%                (36 s; longest-edge branching: 63 s)
if strcmp(o.model, 'firstorder')
    dflt = struct('screening', false, 'Ngrid', 3, 'q0', 8, 'Nhat', 5, 'Nnoscreen', 3, 'branching', 'pseudo', 'refine', 'none');
elseif strcmp(o.model, 'onecomp')
    dflt = struct('screening', true, 'Ngrid', 2, 'q0', 2, 'Nhat', 3, 'Nnoscreen', 3, 'branching', 'pseudo', 'refine', 'third');
else
    dflt = struct('screening', false, 'Ngrid', 3, 'q0', 8, 'Nhat', 5, 'Nnoscreen', 3, 'branching', 'pseudo', 'refine', 'none');
end
for f = fieldnames(dflt)'
    if isempty(o.(f{1})), o.(f{1}) = dflt.(f{1}); end
end
ge.y = ge.y(:); ge.u = ge.u(:); ge.mlmodel = o.model;
E = o.i_emax;
if o.estimate_e0
    % E0 is a second linear parameter: rectangle e0_range x i_emax (4 vertices).
    % Default range from the data: the first samples are the baseline plus noise.
    if isempty(o.e0_range), o.e0_range = max(ge.y) + [-10, 5]; end
    ge.E0range = o.e0_range(:)';
    if o.obbt
        error('identify_ml:obbt', 'Bound tightening (obbt) assumes E0 known: use ''estimate_e0'', false.');
    end
else
    ge.E0range = [];
end
lbopt = struct('lb_safety', o.lb_safety, 'errmodel', o.errmodel);
ubfun = @(zz) ml_object_fun(zz, ge, E);
lbfun = @(zmin, zmax, N, M, yd) ml_lower_bound(zmin, zmax, ge, yd, N, M, lbopt);
switch o.model
    case 'onecomp'
        z0min = [o.i_gamma(1); o.i_ell(1); o.i_w(1); o.i_w(1)];
        z0max = [o.i_gamma(2); o.i_ell(2); o.i_w(2); o.i_w(2)];
        o.discard = @(zmin, zmax) zmin(3) > zmax(4);
    case {'doublepole', 'firstorder'}
        z0min = [o.i_gamma(1); o.i_ell(1); o.i_w(1)];
        z0max = [o.i_gamma(2); o.i_ell(2); o.i_w(2)];
        o.discard = [];
    otherwise
        error('identify_ml:model', 'Unknown model "%s".', o.model);
end
o.ydom0 = E;
if strcmp(o.refine, 'third')
    ropt = struct('lb_safety', o.lb_safety, 'errmodel', 'third');
    o.refine = @(zmin, zmax, yd) ml_lower_bound(zmin, zmax, ge, yd, o.Nnoscreen, o.Nnoscreen, ropt);
else
    o.refine = [];
end
if o.obbt
    o.tighten = @(zmin, zmax, yd, tau) ml_obbt(zmin, zmax, ge, yd, o.Nnoscreen, tau, struct('coords', o.obbt_coords));
else
    o.tighten = [];
end
[z, ylin, stat] = bnb_screening(ubfun, lbfun, z0min, z0max, o);
if o.estimate_e0
    stat.e0 = ylin(1); emax = ylin(2); stat.e0_range = ge.E0range;
else
    stat.e0 = ge.E0; emax = ylin;
end
if strcmp(o.model, 'onecomp') && z(3) > z(4), z([3 4]) = z([4 3]); end  % mirror image
gamma = z(1);
stat.rss = stat.fobj^2;
stat.sigma_hat = stat.fobj / sqrt(numel(ge.y));
poles = 1 - exp(-z(3:end));
stat.pole = poles;
stat.rates = sort(-log(poles(:)') / o.Ts, 'descend');   % [fast, slow] in 1/s
stat.options = o;
end
