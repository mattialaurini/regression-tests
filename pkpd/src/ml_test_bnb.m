function output = ml_test_bnb(varargin)
%ML_TEST_BNB Certified Gaussian-ML identification for one simulated patient.
%   Options: id, noise, seed, and any option of IDENTIFY_ML, including
%   'model' ('onecomp' default: one compartment + effect site; 'doublepole').
%   'progress' (default 5): status line every 5 s during the search; 0 = silent.
d = struct('id', 1, 'noise', 0, 'seed', [], 'omega', 0, 'estimate_e0', true, 'e0_range', [], 'Reltol', 0.01, 'Abstol', 1e-7, ...
    'i_gamma', [1 7], 'i_ell', [-3 3], 'i_w', [3 7.6], 'i_emax', [50 160], ...
    'Ngrid', [], 'q0', [], 'Nhat', [], 'screening', [], 'Nnoscreen', [], ...
    'localsearch', true, 'timeout', inf, 'verbose', false, 'progress', 5, ...
    'model', 'onecomp', 'ls_starts', 20, 'obbt', false, 'obbt_trigger', 0.9, 'obbt_rounds', 2, 'obbt_coords', 'all', ...
    'branching', [], 'rel_branch', 4, 'branch_alpha', 0.25, 'pc_score', 'product', 'split', 'mid', 'split_beta', 0.25, 'errmodel', 'interval', 'refine', [], 'refine_trigger', 0.9);
o = parse_opts(d, varargin{:});
sim = simulate_patient(o.id, o.noise, 'seed', o.seed, 'omega', o.omega);
ge = struct('y', sim.b, 'u', sim.u, 'E0', sim.E0);
pass = rmfield(o, {'id', 'noise', 'seed', 'omega'});
args = [fieldnames(pass)'; struct2cell(pass)'];
[gamma, emax, z, stat] = identify_ml(ge, args{:});
output.total_time = stat.time; output.fobj = stat.fobj; output.rss = stat.rss;
output.sigma_hat = stat.sigma_hat; output.xstar = [gamma, emax]; output.z = z;
output.pole = stat.pole; output.rates = stat.rates; output.iterations = stat.iterations; output.lb = stat.lb; output.counts = stat.counts;
output.converged = stat.converged;
output.gamma_ob = sim.gamma_ob; output.Emax_ob = sim.Emax_ob;
output.err_gamma = abs(gamma - sim.gamma_ob) / sim.gamma_ob;
output.err_Emax = abs(emax - sim.Emax_ob) / sim.Emax_ob;
output.err_norm = norm([output.err_gamma, output.err_Emax]);
output.e0 = stat.e0; output.err_E0 = abs(stat.e0 - sim.E0);     % BIS units
psi = ml_psi(z, sim.u, o.model);
chat = zeros(size(psi)); ok = psi > 0 & psi < 1;
chat(ok) = (psi(ok) ./ (1 - psi(ok))).^(1 / gamma);            % c = (psi/(1-psi))^(1/gamma)
ctrue = sim.Ce / sim.Ce50;
output.nrmse_c = norm(chat - ctrue) / norm(ctrue);
output.rmse_b = sqrt(mean((stat.e0 - emax * psi - sim.b_noiseless).^2));
output.seed = sim.seed; output.b = sim.b; output.u = sim.u; output.E0 = sim.E0;
end
