function out = fd_test_bnb(varargin)
%FD_TEST_BNB Simulated patient (optionally with PK model error) + certified
%   identification with IDENTIFY_ARX. Options: id, noise, omega (log-normal
%   deviation of every PK parameter of the true model from Schnider), seed,
%   and any option of IDENTIFY_ARX. Compares with the truth: err_gamma,
%   err_Emax, nrmse_c (identified concentration Q theta vs Ce/Ce50), rmse_b.
names = varargin(1:2:end); own = ismember(names, {'id', 'noise', 'omega', 'seed'});
k = find(own); idx = sort([2*k-1, 2*k]);
o = parse_opts(struct('id', 1, 'noise', 0.1, 'omega', 0, 'seed', []), varargin{idx});
rest = varargin; rest(idx) = [];
patient_profiles
if isempty(o.seed), o.seed = o.id; end
sim = simulate_patient(o.id, o.noise, 'seed', o.seed, 'omega', o.omega);
u = sim.u; E0 = sim.E0; g = sim.gamma_ob; Em = sim.Emax_ob;
ctrue = sim.Ce / sim.Ce50; b0 = sim.b_noiseless; b = min(sim.b, E0);
ge = struct('y', b, 'u', u, 'E0', E0);
[gh, eh, th, st] = identify_arx(ge, Patients(o.id, 2:5), rest{:});
chat = max(st.Q * th, 0);
bhat = E0 - eh * chat.^gh ./ (1 + chat.^gh);
out.xstar = [gh, eh]; out.theta = th; out.fobj = st.fobj; out.lb = st.lb; out.converged = st.converged;
out.total_time = st.time; out.iterations = st.iterations; out.counts = st.counts;
out.err_gamma = abs(gh - g) / g; out.err_Emax = abs(eh - Em) / Em;
out.nrmse_c = norm(chat - ctrue) / norm(ctrue); out.rmse_b = sqrt(mean((bhat - b0).^2));
out.gamma_ob = g; out.Emax_ob = Em; out.theta_interior = st.theta_interior; out.rho = st.rho; out.rho_attempts = st.rho_attempts; out.dalpha_at_bound = st.dalpha_at_bound;
if isfield(st, 'dalpha_in_prior'), out.dalpha_in_prior = st.dalpha_in_prior; else, out.dalpha_in_prior = NaN; end
end
