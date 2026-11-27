function sim = simulate_patient(id, noise, varargin)
% Reproduction note: draw the seven PK perturbations before measurement noise.
% Keeping that order and seed=id is essential for the shared datasets.
% Input u is a propofol mass flow in mg/s, sampled every Ts seconds.
% Default input: 10 mg/s for 0<=t<10 s, 3 mg/s for 10<=t<35 s, then zero.
% Mismatch perturbs only the seven PK quantities below. The tabulated
% Hill-map parameters remain fixed when constructing the simulated response.
%SIMULATE_PATIENT Simulated induction of patient 'id' (Table data_patients):
%   Schnider PK model, ZOH with period Ts, input of Eq. (input_u),
%   Hill map, additive i.i.d. N(0, noise^2) on the BIS samples,
%   noise realisation rng(seed,'twister') (default seed = id). No clipping.
%   omega > 0: the TRUE PK model is Schnider with every parameter (k10, k12, k13,
%   k21, k31, ke0 = k1e, V1) multiplied by exp(omega*xi), xi ~ N(0,1) (drawn
%   after rng(seed), before the noise).
d = struct('Ts', 1, 'horizon', 150, 'u', [], 'seed', [], 'omega', 0);
o = parse_opts(d, varargin{:});
if isempty(o.seed), o.seed = id; end
patient_profiles
[k10, k12, k13, k21, k31, k1e, ke0, V1, gamma_ob, E0, Emax_ob, Ce50] = ...
    propofol_patient_parameters(id, Patients);
rng(o.seed, 'twister');
pk = [k10, k12, k13, k21, k31, ke0, V1];
if o.omega > 0, pk = pk .* exp(o.omega * randn(1, 7)); end   % (no draw if omega = 0)
k10 = pk(1); k12 = pk(2); k13 = pk(3); k21 = pk(4); k31 = pk(5); ke0 = pk(6); k1e = pk(6); V1 = pk(7);
t = (0:o.horizon - 1)' * o.Ts;
if isempty(o.u), u = 10 * (t < 10) + 3 * (t >= 10 & t < 35); else, u = o.u(:); end
A = [-(k10 + k12 + k13), k21, k31, 0; k12, -k21, 0, 0; k13, 0, -k31, 0; k1e / V1, 0, 0, -ke0];
Pd = c2d(ss(A, [1; 0; 0; 0], [0, 0, 0, 1], 0), o.Ts, 'zoh');
[Ad, Bd, Cd] = ssdata(Pd);
x = zeros(4, 1); Ce = zeros(numel(u), 1);
for k = 1:numel(u)
    Ce(k) = Cd * x;
    x = Ad * x + Bd * u(k);
end
c = Ce / Ce50;
b0 = E0 - Emax_ob * c.^gamma_ob ./ (1 + c.^gamma_ob);
sim.b = b0 + noise * randn(size(b0));
sim.b_noiseless = b0; sim.u = u; sim.E0 = E0; sim.seed = o.seed;
sim.pk = pk; sim.gamma_ob = gamma_ob; sim.Emax_ob = Emax_ob; sim.Ce50 = Ce50; sim.Ce = Ce;
end
