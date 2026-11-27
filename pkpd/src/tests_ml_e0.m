% TESTS_ML_E0  Containment of the ML relaxation with E0 estimated (linear parameters
% (E0, Emax) on a rectangle, four vertices), for the double-pole and one-compartment
% models, with the interval error model. Expected: worst violation ~1e-15.
sim = simulate_patient(1, 0.5, 'omega', 0.2);
rng(21, 'twister'); E = [50 160];
for mdl = {'doublepole', 'onecomp'}
    ge = struct('y', sim.b, 'u', sim.u, 'E0', sim.E0, 'mlmodel', mdl{1}, 'E0range', max(sim.b) + [-10, 5]);
    if strcmp(mdl{1}, 'onecomp'), L0 = [1; -3; 3; 3]; H0 = [7; 3; 7.6; 7.6]; else, L0 = [1; -3; 3]; H0 = [7; 3; 7.6]; end
    W = H0 - L0; worst = -inf; d = numel(L0);
    for trial = 1:15
        zmin = L0 + (W - 1e-3) .* rand(d, 1); zmax = min(zmin + W .* rand(d, 1).^2 + 1e-3, H0);
        N = 2 + floor(2 * rand);
        [~, info] = ml_lower_bound(zmin, zmax, ge, E, N, N, struct('lb_safety', 0, 'return_matrices', true, 'no_solve', true));
        worst = max(worst, ml_containment_check(info, ge, zmin, zmax, N, E));
    end
    fprintf('ML containment with E0 estimated (%s): worst relative violation %.3g\n', mdl{1}, worst);
end
