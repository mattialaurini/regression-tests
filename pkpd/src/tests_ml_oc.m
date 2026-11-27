% TESTS_ML_OC  Checks for the one-compartment + effect-site ML model
% (z = (gamma, ell, w1, w2)). Expected: 0 violations, containment ~1e-15.
sim = simulate_patient(1, 0.5);
ge = struct('y', sim.b, 'u', sim.u, 'E0', sim.E0);
rng(11, 'twister');
zlo = [1; -3; 3; 3]; zhi = [7; 3; 7.6; 7.6]; W = zhi - zlo;
nviol = 0; ntest = 0;
for trial = 1:50
    zmin = zlo + (W - 1e-3) .* rand(4, 1) .* [1; 1; 1; 1];
    zmax = min(zmin + W .* rand(4, 1).^2 + 1e-3, zhi);
    [lo, hi] = ml_d2_bounds_oc(ge.u, zmin, zmax);
    for pt = 1:10
        z = zmin + (zmax - zmin) .* rand(4, 1);
        for dd = 1:4
            hstep = 1e-4 * (zmax(dd) - zmin(dd) + 1);
            e = zeros(4, 1); e(dd) = hstep;
            fd = (ml_psi(z + e, ge.u) - 2*ml_psi(z, ge.u) + ml_psi(z - e, ge.u)) / hstep^2;
            tol = 1e-4 * (1 + abs(fd));
            nviol = nviol + sum(fd < lo(:, dd) - tol | fd > hi(:, dd) + tol);
            ntest = ntest + numel(fd);
        end
    end
end
fprintf('ml_d2_bounds_oc: %d violations out of %d checks\n', nviol, ntest);
E = [50 160];
opts = struct('lb_safety', 0, 'return_matrices', true, 'no_solve', true);
worst = -inf;
for trial = 1:20
    zmin = zlo + (W - 1e-3) .* rand(4, 1);
    zmax = min(zmin + W .* rand(4, 1).^2 + 1e-3, zhi);
    N = 2 + floor(2*rand);
    [~, info] = ml_lower_bound(zmin, zmax, ge, E, N, N, opts);
    worst = max(worst, ml_containment_check(info, ge, zmin, zmax, N, E));
end
fprintf('ML (one compartment + effect site) containment: worst relative violation %.3g\n', worst);
