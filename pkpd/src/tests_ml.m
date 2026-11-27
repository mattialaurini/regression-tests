% TESTS_ML  Checks of the Gaussian-ML relaxation (run in MATLAB with Gurobi).
% 1) ml_d2_bounds encloses the second derivatives of psi (finite differences).
% 2) Separable containment: f(y,z) = G v + e with C- v <= e <= C+ v
%    for random (z, y) in the region, including grid nodes (needs no solver).
% 3) m1(R,N) <= U(z) at random points; m1(R,M) <= m2(R,N,M) (needs Gurobi).
% Expected: 0 violations, worst containment violation <= ~1e-12, two 1's.
% Set skip_solver = true before running to skip part 3.
if ~exist('skip_solver', 'var'), skip_solver = false; end
sim = simulate_patient(1, 0.5);
ge = struct('y', sim.b, 'u', sim.u, 'E0', sim.E0);
rng(7, 'twister');
nviol = 0; ntest = 0;
for trial = 1:60
    zmin = [1 + 5*rand; -3 + 5*rand; 3 + 4*rand];
    wd = [6; 6; 4.6] .* rand(3, 1).^2 .* [1; 1; 1];
    zmax = min(zmin + wd + 1e-3, [7; 3; 7.6]);
    [lo, hi] = ml_d2_bounds(ge.u, zmin, zmax);
    for pt = 1:15
        z = zmin + (zmax - zmin) .* rand(3, 1);
        for dd = 1:3
            hstep = 1e-4 * (zmax(dd) - zmin(dd) + 1);
            e = zeros(3, 1); e(dd) = hstep;
            fd = (ml_psi(z + e, ge.u) - 2*ml_psi(z, ge.u) + ml_psi(z - e, ge.u)) / hstep^2;
            tol = 1e-4 * (1 + abs(fd));
            bad = fd < lo(:, dd) - tol | fd > hi(:, dd) + tol;
            nviol = nviol + sum(bad); ntest = ntest + numel(fd);
        end
    end
end
fprintf('ml_d2_bounds: %d violations out of %d checks\n', nviol, ntest);

E = [50 160]; V = E;
opts = struct('lb_safety', 0, 'return_matrices', true, 'no_solve', true);
worst = -inf;
for trial = 1:25
        zmin = [1 + 5*rand; -3 + 5*rand; 3 + 4*rand];
        zmax = min(zmin + [6; 6; 4.6] .* rand(3, 1).^2 + 1e-3, [7; 3; 7.6]);
        N = 2 + floor(3*rand);
        [~, info] = ml_lower_bound(zmin, zmax, ge, E, N, N, opts);
    worst = max(worst, ml_containment_check(info, ge, zmin, zmax, N, V));
end
fprintf('ML containment: worst relative violation %.3g\n', worst);

if ~skip_solver
    opts = struct('lb_safety', 0);
    ok1 = true; ok4 = true;
    for trial = 1:4
        zmin = [1.5 + 4*rand; -1 + 2*rand; 4 + 1.5*rand];
        zmax = zmin + [0.5; 0.4; 0.3] .* (0.2 + rand(3, 1));
        m1 = ml_lower_bound(zmin, zmax, ge, E, 5, 5, opts);
        for pt = 1:10
            ok1 = ok1 && (m1 <= ml_object_fun(zmin + (zmax - zmin) .* rand(3, 1), ge, E) + 1e-8);
        end
        m2 = ml_lower_bound(zmin, zmax, ge, E, 3, 5, opts);
        ok4 = ok4 && (m1 <= m2 + 1e-8);
    end
    fprintf('Prop. 5(1) m1 <= f: %d,  5(4) m1(R,M) <= m2(R,N,M): %d\n', ok1, ok4);
end

