% TESTS_FD  Containment test for the fixed-denominator relaxation (FD_LOWER_BOUND):
% for random (z, theta), f = F(z)[1;theta] = G v + e with v_ab = lambda_a mu_b(z)
% and Clo v <= e <= Chi v. Expected: worst violation <= ~1e-12.
patient_profiles
sim = simulate_patient(1, 0.5); b = min(sim.b, sim.E0);
ge = struct('y', b, 'u', sim.u, 'E0', sim.E0);
Q = fd_basis(Patients(1, 2:5), sim.u, 4);
R = 2 * norm(hill_inv(b, sim.E0, 2, 100));
V = R * (2 * (dec2bin(0:15) - '0')' - 1); p = 16;
rng(9, 'twister'); worst = -inf; Emn = 1.02 * max(ge.E0 - b);
for trial = 1:25
    g1 = 1 + 5*rand; g2 = min(7, g1 + 0.02 + 2*rand^2);
    E1 = Emn + (160 - Emn) * rand * 0.8; E2 = min(160, E1 + 0.5 + 30*rand^2);
    N = 2 + floor(3*rand);
    [~, info] = fd_lower_bound([g1; E1], [g2; E2], ge, V, Q, N, N, ...
        struct('lb_safety', 0, 'no_solve', true, 'return_matrices', true));
    for pt = 1:20
        z = [g1 + (g2-g1)*rand; E1 + (E2-E1)*rand];
        if mod(pt, 2) == 0, z = [info.gv(randi(N)); info.Ev(randi(N))]; end
        lam = -log(rand(p, 1)); lam = lam / sum(lam);
        f = fd_matrix(ge, z, Q) * [ones(1, p); V] * lam;
        ig = min(N-1, floor((z(1)-g1)/(g2-g1)*(N-1)) + 1); iE = min(N-1, floor((z(2)-E1)/(E2-E1)*(N-1)) + 1);
        tg = (z(1) - info.gv(ig)) / (info.gv(ig+1) - info.gv(ig)); tE = (z(2) - info.Ev(iE)) / (info.Ev(iE+1) - info.Ev(iE));
        mu = zeros(N^2, 1); bi = @(a, c) (c-1)*N + a;
        mu(bi(ig, iE)) = (1-tg)*(1-tE); mu(bi(ig+1, iE)) = tg*(1-tE); mu(bi(ig, iE+1)) = (1-tg)*tE; mu(bi(ig+1, iE+1)) = tg*tE;
        v = kron(mu, lam); e = f - info.G * v;
        worst = max(worst, max([info.Clo*v - e; e - info.Chi*v]) / max(1e-12, max(abs(f))));
    end
end
fprintf('fixed-denominator containment: worst relative violation %.3g\n', worst);
