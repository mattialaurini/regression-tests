% TESTS_ML_THIRD  Validity of the third-order error model: at random points z
% of random boxes, the interpolation error of psi on the grid cell of z must
% satisfy, row by row,
%   | e_psi + sum_d theta_d psi_dd(cell centre) | <= sum_d theta_d sum_j (h_j/2) D3(:,d,j)
% with theta_d = (z_d - a_d)(b_d - z_d)/2. Expected: 0 violations.
sim = simulate_patient(1, 0.5); u = sim.u; rng(8, 'twister');
for cfg = {{[1;-3;3], [7;3;7.6]}, {[1;-3;3;3], [7;3;7.6;7.6]}}
    L0 = cfg{1}{1}; H0 = cfg{1}{2}; d = numel(L0); W = H0 - L0; nv = 0; nt = 0;
    for trial = 1:40
        zmin = L0 + (W - 1e-3) .* rand(d, 1); zmax = min(zmin + W .* rand(d, 1).^3 * 0.3 + 1e-3, H0);
        N = 2 + floor(2 * rand); h = (zmax - zmin) / (N - 1);
        D3 = ml_d3_bounds(u, zmin, zmax);
        for pt = 1:5
            z = zmin + (zmax - zmin) .* rand(d, 1);
            lowc = min(N - 1, floor((z - zmin) ./ h) + 1);
            a = zmin + (lowc - 1) .* h; b = a + h; tt = (z - a) ./ h;
            Li = 0;                                       % multilinear interpolant on the cell
            for c = 0:2^d - 1
                bits = bitget(c, 1:d)';
                wgt = prod(bits .* tt + (1 - bits) .* (1 - tt));
                Li = Li + wgt * ml_psi(a + bits .* h, u);
            end
            ep = ml_psi(z, u) - Li;
            theta = (z - a) .* (b - z) / 2;
            Mpt = ml_point_d2((a + b) / 2, u);
            lhs = abs(ep + Mpt * theta);
            rhs = zeros(numel(u), 1);
            for dd = 1:d, rhs = rhs + theta(dd) * (reshape(D3(:, dd, :), [], d) * (h / 2)); end
            nv = nv + sum(lhs > rhs + 1e-10 * (1 + abs(ep))); nt = nt + numel(lhs);
        end
    end
    fprintf('third-order model, d = %d: %d violations out of %d checks\n', d, nv, nt);
end
