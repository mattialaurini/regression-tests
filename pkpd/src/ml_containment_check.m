function worst = ml_containment_check(info, ge, zmin, zmax, N, V)
%ML_CONTAINMENT_CHECK Worst relative violation of C- v <= e <= C+ v, with
%   v_ab = lambda_a mu_b(z), e = f(y,z) - G v, at random points and grid nodes
%   (any number d of nonlinear parameters). With ge.E0range, the linear
%   parameters are (E0, Emax) on a rectangle (four vertices).
worst = -inf; d = numel(zmin);
if ~isfield(ge, 'E0range') || isempty(ge.E0range), s = ge.E0 - ge.y(:); end
for pt = 1:20
    z = zmin + (zmax - zmin) .* rand(d, 1);
    if mod(pt, 2) == 0
        for dd = 1:d, z(dd) = info.grid{dd}(randi(N)); end
    end
    mdl = ''; if isfield(ge, 'mlmodel'), mdl = ge.mlmodel; end
    lam = rand; y = lam*V(1) + (1-lam)*V(2);
    if isfield(ge, 'E0range') && ~isempty(ge.E0range)
        % (E0, Emax) in the rectangle: bilinear weights on the vertices ordered
        % (E0-,E-), (E0-,E+), (E0+,E-), (E0+,E+) as in ML_LOWER_BOUND
        a0 = rand; e0 = (1-a0)*ge.E0range(1) + a0*ge.E0range(2);
        lw = kron([1-a0; a0], [lam; 1-lam]);
        f = (e0 - ge.y(:)) - y * ml_psi(z, ge.u, mdl);
    else
        lw = [lam; 1-lam];
        f = s - y * ml_psi(z, ge.u, mdl);
    end
    idx = zeros(d, 1); tt = zeros(d, 1);
    for dd = 1:d
        g = info.grid{dd};
        idx(dd) = min(N-1, floor((z(dd) - g(1)) / (g(end) - g(1)) * (N-1)) + 1);
        tt(dd) = (z(dd) - g(idx(dd))) / (g(idx(dd)+1) - g(idx(dd)));
    end
    mu = zeros(N^d, 1);
    for c = 0:2^d-1
        bits = bitget(c, 1:d)';
        ii = idx + bits;
        wgt = prod(bits .* tt + (1 - bits) .* (1 - tt));
        b = 1 + sum((ii - 1) .* (N .^ (0:d-1))');
        mu(b) = mu(b) + wgt;
    end
    v = kron(mu, lw);
    e = f - info.G * v;
    worst = max(worst, max([info.Clo*v - e; e - info.Chi*v]) / max(abs(f)));
end
end
