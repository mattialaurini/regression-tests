function D2 = ml_point_d2(z, u)
%ML_POINT_D2 Exact (up to rounding) pure second derivatives of psi(k; z) at a
%   point z, by forward-mode differentiation of the filter recursions.
%   z = (gamma, ell, w) [double pole] or (gamma, ell, w1, w2) [one compartment +
%   effect site]. Output n x d: d2psi/dz_d^2 for every sample k.
%   psi = sigma(t), t = gamma (ell + log x):
%     psi_gg = sigma'' L^2, psi_ll = sigma'' gamma^2,
%     psi_ww = sigma'' (gamma x'/x)^2 + sigma' gamma (x''/x - (x'/x)^2).
u = u(:); n = numel(u); d = numel(z);
g = z(1); ell = z(2);
if d == 3
    p = 1 - exp(-z(3)); dp = 1 - p; ddp = -(1 - p);
    x1 = 0; a1 = 0; b1 = 0; x = 0; a = 0; b = 0;     % values, d/dw, d2/dw2
    X = zeros(n, 1); X1 = X; X2 = X;
    for k = 1:n
        X(k) = x; X1(k) = a; X2(k) = b;
        x1n = u(k) + p * (x1 - u(k));
        a1n = dp * (x1 - u(k)) + p * a1;
        b1n = ddp * (x1 - u(k)) + 2 * dp * a1 + p * b1;
        xn = x1 + p * (x - x1);
        an = a1 + dp * (x - x1) + p * (a - a1);
        bn = b1 + ddp * (x - x1) + 2 * dp * (a - a1) + p * (b - b1);
        x1 = x1n; a1 = a1n; b1 = b1n; x = xn; a = an; b = bn;
    end
    Dw = {X1}; DDw = {X2};
else
    p1 = 1 - exp(-z(3)); dp1 = 1 - p1; ddp1 = -(1 - p1);
    p2 = 1 - exp(-z(4)); dp2 = 1 - p2; ddp2 = -(1 - p2);
    cp = 0; c1 = 0; c11 = 0;                 % plasma: value, d/dw1, d2/dw1^2
    x = 0; x1 = 0; x11 = 0; x2 = 0; x22 = 0; % effect site and derivatives
    X = zeros(n, 1); XA = X; XAA = X; XB = X; XBB = X;
    for k = 1:n
        X(k) = x; XA(k) = x1; XAA(k) = x11; XB(k) = x2; XBB(k) = x22;
        cpn = u(k) + p1 * (cp - u(k));
        c1n = dp1 * (cp - u(k)) + p1 * c1;
        c11n = ddp1 * (cp - u(k)) + 2 * dp1 * c1 + p1 * c11;
        xn = cp + p2 * (x - cp);
        x1n = (1 - p2) * c1 + p2 * x1;
        x11n = (1 - p2) * c11 + p2 * x11;
        x2n = dp2 * (x - cp) + p2 * x2;
        x22n = ddp2 * (x - cp) + 2 * dp2 * x2 + p2 * x22;
        cp = cpn; c1 = c1n; c11 = c11n; x = xn; x1 = x1n; x11 = x11n; x2 = x2n; x22 = x22n;
    end
    Dw = {XA, XB}; DDw = {XAA, XBB};
end
D2 = zeros(n, d);
pos = X > 0;
L = ell + log(X(pos));
t = g * L;
sg = 1 ./ (1 + exp(-t)); s1 = sg .* (1 - sg); s2 = s1 .* (1 - 2 * sg);
D2(pos, 1) = s2 .* L.^2;
D2(pos, 2) = s2 * g^2;
for i = 1:numel(Dw)
    r1 = Dw{i}(pos) ./ X(pos); r2 = DDw{i}(pos) ./ X(pos);
    D2(pos, 2 + i) = s2 .* (g * r1).^2 + s1 * g .* (r2 - r1.^2);
end
end
