function [Ahat, Bhat] = schnider_ahat(demo, nc, u, Ts, pk)
%SCHNIDER_AHAT Prefilter denominator from the Schnider PK model.
%   demo = [age, height, weight, gender] (gender: 1 = m, 2 = f), nc = ARX order,
%   u = input samples, Ts = sampling period. The effect-site response of the
%   Schnider model (ZOH) to u is simulated without noise and an output-error
%   model B/A of order nc is fitted to it (equation-error start, Nelder-Mead
%   on the output error); Ahat = [1, a_1, ..., a_nc] is its (stabilised)
%   denominator and Bhat = [0, b_1, ..., b_nc] the numerator (response normalised
%   to unit peak); for nc = 4 the exact Schnider transfer function is returned.
%   Optional pk = [k10 k12 k13 k21 k31 ke0 V1] replaces the
%   demographic Schnider parameters (used for perturbed models).
%   Only the population model and the input are used, not the BIS.
if nargin < 4, Ts = 1; end
u = u(:); n = numel(u);
if nargin >= 5 && ~isempty(pk)
    k10 = pk(1); k12 = pk(2); k13 = pk(3); k21 = pk(4); k31 = pk(5); ke0 = pk(6); k1e = pk(6); V1 = pk(7);
else
    row = [1, demo(1), demo(2), demo(3), demo(4), 1, 1, 100, 100];   % PD values irrelevant
    [k10, k12, k13, k21, k31, k1e, ke0, V1] = propofol_patient_parameters(1, row);
end
A = [-(k10 + k12 + k13), k21, k31, 0; k12, -k21, 0, 0; k13, 0, -k31, 0; k1e / V1, 0, 0, -ke0];
sysd = c2d(ss(A, [1; 0; 0; 0], [0, 0, 0, 1], 0), Ts, 'zoh');
[Ad, Bd, Cd] = ssdata(sysd);
if nc == 4
    % order 4: the Schnider model itself (exact discrete transfer function)
    [num, den] = tfdata(tf(sysd), 'v');
    num = [zeros(1, numel(den) - numel(num)), num(:)'];  % align with den (strictly proper)
    cpk = max(filter(num, den, u(:)));
    Ahat = den(:)' / den(1); Bhat = num / den(1) / cpk; Bhat(1) = 0;
    return
end
x = zeros(4, 1); c = zeros(n, 1);
for k = 1:n, c(k) = Cd * x; x = Ad * x + Bd * u(k); end
c = c / max(c);
X = zeros(n, 2 * nc);
for i = 1:nc
    X(i+1:end, i) = -c(1:end-i); X(i+1:end, nc+i) = u(1:end-i);
end
th = X \ c;                                             % equation-error start
% Steiglitz-McBride iterations (linear least squares on prefiltered data)
for it = 1:20
    a = stab([1; th(1:nc)]);
    cf = filter(1, a, c); uf = filter(1, a, u);
    Xf = zeros(n, 2 * nc);
    for i = 1:nc
        Xf(i+1:end, i) = -cf(1:end-i); Xf(i+1:end, nc+i) = uf(1:end-i);
    end
    thn = Xf \ cf;
    if norm(thn - th) <= 1e-12 * norm(th), th = thn; break; end
    th = thn;
end
% short output-error polish; stability by the Schur-Cohn test (no roots)
oe = @(t) oe_cost(t, nc, u, c);
opt = optimset('Display', 'off', 'TolX', 1e-14, 'TolFun', 1e-20, 'MaxFunEvals', 3000, 'MaxIter', 3000);
if ~schur_stable([1; th(1:nc)]), as = stab([1; th(1:nc)]); th(1:nc) = as(2:end); end
th = fminsearch(oe, th, opt);
Ahat = stab([1; th(1:nc)]);
Bhat = [0, th(nc+1:end)'];
end

function v = oe_cost(t, nc, u, c)
a = [1; t(1:nc)];
if ~schur_stable(a), v = 1e10; return; end
v = sum((filter([0; t(nc+1:end)]', a', u) - c).^2);
end

function a = stab(a)
r = roots(a(:)'); out = abs(r) >= 0.9999; r(out) = 0.9999 ./ conj(r(out));
a = real(poly(r));
end
