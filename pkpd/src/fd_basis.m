function [Q, den] = fd_basis(demo, u, nb, Ts)
%FD_BASIS Orthonormal basis of the fixed-denominator (Schnider) model.
%   The denominator den is that of the full 4th-order Schnider model of the
%   demographics demo = [age, height, weight, gender] (ZOH, period Ts). The
%   effect-site concentration is modelled as c = (B(q)/den(q)) u with a free
%   numerator of nb terms; the columns (1/den) q^{-i} u, i = 1..nb, are nearly
%   collinear (poles close to 1), so they are orthonormalised over the data
%   window (SVD, dropping directions with singular value <= 1e-9 s_max):
%   c = Q theta. Same model space, well-conditioned theta.
if nargin < 4, Ts = 1; end
u = u(:); n = numel(u);
row = [1, demo(1), demo(2), demo(3), demo(4), 1, 1, 100, 100];
[k10, k12, k13, k21, k31, k1e, ke0, V1] = propofol_patient_parameters(1, row);
A = [-(k10 + k12 + k13), k21, k31, 0; k12, -k21, 0, 0; k13, 0, -k31, 0; k1e / V1, 0, 0, -ke0];
den = real(poly(eig(expm(A * Ts))));                 % discrete-time poles (ZOH)
Psi = zeros(n, nb);
for i = 1:nb
    Psi(:, i) = filter(1, den, [zeros(i, 1); u(1:n-i)]);
end
[Uq, Sq] = svd(Psi, 0); sv = diag(Sq);
Q = Uq(:, sv > 1e-9 * sv(1));          % orthonormal, rank-revealing
end
