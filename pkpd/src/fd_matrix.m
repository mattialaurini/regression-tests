function F = fd_matrix(ge, z, Q)
%FD_MATRIX Regression matrix of the fixed-denominator likelihood approximation:
%   r(theta, z) = F(z) [1; theta] = w(z) .* (c(z) - Q theta),
%   c = g_z^{-1}(b), w = slope of the Hill map at c (BIS linearization).
s = ge.E0 - ge.y(:);
c = hill_inv(ge.y, ge.E0, z(1), z(2));
w = zeros(size(s)); pos = s > 0;
w(pos) = z(1) * s(pos) .* (z(2) - s(pos)) ./ (z(2) * c(pos));
F = [w .* c, -w .* Q];
end
