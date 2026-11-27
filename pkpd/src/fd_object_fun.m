function [U, theta] = fd_object_fun(z, ge, V, Q)
%FD_OBJECT_FUN min over theta in C = conv(V) of ||F(z)[1; theta]||_2.
F = fd_matrix(ge, z, Q);
[U, lam, ok] = norm_min_lp(F * [ones(1, size(V, 2)); V], [], [], '2');
if ok, theta = V * lam; else, U = inf; theta = nan(size(V, 1), 1); end
end
