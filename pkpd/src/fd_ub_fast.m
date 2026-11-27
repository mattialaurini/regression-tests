function [U, theta, lb_ls] = fd_ub_fast(z, ge, V, Q, Vinv)
%FD_UB_FAST Exact min over theta in C = conv(V) of ||F(z)[1; theta]||_2, cheap when
%   the unconstrained least-squares solution lies in C:
%   theta_ls = argmin over all theta; if theta_ls is in the simplex C (barycentric
%   coordinates Vinv*[theta; 1] >= 0, V with n+1 columns) it is also the minimiser
%   over C; otherwise the QP is solved (FD_OBJECT_FUN).
%   lb_ls = ||F(z)[1; theta_ls]|| <= U (minimum over a larger set): a valid lower
%   bound of U, used to skip exact evaluations that cannot improve the incumbent.
F = fd_matrix(ge, z, Q);
th = -pinv(F(:, 2:end)) * F(:, 1);
lb_ls = norm(F(:, 1) + F(:, 2:end) * th);
if ~isempty(Vinv) && all(Vinv * [th; 1] >= -1e-12)
    U = lb_ls; theta = th;
else
    [U, theta] = fd_object_fun(z, ge, V, Q);
end
end
