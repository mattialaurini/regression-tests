function x = ml_sim_x(u, w)
%ML_SIM_X Double-pole filtered input (unit DC gain), pole p = 1 - exp(-w):
%   x1(k+1) = p x1(k) + (1-p) u(k),  x(k+1) = p x(k) + (1-p) x1(k),  x(0) = x1(0) = 0.
%   Equivalently x(k) = sum_{j>=1} h_j u(k-1-j),  h_j = j (1-p)^2 p^(j-1).
% (implemented with FILTER: two cascaded first-order recursions)
u = u(:); p = 1 - exp(-w);
x = filter([0, 1 - p], [1, -p], filter([0, 1 - p], [1, -p], u));
end
