function x = ml_sim_x1(u, w)
%ML_SIM_X1 First-order filtered input (unit DC gain), pole p = 1 - exp(-w):
%   x(k+1) = p x(k) + (1-p) u(k), x(0) = 0;
%   x(k) = sum_{j>=0} f_j(w) u(k-1-j),  f_j = e^{-w} (1 - e^{-w})^j.
% (implemented with FILTER)
u = u(:); p = 1 - exp(-w);
x = filter([0, 1 - p], [1, -p], u);
end
