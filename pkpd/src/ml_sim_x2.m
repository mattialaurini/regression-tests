function x = ml_sim_x2(u, w1, w2)
%ML_SIM_X2 One-compartment PK + effect-site filter, discrete time (ZOH),
%   each block with unit DC gain:
%     Cp(k+1) = p1 Cp(k) + (1-p1) u(k)      plasma,      p1 = 1 - exp(-w1)
%     x(k+1)  = p2 x(k)  + (1-p2) Cp(k)     effect site, p2 = 1 - exp(-w2)
%   x(0) = Cp(0) = 0. The gain (volume, Ce50) is a separate parameter.
%   Equivalently x(k) = sum_{i,l>=0} u(k-2-i-l) (1-p1) p1^i (1-p2) p2^l,
%   which is symmetric in (p1, p2).
% (implemented with FILTER: compartment, then effect site)
u = u(:); p1 = 1 - exp(-w1); p2 = 1 - exp(-w2);
x = filter([0, 1 - p2], [1, -p2], filter([0, 1 - p1], [1, -p1], u));
end
