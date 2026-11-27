function psi = ml_psi(z, u, model)
%ML_PSI Hill response psi(k) = c^gamma/(1+c^gamma), c = exp(ell)*x(k),
%   written as the logistic of t = gamma*(ell + log x).
%   z = (gamma, ell, w)       : double pole (ML_SIM_X), or first order (ML_SIM_X1)
%                               if model = 'firstorder'
%   z = (gamma, ell, w1, w2)  : one compartment + effect-site filter (ML_SIM_X2)
if nargin < 3, model = ''; end
if strcmp(model, 'firstorder')
    x = ml_sim_x1(u, z(3));
elseif numel(z) == 3
    x = ml_sim_x(u, z(3));
else
    x = ml_sim_x2(u, z(3), z(4));
end
psi = zeros(size(x));
pos = x > 0;
t = z(1) * (z(2) + log(x(pos)));
psi(pos) = 1 ./ (1 + exp(-t));
end
