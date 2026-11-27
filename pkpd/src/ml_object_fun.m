function [U, y] = ml_object_fun(z, ge, Erange)
%ML_OBJECT_FUN Gaussian ML criterion (output error in BIS space) at z, minimised over
%   the linear parameters on the same domain C used by ML_LOWER_BOUND:
%     E0 known:     U = min_{Emax in [E-, E+]} || (E0 - b) - Emax psi(z) ||_2,  y = Emax
%                   (projected least-squares value, exact);
%     E0 estimated (ge.E0range): U = min over (E0, Emax) in ge.E0range x [E-, E+],
%                   y = [E0; Emax] (exact, closed form: ML_LS_E0).
b = ge.y(:);
mdl = ''; if isfield(ge, 'mlmodel'), mdl = ge.mlmodel; end
psi = ml_psi(z, ge.u, mdl);
if isfield(ge, 'E0range') && ~isempty(ge.E0range)
    lo = [ge.E0range(1); Erange(1)]; hi = [ge.E0range(2); Erange(2)];
    [y, U] = ml_ls_e0(psi, b, lo, hi);                 % exact, closed form
else
    s = ge.E0 - b; pp = psi' * psi;
    if pp > 0, y = min(max((psi' * s) / pp, Erange(1)), Erange(2)); else, y = Erange(1); end
    U = norm(s - y * psi);
end
end
