function E = lw_enclosures(s, g1, g2, E1, E2, Ksub)
%LW_ENCLOSURES Interval enclosures over gamma in [g1,g2], Emax in [E1,E2] of
%   c = (s/(E-s))^(1/gamma) and of the BIS-linearization weight
%   w = gamma s (E-s)/(E c), with their first and pure second derivatives
%   (fields C, Cg, CE, Cgg, CEE, W, Wg, Wgg, WE, WEE; n x 2 each; zero where s = 0):
%     c_g = -c L/g^2, c_E = -c/(g(E-s)), c_gg, c_EE (HILL_INV_D2_BOUNDS),
%     w_g = w (g+L)/g^2, w_gg = w L^2/g^4, w_E = w lam, w_EE = w (lam^2 + mu),
%     L = log(s/(E-s)), lam = (s g + E)/(g E (E-s)),
%     mu = -s(2E-s)/(E^2 (E-s)^2) - 1/(g (E-s)^2).
if nargin < 6, Ksub = 4; end
s = s(:); n = numel(s); Z2 = zeros(n, 2);
E = struct('C', Z2, 'Cg', Z2, 'CE', Z2, 'Cgg', Z2, 'CEE', Z2, 'W', Z2, 'Wg', Z2, 'Wgg', Z2, 'WE', Z2, 'WEE', Z2);
pos = s > 0; sp = s(pos); m = numel(sp);
if m == 0, return; end
rlo = sp ./ (E2 - sp); rhi = sp ./ (E1 - sp);
cc = [rlo.^(1/g1), rlo.^(1/g2), rhi.^(1/g1), rhi.^(1/g2)];
clo = min(cc, [], 2); chi = max(cc, [], 2);
Lr = [log(rlo), log(rhi)]; Gi = repmat([g1, g2], m, 1);
E.C(pos, :) = [clo, chi];
E.Cg(pos, :) = imul(imul([clo, chi], -fliplr(Lr)), repmat([1/g2^2, 1/g1^2], m, 1));
E.CE(pos, :) = -fliplr(imul(imul([clo, chi], repmat([1/g2, 1/g1], m, 1)), [1 ./ (E2 - sp), 1 ./ (E1 - sp)]));
[ggl, ggh, eel, eeh] = hill_inv_d2_bounds(s, g1, g2, E1, E2, Ksub);
E.Cgg(pos, :) = [ggl(pos), ggh(pos)]; E.CEE(pos, :) = [eel(pos), eeh(pos)];
w0 = imul(imul(Gi, [sp .* (1 - sp / E1), sp .* (1 - sp / E2)]), [1 ./ chi, 1 ./ clo]);
E.W(pos, :) = w0;
E.Wg(pos, :) = imul(imul(w0, Lr + Gi), repmat([1/g2^2, 1/g1^2], m, 1));
E.Wgg(pos, :) = imul(imul(w0, isq(Lr)), repmat([1/g2^4, 1/g1^4], m, 1));
lam = [(sp * g1 + E1) ./ (g2 * E2 * (E2 - sp)), (sp * g2 + E2) ./ (g1 * E1 * (E1 - sp))];
t1 = [sp .* (2*E1 - sp) ./ (E2^2 * (E2 - sp).^2), sp .* (2*E2 - sp) ./ (E1^2 * (E1 - sp).^2)];
t2 = [1 ./ (g2 * (E2 - sp).^2), 1 ./ (g1 * (E1 - sp).^2)];
mu = -fliplr(t1 + t2);
E.WE(pos, :) = imul(w0, lam);
E.WEE(pos, :) = imul(w0, isq(lam) + mu);
end

function r = imul(a, b)
P = [a(:,1).*b(:,1), a(:,1).*b(:,2), a(:,2).*b(:,1), a(:,2).*b(:,2)];
r = [min(P, [], 2), max(P, [], 2)];
end
function r = isq(a)
q = [a(:,1).^2, a(:,2).^2]; r = [min(q, [], 2), max(q, [], 2)];
r(a(:,1) <= 0 & a(:,2) >= 0, 1) = 0;
end
