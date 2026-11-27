function [y, U] = ml_ls_e0(psi, b, lo, hi)
%ML_LS_E0 Exact  min ||E0*1 - b - Emax*psi||_2  over the rectangle lo <= (E0, Emax) <= hi,
%   in closed form. Unconstrained solution: Emax = -psi_c'b_c / psi_c'psi_c (centred
%   vectors), E0 = mean(b) + Emax mean(psi). If it is outside the rectangle, the
%   minimiser of this convex quadratic lies on an edge; on each edge the problem is a
%   one-dimensional least squares, solved and projected exactly. y = [E0; Emax].
b = b(:); psi = psi(:); n = numel(b);
f = @(e0, em) norm(e0 - b - em * psi);
pm = sum(psi) / n; bm = sum(b) / n; pc = psi - pm; bc = b - bm; pp = pc' * pc;
if pp > 0
    em = -(pc' * bc) / pp; e0 = bm + em * pm;
    if e0 >= lo(1) && e0 <= hi(1) && em >= lo(2) && em <= hi(2)
        y = [e0; em]; U = f(e0, em); return
    end
end
cand = zeros(2, 4); ps = psi' * psi;
for i = 1:2                                  % E0 on its bounds
    e0 = [lo(1), hi(1)]; e0 = e0(i);
    em = min(max(-(psi' * (b - e0)) / max(ps, realmin), lo(2)), hi(2));
    cand(:, i) = [e0; em];
end
for i = 1:2                                  % Emax on its bounds
    em = [lo(2), hi(2)]; em = em(i);
    cand(:, 2 + i) = [min(max(bm + em * pm, lo(1)), hi(1)); em];
end
vals = arrayfun(@(i) f(cand(1, i), cand(2, i)), 1:4);
[U, k] = min(vals); y = cand(:, k);
end
