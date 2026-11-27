function [x, val] = bounded_ls(Mm, dd, lo, hi)
%BOUNDED_LS  min ||dd + Mm x||_2  s.t.  lo <= x <= hi (bounded coordinates are the
%   finite entries of lo/hi), solved exactly by enumerating, for every bounded
%   coordinate, the three cases lower bound / upper bound / free. The free part is
%   solved in the minimum-norm least-squares sense, so the solution is well defined
%   also when Mm is rank deficient.
%   Implementation: one thin QR factorisation Mm = Q R per call; then
%   ||dd + Mm x||^2 = ||dd||^2 - ||c||^2 + ||c + R x||^2 with c = Q'dd, so every case
%   works with the small factor R (size k x k) instead of the n x k matrix, without
%   squaring the condition number. Fast path: if the unconstrained solution satisfies
%   the bounds it is returned at once.
nv = size(Mm, 2); bd = find(~isinf(lo) | ~isinf(hi));
[Qm, R] = qr(Mm, 0); c = Qm' * dd;
r0 = max(dd' * dd - c' * c, 0);                    % part of the residual orthogonal to range(Mm)
x0 = -pinv(R) * c;
if all(x0 >= lo - 1e-12 & x0 <= hi + 1e-12)
    x = x0; val = sqrt(r0 + norm(c + R * x0)^2); return
end
val = inf; x = zeros(nv, 1);
for code = 0:3^numel(bd) - 1
    st = mod(floor(code ./ 3.^(0:numel(bd)-1)), 3);     % 0 free, 1 lower, 2 upper
    xx = zeros(nv, 1); fixd = false(nv, 1);
    for j = 1:numel(bd)
        if st(j) == 1, xx(bd(j)) = lo(bd(j)); fixd(bd(j)) = true; end
        if st(j) == 2, xx(bd(j)) = hi(bd(j)); fixd(bd(j)) = true; end
    end
    fr = ~fixd;
    if any(fr), xx(fr) = -pinv(R(:, fr)) * (c + R(:, fixd) * xx(fixd)); end
    if all(xx >= lo - 1e-12 & xx <= hi + 1e-12)
        v = sqrt(r0 + norm(c + R * xx)^2);
        if v < val, val = v; x = xx; end
    end
end
end
