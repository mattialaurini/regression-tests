function [z, U, y] = local_search(fun, z0, zmin, zmax, maxevals)
%LOCAL_SEARCH Optional local search of the routine r(R) (Section 5).
%   Nelder-Mead in coordinates normalised to the region box [zmin, zmax],
%   started from z0, with tight tolerances. fun(z) returns [U, y].
%   Points outside the box are evaluated at their projection, plus a
%   penalty, so every returned point is feasible (z in [zmin, zmax]).
%   The default fminsearch tolerances (1e-4) are of the same order as the
%   PK/PD objective values and stopped the search prematurely.
if nargin < 5 || isempty(maxevals), maxevals = 300; end
zmin = zmin(:); zmax = zmax(:); w = zmax - zmin;
t0 = (z0(:) - zmin) ./ w;
pen = @(t) 1e3 * sum(max(0, -t) + max(0, t - 1));
g = @(t) fun(zmin + min(max(t(:), 0), 1) .* w) + pen(t(:));
opt = optimset('Display', 'off', 'TolX', 1e-10, 'TolFun', 1e-14, ...
    'MaxFunEvals', maxevals, 'MaxIter', maxevals, 'OutputFcn', []);
t = fminsearch(g, t0, opt);
z = zmin + min(max(t(:), 0), 1) .* w;
[U, y] = fun(z);
end
