function v = paper_prctile(x, p)
%PAPER_PRCTILE Percentile with linear interpolation (no toolbox needed).
x = sort(x(:)); n = numel(x);
if n == 0, v = NaN; return; end
r = 1 + (n - 1) * p / 100; f = floor(r);
v = x(f) + (r - f) * (x(min(n, f + 1)) - x(f));
end
