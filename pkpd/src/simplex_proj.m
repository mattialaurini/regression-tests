function x = simplex_proj(y)
%SIMPLEX_PROJ Euclidean projection of y onto the unit simplex {x >= 0, sum x = 1}.
y = y(:); n = numel(y);
u = sort(y, 'descend'); cs = cumsum(u);
k = find(u - (cs - 1) ./ (1:n)' > 0, 1, 'last');
tau = (cs(k) - 1) / k;
x = max(y - tau, 0);
end
