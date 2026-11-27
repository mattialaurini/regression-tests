function P = halton_points(K, d)
%HALTON_POINTS First K points (skipping 0) of the d-dimensional Halton
%   sequence in [0,1]^d (rows). Deterministic multi-start points.
pr = [2 3 5 7 11 13 17 19 23 29];
P = zeros(K, d);
for j = 1:d
    b = pr(j);
    for i = 1:K
        f = 1; r = 0; n = i;
        while n > 0
            f = f / b; r = r + f * mod(n, b); n = floor(n / b);
        end
        P(i, j) = r;
    end
end
end
