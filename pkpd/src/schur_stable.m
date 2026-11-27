function st = schur_stable(a)
%SCHUR_STABLE True if all roots of the polynomial a(1) z^n + ... + a(n+1) lie strictly
%   inside the unit circle (Schur-Cohn recursion on the reflection coefficients;
%   O(n^2), no root computation).
a = a(:)' / a(1);
while numel(a) > 1
    k = a(end);
    if abs(k) >= 1, st = false; return; end
    a = (a - k * fliplr(a)) / (1 - k^2);
    a = a(1:end-1);
end
st = true;
end
