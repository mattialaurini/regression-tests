function c = hill_inv(b, E0, gamma, Emax)
%HILL_INV Inverse Hill map of the paper, Eq. (inverse of g):
%   c(k) = ((E0 - b(k)) / (Emax - E0 + b(k)))^(1/gamma),
%   with c(k) = 0 for samples with b(k) = E0.
%   c is the effect-site concentration normalised by Ce50.
s = E0 - b(:);
if any(s < 0)
    error('hill_inv:range', 'All samples must satisfy b(k) <= E0.');
end
den = Emax - s;
pos = s > 0;
if any(den(pos) <= 0)
    error('hill_inv:range', 'Requires Emax - E0 + min_k b(k) > 0.');
end
c = zeros(size(s));
c(pos) = (s(pos) ./ den(pos)).^(1/gamma);
end
