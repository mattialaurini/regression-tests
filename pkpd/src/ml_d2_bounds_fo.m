function [lo, hi, ex1] = ml_d2_bounds_fo(u, zmin, zmax, K)
%ML_D2_BOUNDS_FO Enclosures of the pure second derivatives of psi(k; z) for the
%   first-order model, z = (gamma, ell, w) (columns as in ML_D2_BOUNDS).
%   x = sum_j u(k-1-j) f_j(w), f_j = e^{-w}(1-e^{-w})^j > 0 (exact range of each
%   term); d log x/dw = E[a_j], a_j = -1 + j s(w); d2 log x/dw2 = E[a_j'] + Var(a_j),
%   a_j' = -j phi(w); weighted means over the active lags j in [jmin, jmax],
%   Var(a) <= ((jmax - jmin) s)^2/4. Chain rule: ML_CHAIN_BOUNDS.
if nargin < 4 || isempty(K), K = 4; end
u = u(:); n = numel(u);
if any(u < 0), error('ml_d2_bounds_fo:input', 'Requires a nonnegative input.'); end
wa = zmin(3); wb = zmax(3); G = [zmin(1), zmax(1)];
s = @(w) 1 ./ (exp(w) - 1); phi = @(w) exp(w) ./ (exp(w) - 1).^2;
sl = s(wb); sh = s(wa); fl = phi(wb); fh = phi(wa);
lo = zeros(n, 3); hi = zeros(n, 3); ex1.psi = zeros(n, 2); ex1.d1 = zeros(n, 3, 2);
LA = nan(n, 2); RA = LA; QA = LA;
for kk = 2:n
    k = kk - 1; j = (0:k-1)';
    ui = u(k - j); act = ui > 0;                 % u(k-1-j) -> index k-j
    if ~any(act), continue; end
    j = j(act); ua = ui(act);
    f = @(w) exp(-w) .* (1 - exp(-w)).^j;
    fa = f(wa); fb = f(wb); fmin = min(fa, fb); fmax = max(fa, fb);
    ws = log(j + 1); in = ws > wa & ws < wb;      % stationary point of f_j
    if any(in), fs = f(ws); fmax(in) = max(fmax(in), fs(in)); end
    xlo = ua' * fmin; xhi = ua' * fmax;
    j0 = min(j); j1 = max(j);
    RA(kk, :) = [-1 + j0 * sl, -1 + j1 * sh];
    QA(kk, :) = [-j1 * fh, -j0 * fl + ((j1 - j0) * sh)^2 / 4];
    LA(kk, :) = [zmin(2) + log(xlo), zmax(2) + log(xhi)];
end
act = find(~isnan(LA(:, 1)));
if isempty(act), return; end
[d2g, d2l, d2w, ex] = ml_chain_bounds(G, LA(act, :), {RA(act, :)}, {QA(act, :)}, K);
lo(act, :) = [d2g(:,1), d2l(:,1), d2w{1}(:,1)];
hi(act, :) = [d2g(:,2), d2l(:,2), d2w{1}(:,2)];
ex1.psi(act, :) = ex.psi;
D = {ex.d1g, ex.d1l, ex.d1w{1}};
for jj = 1:3, ex1.d1(act, jj, 1) = D{jj}(:, 1); ex1.d1(act, jj, 2) = D{jj}(:, 2); end
end
