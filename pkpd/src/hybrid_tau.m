function t = hybrid_tau(newt)
%HYBRID_TAU Current pruning threshold for the 'hybrid' solver mode (set by
%   BNB_SCREENING at every iteration; NaN = unknown, then no shortcut is taken).
persistent tau
if nargin > 0, tau = newt; end
if isempty(tau), tau = NaN; end
t = tau;
end
