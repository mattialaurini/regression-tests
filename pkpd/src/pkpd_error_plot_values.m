function [displayed, zeroLevel] = pkpd_error_plot_values(values, zeroLevel)
%PKPD_ERROR_PLOT_VALUES Move exact zeros to a positive level for log plotting.
% Positive errors are returned unchanged; statistics must use the original data.
% Pass a shared zeroLevel to keep all panels on the same plotting convention.
if nargin < 2
    positive = values(isfinite(values) & values > 0);
    if isempty(positive), zeroLevel = 1e-12;
    else, zeroLevel = min(positive) / 10; end
end
assert(isscalar(zeroLevel) && isfinite(zeroLevel) && zeroLevel > 0, ...
    'The zero-error plotting level must be finite and strictly positive.');
displayed = values;
displayed(values == 0) = zeroLevel;
end
