function [toleranceOnly,atIncumbent,eligible] = nist_tolerance_prune_masks(lower,upper,relative,absolute)
%NIST_TOLERANCE_PRUNE_MASKS Classify with the incumbent AT DISCARD TIME.
% Counts are based on valid lower bounds, never on m2 priorities.
% The upper endpoint of the tolerance-only interval is strictly excluded.
arguments
    lower double
    upper double
    relative (1,1) double {mustBeNonnegative,mustBeLessThan(relative,1)}
    absolute (1,1) double {mustBeNonnegative}
end
threshold=(1-relative).*upper-absolute;
eligible=isfinite(upper) & ~isnan(lower) & lower>=threshold;
toleranceOnly=eligible & lower<upper;
atIncumbent=eligible & lower>=upper;
end
