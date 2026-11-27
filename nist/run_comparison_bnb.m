function [results,summary]=run_comparison_bnb(varargin)
%RUN_COMPARISON_BNB The 24 BnB runs used for the SCIP comparison.
[results,summary]=reproduce_nist("bnb",varargin{:});
end
