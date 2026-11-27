function [results,summary]=run_comparison_scip(varargin)
%RUN_COMPARISON_SCIP The 24 SCIP runs, independently of BnB execution.
[results,summary]=reproduce_nist("scip",varargin{:});
end
