function [results,summary]=run_ablation(varargin)
%RUN_ABLATION Matched grids, M=Nhat=3/5/9/17, with and without m2.
[results,summary]=reproduce_nist("ablation",varargin{:});
end
