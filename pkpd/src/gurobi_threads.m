function n = gurobi_threads(newn)
%GUROBI_THREADS Number of threads used by every Gurobi solve (0 = Gurobi's automatic
%   multi-threading).
%   Default: 0 (all cores) in a normal MATLAB session, where the runs execute one after
%   another; 1 inside a parallel worker (parfor), because Gurobi's default in every
%   worker oversubscribes the machine (identical runs were measured 3.3 times slower,
%   median). gurobi_threads(k) sets the value for the current session (or worker).
persistent nt
if nargin > 0, nt = newn; end
if isempty(nt)
    nt = 0;
    try
        if ~isempty(getCurrentTask()), nt = 1; end      % inside a parallel worker
    catch
    end
end
n = nt;
end
