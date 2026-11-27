function v = gurobi_presolve(newv)
%GUROBI_PRESOLVE Gurobi Presolve parameter used by every solve of the bounds.
%   Default -1 (Gurobi's automatic choice). Measured with BENCH_PRESOLVE (MATLAB,
%   Gurobi 12): presolve off (0) changes the time per QP by less than 1% in total and the
%   complete searches by -6% to +8% (one-compartment ML model +8%), with identical
%   certified results: no consistent gain. The cost of a small QP (about 20 ms even with
%   3-7 columns) is mostly the overhead of the call, which presolve does not affect; the
%   hybrid solver mode reduces the number of calls instead.
%   gurobi_presolve(0) disables presolve; 1 and 2 are conservative and aggressive.
persistent pv
if nargin > 0, pv = newv; end
if isempty(pv), pv = -1; end
v = pv;
end
