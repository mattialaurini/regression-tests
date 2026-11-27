function v = getfield_or(s, f, default)
%GETFIELD_OR Field f of struct s, or default if s has no such field (e.g. results files
%   written before the field existed: presolve was then Gurobi's automatic choice, -1).
if isfield(s, f), v = s.(f); else, v = default; end
end
