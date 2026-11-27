function o=nist_diagonal_options(base,spec)
o=base;
o.InitialApproxGridPoints=spec.M;
o.MaximumFineGridPoints=spec.Nhat;
o.CoarseGridPoints=2; % Unused by the no-m2 solver.
o.EpsilonRelative=spec.epsilon_relative;
o.EpsilonAbsolute=spec.epsilon_absolute;
o.TimeLimit=spec.target_time_limit;
end
