function options = nist_bnb_options(varargin)
%NIST_BNB_OPTIONS Options for the interpolation-based NIST BnB.
p = inputParser;
addParameter(p,"TimeLimit",600);
addParameter(p,"NodeLimit",Inf);
addParameter(p,"EpsilonRelative",0.01);
addParameter(p,"EpsilonAbsolute",1e-8);
addParameter(p,"CoarseGridPoints",2);
addParameter(p,"InitialApproxGridPoints",5);
addParameter(p,"MaximumFineGridPoints",3);
addParameter(p,"StrongBranchingTieBreak",false);
addParameter(p,"BranchingRule","normalized_longest_edge");
addParameter(p,"TieTolerance",1e-12);
addParameter(p,"RecordRelativeTargets",[0.1 0.05 0.01]);
addParameter(p,"Verbosity",1);
addParameter(p,"ProgressIntervalSeconds",20);
addParameter(p,"QPOptions",nist_qp_options());

% Budgeted reduced-space primal heuristic.
addParameter(p,"EnableLocalSearch",true);
addParameter(p,"LocalSearchAtRoot",true);
addParameter(p,"LocalSearchStagnationNodes",500);
addParameter(p,"LocalSearchCooldownNodes",500);
addParameter(p,"LocalSearchPoorUBTrigger",true);
addParameter(p,"LocalSearchPoorUBRatio",0.10);
addParameter(p,"LocalSearchPoorUBMinNodes",250);
addParameter(p,"LocalSearchMaxCalls",12);
addParameter(p,"LocalSearchMaxFunctionEvaluations",80);
addParameter(p,"LocalSearchMaxIterations",30);
addParameter(p,"LocalSearchMaxTimePerCall",2.0);
addParameter(p,"LocalSearchMaxTotalTimeFraction",0.05);
addParameter(p,"LocalSearchOptimalityTolerance",1e-6);
addParameter(p,"LocalSearchStepTolerance",1e-8);
addParameter(p,"LocalSearchConstraintTolerance",1e-8);
parse(p,varargin{:});
options = p.Results;
options.BranchingRule = string(options.BranchingRule);
if ~ismember(options.BranchingRule,["normalized_longest_edge","absolute_longest_edge"])
    error("BranchingRule must be normalized_longest_edge or absolute_longest_edge.");
end
N = options.CoarseGridPoints; M = options.InitialApproxGridPoints;
if N<2 || M<N || mod(M-1,N-1)~=0
    error("InitialApproxGridPoints-1 must be a multiple of CoarseGridPoints-1.");
end
if options.MaximumFineGridPoints<N
    error("MaximumFineGridPoints must be at least CoarseGridPoints.");
end
if options.ProgressIntervalSeconds<15 && isfinite(options.ProgressIntervalSeconds)
    error("ProgressIntervalSeconds must be at least 15 seconds, or Inf to disable periodic progress output.");
end
if options.LocalSearchStagnationNodes<0 || options.LocalSearchCooldownNodes<0
    error("Local-search node intervals must be nonnegative.");
end
if options.LocalSearchPoorUBRatio<0 || options.LocalSearchPoorUBRatio>=1
    error("LocalSearchPoorUBRatio must satisfy 0 <= ratio < 1.");
end
if options.LocalSearchMaxCalls<0 || options.LocalSearchMaxFunctionEvaluations<1 || ...
        options.LocalSearchMaxIterations<1 || options.LocalSearchMaxTimePerCall<=0
    error("Invalid local-search budget.");
end
if options.LocalSearchMaxTotalTimeFraction<0 || options.LocalSearchMaxTotalTimeFraction>1
    error("LocalSearchMaxTotalTimeFraction must be in [0,1].");
end
end
