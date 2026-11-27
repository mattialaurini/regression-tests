function result = solve_nist_bnb_diagonal_counted(problem,options)
%SOLVE_NIST_BNB_TOLERANCE_STUDY Reviewed v1.0.6 solver with pruning counters.
% Changes: unique visits; inclusive pruning; exact configured gap test;
% separate counts for explicit pruning and queue closure at global success.
% Branch fix: relative ties, strictly progressing splits, explicit precision limit.
%
% The linear coefficient box is represented by all its vertices and is kept
% fixed throughout the search. Only the nonlinear box is subdivided.
arguments
    problem (1,1) struct
    options (1,1) struct = nist_bnb_options()
end
problem = validate_nist_problem(problem);
options = localNormalizeOptions(options);
assert(options.InitialApproxGridPoints==options.MaximumFineGridPoints, ...
    "This instrumented entry point is restricted to M=Nhat.");
startTimer = tic;
lastProgressPrintTime = 0;
rootRegion = struct("zLower",problem.zLower,"zUpper",problem.zUpper);

stats = localInitialStats(problem.nNonlinear);
if isfield(problem,"initial") && ~isempty(problem.initial)
    incumbent.linear = problem.initial.linear;
    incumbent.z = problem.initial.z;
    [incumbent.objective,incumbent.residual,incumbent.prediction] = ...
        evaluate_nist_objective(problem,incumbent.linear,incumbent.z);
    incumbent.source = "initial_point";
else
    upper = solve_nist_linear_block_qp(problem,(problem.zLower+problem.zUpper)/2,options.QPOptions);
    if ~upper.hasSolution, error("Could not construct an initial feasible point."); end
    incumbent.linear = upper.linear; incumbent.z = (problem.zLower+problem.zUpper)/2;
    incumbent.objective = upper.objective; incumbent.residual = upper.residual;
    incumbent.prediction = upper.prediction; incumbent.source = "root_center";
end
stats.initialIncumbentObjective=incumbent.objective;
stats.lastIncumbentImprovementNode=0;
stats.lastIncumbentImprovementTime=0;

if options.EnableLocalSearch && options.LocalSearchAtRoot
    rootStart=(problem.zLower+problem.zUpper)/2;
    [incumbent,stats]=localApplyLocalSearch(problem,rootRegion,rootStart,incumbent, ...
        options,startTimer,stats,"root");
end

heap = localHeapCreate();
rootBranchDepth=zeros(1,problem.nNonlinear);
rootNode = localMakeNode(rootRegion,options.InitialApproxGridPoints,0,0,-Inf,rootBranchDepth);
heap = localHeapPush(heap,rootNode,-Inf); stats.nodesCreated = 1;
history = localHistoryCreate(options.RecordRelativeTargets,options.EpsilonRelative,options.EpsilonAbsolute);
history = localHistoryRecord(history,0,incumbent.objective,0,0);
status = "unknown"; lastGlobalLowerBound = 0;

while heap.count>0
    elapsed = toc(startTimer);
    globalLowerBound = min(localHeapMinimumValidLowerBound(heap),stats.minimumPrunedLowerBound);
    lastGlobalLowerBound = max(lastGlobalLowerBound,globalLowerBound);
    history = localHistoryRecord(history,elapsed,incumbent.objective,globalLowerBound,stats.nodesProcessed);
    if localCertificate(globalLowerBound,incumbent.objective,options), status="paper_gap"; break, end
    if elapsed>=options.TimeLimit, status="timelimit"; break, end
    if stats.nodesProcessed>=options.NodeLimit, status="nodelimit"; break, end

    [heap,node] = localHeapPop(heap); stats.nodesProcessed = stats.nodesProcessed+1;
    if ~node.visitedOnce
        stats.nodesVisited=stats.nodesVisited+1;
        node.visitedOnce=true;
    end
    stats.maxNodeDepth=max(stats.maxNodeDepth,node.depth);
    stats.maxProcessedBranchDepthByDimension=max(stats.maxProcessedBranchDepthByDimension,node.branchDepthByDimension);
    upperStart = tic;
    zCenter = (node.region.zLower+node.region.zUpper)/2;
    upper = solve_nist_linear_block_qp(problem,zCenter, ...
        localRemainingQPOptions(options.QPOptions,startTimer,options.TimeLimit));
    stats.upperBoundEvaluations = stats.upperBoundEvaluations+1;
    stats.upperBoundTime = stats.upperBoundTime+toc(upperStart);
    stats.qpTime = stats.qpTime+upper.wallTime;
    if upper.hasSolution && upper.objective<incumbent.objective
        incumbent.linear = upper.linear; incumbent.z = zCenter;
        incumbent.objective = upper.objective; incumbent.residual = upper.residual;
        incumbent.prediction = upper.prediction; incumbent.source = "node_center";
        stats.incumbentUpdates=stats.incumbentUpdates+1;
        stats.centerIncumbentUpdates=stats.centerIncumbentUpdates+1;
        stats.lastIncumbentImprovementNode=stats.nodesProcessed;
        stats.lastIncumbentImprovementTime=toc(startTimer);
    end
    if toc(startTimer)>=options.TimeLimit
        heap=localHeapPush(heap,node,node.priority); status="timelimit"; break
    end

    if node.hasCachedM2 && node.cachedM2M==node.M
        m2Value=node.cachedM2; node.hasCachedM2=false;
    else
        t=tic; m2=compute_nist_m2_screening_stable(problem,node.region, ...
            options.CoarseGridPoints,node.M, ...
            localRemainingQPOptions(options.QPOptions,startTimer,options.TimeLimit));
        m2Value=m2.value; stats.m2Evaluations=stats.m2Evaluations+1;
        stats.m2UnavailableEvaluations=stats.m2UnavailableEvaluations+double(~m2.screeningAvailable);
        stats.m2QpRetries=stats.m2QpRetries+m2.qpRetryCount;
        stats.m2Time=stats.m2Time+toc(t); stats.qpTime=stats.qpTime+m2.wallTime;
    end

    [runLocalSearch,localSearchReason,promiseScore]= ...
        localPeriodicLocalSearchTrigger(m2Value,incumbent.objective,stats,options,startTimer);
    stats.maxPoorUBPromiseScore=max(stats.maxPoorUBPromiseScore,promiseScore);
    if runLocalSearch
        [incumbent,stats]=localApplyLocalSearch(problem,node.region,zCenter,incumbent, ...
            options,startTimer,stats,localSearchReason);
    end
    if toc(startTimer)>=options.TimeLimit
        heap=localHeapPush(heap,node,node.priority); status="timelimit"; break
    end

    threshold=(1-options.EpsilonRelative)*incumbent.objective-options.EpsilonAbsolute;
    shouldBranch=false;
    if m2Value<threshold
        shouldBranch=true; stats.m1SkippedByScreening=stats.m1SkippedByScreening+1;
    elseif node.M<=options.MaximumFineGridPoints
        [node,stats]=localMarkM1(node,stats);
        t=tic; m1=compute_nist_m1_bound_stable(problem,node.region,node.M, ...
            localRemainingQPOptions(options.QPOptions,startTimer,options.TimeLimit));
        m1Value=m1.value; stats.m1Evaluations=stats.m1Evaluations+1;
        stats.m1UnavailableBounds=stats.m1UnavailableBounds+double(~m1.boundAvailable);
        stats.m1QpRetries=stats.m1QpRetries+m1.qpRetryCount;
        stats.m1Time=stats.m1Time+toc(t); stats.qpTime=stats.qpTime+m1.wallTime;
        node.validLowerBound=max(node.validLowerBound,m1Value);
        if m1Value<threshold
            shouldBranch=true;
        else
            stats.nodesPruned=stats.nodesPruned+1;
            [toleranceOnly,ordinary,eligible]=nist_tolerance_prune_masks( ...
                node.validLowerBound,incumbent.objective, ...
                options.EpsilonRelative,options.EpsilonAbsolute);
            assert(eligible,'Pruning classification inconsistent with the gap test.');
            stats.prunedToleranceOnlyExplicit=stats.prunedToleranceOnlyExplicit+toleranceOnly;
            stats.prunedAtIncumbentExplicit=stats.prunedAtIncumbentExplicit+ordinary;
            stats.minimumPrunedLowerBound=min(stats.minimumPrunedLowerBound,node.validLowerBound);
        end
    else
        qNew=max(1,floor((node.M-1)/(2*(options.CoarseGridPoints-1))));
        MNew=1+qNew*(options.CoarseGridPoints-1);
        if MNew>=node.M, error("Internal grid update did not reduce M."); end
        node.M=MNew; node.priority=m2Value; node.hasCachedM2=false;
        heap=localHeapPush(heap,node,node.priority); stats.nodesReinserted=stats.nodesReinserted+1;
    end

    if shouldBranch
        t=tic;
        [child1,child2,p1,p2,branchStats,selectedDimension]=localDivide(problem,node,options,startTimer);
        stats.branchingTime=stats.branchingTime+toc(t);
        stats.m2Evaluations=stats.m2Evaluations+branchStats.m2Evaluations;
        stats.m2UnavailableEvaluations=stats.m2UnavailableEvaluations+branchStats.m2UnavailableEvaluations;
        stats.m2QpRetries=stats.m2QpRetries+branchStats.m2QpRetries;
        stats.m2Time=stats.m2Time+branchStats.m2Time;
        stats.qpTime=stats.qpTime+branchStats.qpTime;
        if isempty(child1)
            % No representable bisection exists. Screening cannot justify
            % pruning: evaluate a valid m1 bound before deciding how to stop.
            if toc(startTimer)>=options.TimeLimit
                heap=localHeapPush(heap,node,node.priority); status="timelimit"; break
            end
            M=node.M;
            while M>options.MaximumFineGridPoints
                q=max(1,floor((M-1)/(2*(options.CoarseGridPoints-1))));
                nextM=1+q*(options.CoarseGridPoints-1);
                assert(nextM<M,'Invalid precision-fallback grid reduction.');
                M=nextM;
            end
            [node,stats]=localMarkM1(node,stats);
            t=tic;
            m1=compute_nist_m1_bound_stable(problem,node.region,M, ...
                localRemainingQPOptions(options.QPOptions,startTimer,options.TimeLimit));
            stats.m1Evaluations=stats.m1Evaluations+1;
            if m2Value<threshold
                stats.m1SkippedByScreening=stats.m1SkippedByScreening-1;
            end
            stats.precisionFallbackM1Evaluations=stats.precisionFallbackM1Evaluations+1;
            stats.m1UnavailableBounds=stats.m1UnavailableBounds+double(~m1.boundAvailable);
            stats.m1QpRetries=stats.m1QpRetries+m1.qpRetryCount;
            stats.m1Time=stats.m1Time+toc(t); stats.qpTime=stats.qpTime+m1.wallTime;
            node.M=M; node.validLowerBound=max(node.validLowerBound,m1.value);
            threshold=(1-options.EpsilonRelative)*incumbent.objective-options.EpsilonAbsolute;
            if node.validLowerBound>=threshold
                stats.nodesPruned=stats.nodesPruned+1;
                [toleranceOnly,ordinary,eligible]=nist_tolerance_prune_masks( ...
                    node.validLowerBound,incumbent.objective, ...
                    options.EpsilonRelative,options.EpsilonAbsolute);
                assert(eligible,'Precision-fallback pruning classification failed.');
                stats.prunedToleranceOnlyExplicit=stats.prunedToleranceOnlyExplicit+toleranceOnly;
                stats.prunedAtIncumbentExplicit=stats.prunedAtIncumbentExplicit+ordinary;
                stats.minimumPrunedLowerBound=min(stats.minimumPrunedLowerBound,node.validLowerBound);
            else
                % Preserve the unresolved region and its bound. Never turn a
                % precision limit into a global optimality certificate.
                heap=localHeapPush(heap,node,node.priority);
                stats.precisionLimitedNodes=stats.precisionLimitedNodes+1;
                status="numerical_resolution_limit";
                if toc(startTimer)>=options.TimeLimit, status="timelimit"; end
                break
            end
        else
            stats.nodesBranched=stats.nodesBranched+1; stats.nodesCreated=stats.nodesCreated+2;
            stats.branchCountByDimension(selectedDimension)=stats.branchCountByDimension(selectedDimension)+1;
            stats.maxBranchDepthByDimension=max(stats.maxBranchDepthByDimension,child1.branchDepthByDimension);
            stats.maxBranchDepthByDimension=max(stats.maxBranchDepthByDimension,child2.branchDepthByDimension);
            heap=localHeapPush(heap,child1,p1); heap=localHeapPush(heap,child2,p2);
        end
    end

    if options.Verbosity>=2 && isfinite(options.ProgressIntervalSeconds)
        progressElapsed = toc(startTimer);
        if progressElapsed-lastProgressPrintTime >= options.ProgressIntervalSeconds
            fprintf("NIST BnB %s: visited=%d, extracts=%d, queue=%d, UB=%.8g, valid-LB=%.8g, time=%.1fs\n", ...
                problem.id,stats.nodesVisited,stats.nodesProcessed,heap.count,incumbent.objective, ...
                localHeapMinimumValidLowerBound(heap),progressElapsed);
            lastProgressPrintTime = progressElapsed;
        end
    end
end

elapsed=toc(startTimer);
if heap.count==0 && status=="unknown"
    status="paper_gap";
    lastGlobalLowerBound=max(lastGlobalLowerBound,min(stats.minimumPrunedLowerBound,incumbent.objective));
elseif status=="unknown"
    status="stopped";
end
if heap.count>0
    lastGlobalLowerBound=max(lastGlobalLowerBound, ...
        min(localHeapMinimumValidLowerBound(heap),stats.minimumPrunedLowerBound));
end
lastGlobalLowerBound=min(lastGlobalLowerBound,incumbent.objective);
certificate=localCertificate(lastGlobalLowerBound,incumbent.objective,options);
if certificate && heap.count>0
    % These are existing queued regions closed by the final certificate,
    % not an estimate of descendants that would otherwise have been created.
    bounds=zeros(heap.count,1);
    for j=1:heap.count, bounds(j)=heap.nodes{j}.validLowerBound; end
    [toleranceOnly,ordinary,eligible]=nist_tolerance_prune_masks( ...
        bounds,incumbent.objective,options.EpsilonRelative,options.EpsilonAbsolute);
    assert(all(eligible),'A global certificate does not cover all queued regions.');
    stats.closedByGlobalGap=heap.count;
    stats.prunedToleranceOnlyGlobal=sum(toleranceOnly);
    stats.prunedAtIncumbentGlobal=sum(ordinary);
end
stats.nodesWithoutM1=stats.nodesVisited-stats.nodesWithM1;
assert(stats.nodesWithM1>=0 && stats.nodesWithoutM1>=0);
stats.prunedToleranceOnly=stats.prunedToleranceOnlyExplicit+stats.prunedToleranceOnlyGlobal;
stats.prunedAtIncumbent=stats.prunedAtIncumbentExplicit+stats.prunedAtIncumbentGlobal;
stats.prunedTotal=stats.nodesPruned+stats.closedByGlobalGap;
assert(stats.prunedTotal==stats.prunedToleranceOnly+stats.prunedAtIncumbent, ...
    'Pruning counters do not partition discarded nodes.');
history=localHistoryRecord(history,elapsed,incumbent.objective,lastGlobalLowerBound,stats.nodesProcessed);

[checkObjective,checkResidual,checkPrediction]=evaluate_nist_objective( ...
    problem,incumbent.linear,incumbent.z);
validationPassed=abs(checkObjective-incumbent.objective)<=1e-7*(1+abs(checkObjective)) && ...
    max(abs(checkResidual-incumbent.residual))<=1e-7*(1+max(abs(checkResidual)));
if ~validationPassed, error("NIST BnB result failed independent validation."); end

result=struct(); result.solver="BnB"; result.status=status; result.has_solution=true;
result.solution_source="bnb_incumbent"; result.incumbent_source=string(incumbent.source);
result.objective=incumbent.objective;
result.primal_bound=incumbent.objective; result.dual_bound=lastGlobalLowerBound;
result.paper_certificate_satisfied=certificate;
result.paper_certificate_margin=lastGlobalLowerBound- ...
    ((1-options.EpsilonRelative)*incumbent.objective-options.EpsilonAbsolute);
result.linear=incumbent.linear; result.z=incumbent.z;
result.residual=incumbent.residual; result.prediction=incumbent.prediction;
result.residual_original=incumbent.residual*problem.residualScale;
result.prediction_original=incumbent.prediction*problem.residualScale;
result.wall_time=elapsed; result.nodes=stats.nodesVisited; result.total_nodes=stats.nodesCreated;
result.queue_extractions=stats.nodesProcessed;
result.nodes_definition="distinct tree regions extracted at least once";
result.stats=stats; result.history=history; result.options=options;
result.branching_rule=string(options.BranchingRule);
result.problem_id=string(problem.id); result.instance_id=string(problem.instance_id);
result.n_linear=problem.nLinear; result.n_nonlinear=problem.nNonlinear;
result.n_residual=problem.nResidual; result.domain_profile=string(problem.domainProfile);
result.certified_rss=problem.certifiedRSS;
result.objective_scaled=result.objective;
result.dual_bound_scaled=result.dual_bound;
result.rss_original=result.objective*problem.objectiveScale;
result.dual_bound_original=result.dual_bound*problem.objectiveScale;
result.reference_relative_error=abs(result.rss_original-problem.certifiedRSS)/ ...
    max(problem.certifiedRSS,realmin);
result.objective_scale=problem.objectiveScale;
result.matlab_validation=struct("objective_recomputed",checkObjective, ...
    "objective_abs_difference",abs(checkObjective-incumbent.objective),"passed",validationPassed);
% Include final classification and independent validation in measured wall time.
result.wall_time=toc(startTimer);
end

function options=localNormalizeOptions(options)
defaults=nist_bnb_options(); names=fieldnames(defaults);
for k=1:numel(names), if ~isfield(options,names{k}), options.(names{k})=defaults.(names{k}); end, end
options=nist_bnb_options("TimeLimit",options.TimeLimit,"NodeLimit",options.NodeLimit, ...
    "EpsilonRelative",options.EpsilonRelative,"EpsilonAbsolute",options.EpsilonAbsolute, ...
    "CoarseGridPoints",options.CoarseGridPoints, ...
    "InitialApproxGridPoints",options.InitialApproxGridPoints, ...
    "MaximumFineGridPoints",options.MaximumFineGridPoints, ...
    "StrongBranchingTieBreak",options.StrongBranchingTieBreak, ...
    "BranchingRule",options.BranchingRule, ...
    "TieTolerance",options.TieTolerance,"RecordRelativeTargets",options.RecordRelativeTargets, ...
    "Verbosity",options.Verbosity,"ProgressIntervalSeconds",options.ProgressIntervalSeconds, ...
    "QPOptions",options.QPOptions, ...
    "EnableLocalSearch",options.EnableLocalSearch, ...
    "LocalSearchAtRoot",options.LocalSearchAtRoot, ...
    "LocalSearchStagnationNodes",options.LocalSearchStagnationNodes, ...
    "LocalSearchCooldownNodes",options.LocalSearchCooldownNodes, ...
    "LocalSearchPoorUBTrigger",options.LocalSearchPoorUBTrigger, ...
    "LocalSearchPoorUBRatio",options.LocalSearchPoorUBRatio, ...
    "LocalSearchPoorUBMinNodes",options.LocalSearchPoorUBMinNodes, ...
    "LocalSearchMaxCalls",options.LocalSearchMaxCalls, ...
    "LocalSearchMaxFunctionEvaluations",options.LocalSearchMaxFunctionEvaluations, ...
    "LocalSearchMaxIterations",options.LocalSearchMaxIterations, ...
    "LocalSearchMaxTimePerCall",options.LocalSearchMaxTimePerCall, ...
    "LocalSearchMaxTotalTimeFraction",options.LocalSearchMaxTotalTimeFraction, ...
    "LocalSearchOptimalityTolerance",options.LocalSearchOptimalityTolerance, ...
    "LocalSearchStepTolerance",options.LocalSearchStepTolerance, ...
    "LocalSearchConstraintTolerance",options.LocalSearchConstraintTolerance);
end

function stats=localInitialStats(nNonlinear)
stats=struct("nodesWithM1",0,"nodesWithoutM1",0,"nodesProcessed",0,"nodesVisited",0,"nodesCreated",0,"nodesBranched",0, ...
    "prunedToleranceOnlyExplicit",0,"prunedAtIncumbentExplicit",0, ...
    "closedByGlobalGap",0,"prunedToleranceOnlyGlobal",0, ...
    "prunedAtIncumbentGlobal",0,"prunedToleranceOnly",0, ...
    "prunedAtIncumbent",0,"prunedTotal",0, ...
    "nodesPruned",0,"nodesReinserted",0,"m1Evaluations",0,"m2Evaluations",0, ...
    "m1SkippedByScreening",0,"upperBoundEvaluations",0,"incumbentUpdates",0, ...
    "centerIncumbentUpdates",0,"localSearchIncumbentUpdates",0, ...
    "m1Time",0,"m2Time",0,"upperBoundTime",0,"branchingTime",0,"qpTime",0, ...
    "minimumPrunedLowerBound",Inf, ...
    "branchCountByDimension",zeros(1,nNonlinear), ...
    "maxBranchDepthByDimension",zeros(1,nNonlinear), ...
    "maxProcessedBranchDepthByDimension",zeros(1,nNonlinear), ...
    "maxNodeDepth",0,"precisionFallbackM1Evaluations",0,"precisionLimitedNodes",0, ...
    "m1UnavailableBounds",0,"m2UnavailableEvaluations",0,"m1QpRetries",0,"m2QpRetries",0, ...
    "localSearchCalls",0,"localSearchRootCalls",0,"localSearchPeriodicCalls",0, ...
    "localSearchStagnationTriggerCalls",0,"localSearchPoorUBTriggerCalls",0, ...
    "localSearchSuccessfulCalls",0,"localSearchFunctionEvaluations",0, ...
    "localSearchTime",0,"localSearchQPTime",0,"localSearchDisabled",false, ...
    "lastLocalSearchNode",-Inf,"lastIncumbentImprovementNode",0, ...
    "lastIncumbentImprovementTime",0,"initialIncumbentObjective",NaN, ...
    "bestLocalSearchAbsoluteImprovement",0,"bestLocalSearchRelativeImprovement",0, ...
    "maxPoorUBPromiseScore",0);
end

function [run,reason,promiseScore]=localPeriodicLocalSearchTrigger(m2Value,upper,stats,options,startTimer)
run=false; reason=""; promiseScore=0;
if ~options.EnableLocalSearch || stats.localSearchDisabled
    return
end
if stats.localSearchCalls>=options.LocalSearchMaxCalls
    return
end
if ~localLocalSearchBudgetAvailable(stats,options,startTimer)
    return
end
cooldown=stats.nodesProcessed-stats.lastLocalSearchNode>=options.LocalSearchCooldownNodes;
if ~cooldown
    return
end
stagnation=stats.nodesProcessed-stats.lastIncumbentImprovementNode>=options.LocalSearchStagnationNodes;
poorUB=false;
if options.LocalSearchPoorUBTrigger && stats.nodesProcessed>=options.LocalSearchPoorUBMinNodes && ...
        isfinite(upper) && upper>0 && isfinite(m2Value)
    promiseScore=max(0,min(1,1-m2Value/max(upper,realmin)));
    poorUB=m2Value<=options.LocalSearchPoorUBRatio*upper;
end
if stagnation || poorUB
    run=true;
    if stagnation && poorUB
        reason="stagnation+poor_ub";
    elseif stagnation
        reason="stagnation";
    else
        reason="poor_ub";
    end
end
end

function yes=localLocalSearchBudgetAvailable(stats,options,startTimer)
if stats.localSearchCalls>=options.LocalSearchMaxCalls
    yes=false; return
end
remainingGlobal=options.TimeLimit-toc(startTimer);
if isfinite(options.TimeLimit) && remainingGlobal<=0
    yes=false; return
end
if isfinite(options.TimeLimit)
    totalBudget=options.LocalSearchMaxTotalTimeFraction*options.TimeLimit;
else
    totalBudget=Inf;
end
yes=stats.localSearchTime<totalBudget;
end

function [incumbent,stats]=localApplyLocalSearch(problem,region,startZ,incumbent,options,startTimer,stats,reason)
if ~localLocalSearchBudgetAvailable(stats,options,startTimer)
    return
end
remainingGlobal=options.TimeLimit-toc(startTimer);
if isfinite(options.TimeLimit)
    totalBudget=options.LocalSearchMaxTotalTimeFraction*options.TimeLimit;
    remainingHeuristic=max(0,totalBudget-stats.localSearchTime);
else
    remainingHeuristic=Inf;
end
callLimit=min([options.LocalSearchMaxTimePerCall,remainingGlobal,remainingHeuristic]);
if ~isfinite(callLimit), callLimit=options.LocalSearchMaxTimePerCall; end
if callLimit<=1e-6
    return
end
lsOptions=struct("MaxIterations",options.LocalSearchMaxIterations, ...
    "MaxFunctionEvaluations",options.LocalSearchMaxFunctionEvaluations, ...
    "OptimalityTolerance",options.LocalSearchOptimalityTolerance, ...
    "StepTolerance",options.LocalSearchStepTolerance, ...
    "ConstraintTolerance",options.LocalSearchConstraintTolerance, ...
    "TimeLimit",callLimit, ...
    "QPOptions",localRemainingQPOptions(options.QPOptions,startTimer,options.TimeLimit));
oldObjective=incumbent.objective;
ls=solve_nist_reduced_local_search(problem,region,startZ,lsOptions);
stats.localSearchCalls=stats.localSearchCalls+1;
stats.lastLocalSearchNode=stats.nodesProcessed;
stats.localSearchTime=stats.localSearchTime+ls.wallTime;
stats.localSearchQPTime=stats.localSearchQPTime+ls.qpTime;
stats.qpTime=stats.qpTime+ls.qpTime;
stats.localSearchFunctionEvaluations=stats.localSearchFunctionEvaluations+ls.functionEvaluations;
if reason=="root"
    stats.localSearchRootCalls=stats.localSearchRootCalls+1;
else
    stats.localSearchPeriodicCalls=stats.localSearchPeriodicCalls+1;
    if contains(reason,"stagnation")
        stats.localSearchStagnationTriggerCalls=stats.localSearchStagnationTriggerCalls+1;
    end
    if contains(reason,"poor_ub")
        stats.localSearchPoorUBTriggerCalls=stats.localSearchPoorUBTriggerCalls+1;
    end
end
if ~ls.available
    stats.localSearchDisabled=true;
    return
end
if ls.hasSolution
    stats.localSearchSuccessfulCalls=stats.localSearchSuccessfulCalls+1;
end
if ls.hasSolution && ls.objective<incumbent.objective
    incumbent.linear=ls.linear; incumbent.z=ls.z;
    incumbent.objective=ls.objective; incumbent.residual=ls.residual;
    incumbent.prediction=ls.prediction;
    if reason=="root", incumbent.source="root_local_search"; else, incumbent.source="periodic_local_search"; end
    stats.incumbentUpdates=stats.incumbentUpdates+1;
    stats.localSearchIncumbentUpdates=stats.localSearchIncumbentUpdates+1;
    stats.lastIncumbentImprovementNode=stats.nodesProcessed;
    stats.lastIncumbentImprovementTime=toc(startTimer);
    absoluteImprovement=oldObjective-ls.objective;
    relativeImprovement=absoluteImprovement/max(abs(oldObjective),realmin);
    stats.bestLocalSearchAbsoluteImprovement=max(stats.bestLocalSearchAbsoluteImprovement,absoluteImprovement);
    stats.bestLocalSearchRelativeImprovement=max(stats.bestLocalSearchRelativeImprovement,relativeImprovement);
end
end

function node=localMakeNode(region,M,validLowerBound,depth,priority,branchDepthByDimension)
node=struct("region",region,"M",M,"validLowerBound",validLowerBound, ...
    "visitedOnce",false,"didRunM1",false, ...
    "depth",depth,"branchDepthByDimension",branchDepthByDimension, ...
    "priority",priority,"hasCachedM2",false, ...
    "cachedM2",NaN,"cachedM2M",NaN);
end

function yes=localCertificate(lower,upper,options)
yes=isfinite(upper) && lower>=(1-options.EpsilonRelative)*upper-options.EpsilonAbsolute;
end

function [child1,child2,priority1,priority2,stats,selected]=localDivide(problem,parent,options,startTimer)
region=parent.region;
candidates=nist_branch_candidates(region,problem,options);
stats=struct("m2Evaluations",0,"m2Time",0,"qpTime",0,"m2UnavailableEvaluations",0,"m2QpRetries",0);
if isempty(candidates)
    child1=[]; child2=[]; priority1=0; priority2=0; selected=0;
    return
end
if numel(candidates)==1 || ~options.StrongBranchingTieBreak
    selected=candidates(1);
    [region1,region2]=localSplitRegion(region,selected);
    [priority1,qp1,time1,unavailable1,retries1]=localChildPriority(problem,region1,options,startTimer);
    [priority2,qp2,time2,unavailable2,retries2]=localChildPriority(problem,region2,options,startTimer);
    stats.m2QpRetries=retries1+retries2;
    stats.m2UnavailableEvaluations=unavailable1+unavailable2;
    stats.m2Evaluations=2; stats.m2Time=time1+time2; stats.qpTime=qp1+qp2;
else
    best=-Inf; region1=[]; region2=[]; priority1=0; priority2=0;
    for j=candidates(:).'
        [r1,r2]=localSplitRegion(region,j);
        [p1,qp1,time1,unavailable1,retries1]=localChildPriority(problem,r1,options,startTimer);
        [p2,qp2,time2,unavailable2,retries2]=localChildPriority(problem,r2,options,startTimer);
        stats.m2QpRetries=stats.m2QpRetries+retries1+retries2;
        stats.m2UnavailableEvaluations=stats.m2UnavailableEvaluations+unavailable1+unavailable2;
        stats.m2Evaluations=stats.m2Evaluations+2;
        stats.m2Time=stats.m2Time+time1+time2; stats.qpTime=stats.qpTime+qp1+qp2;
        score=min(p1,p2);
        if score>best, best=score; selected=j; region1=r1; region2=r2; priority1=p1; priority2=p2; end
    end
end
childBranchDepth=parent.branchDepthByDimension;
childBranchDepth(selected)=childBranchDepth(selected)+1;
child1=localMakeNode(region1,options.InitialApproxGridPoints,parent.validLowerBound,parent.depth+1,priority1,childBranchDepth);
child2=localMakeNode(region2,options.InitialApproxGridPoints,parent.validLowerBound,parent.depth+1,priority2,childBranchDepth);
child1.hasCachedM2=true; child1.cachedM2=priority1; child1.cachedM2M=options.InitialApproxGridPoints;
child2.hasCachedM2=true; child2.cachedM2=priority2; child2.cachedM2M=options.InitialApproxGridPoints;
end

function [r1,r2]=localSplitRegion(region,j)
r1=region; r2=region;
lo=region.zLower(j); hi=region.zUpper(j);
midpoint=lo+(hi-lo)/2;
if ~isfinite(midpoint), midpoint=lo/2+hi/2; end
assert(lo<midpoint && midpoint<hi,'NISTBranch:NoProgress', ...
    'The selected side has no strictly interior floating-point midpoint.');
r1.zUpper(j)=midpoint; r2.zLower(j)=midpoint;
end

function [priority,qpTime,totalTime,unavailable,retries]=localChildPriority(problem,region,options,startTimer)
t=tic; unavailable=0; retries=0;
if toc(startTimer)>=options.TimeLimit, priority=0; qpTime=0; totalTime=0; return, end
bound=compute_nist_m2_screening_stable(problem,region,options.CoarseGridPoints, ...
    options.InitialApproxGridPoints,localRemainingQPOptions(options.QPOptions,startTimer,options.TimeLimit));
priority=bound.value; qpTime=bound.wallTime; totalTime=toc(t);
unavailable=double(~bound.screeningAvailable); retries=bound.qpRetryCount;
end

function history=localHistoryCreate(targets,finalTarget,epsilonAbsolute)
targets=unique([double(targets(:).'),double(finalTarget)]);
history=struct("time",zeros(0,1),"upper",zeros(0,1),"lower",zeros(0,1), ...
    "nodes",zeros(0,1),"targets",targets,"targetTimes",nan(size(targets)), ...
    "targetNodes",nan(size(targets)),"epsilonAbsolute",epsilonAbsolute);
end

function history=localHistoryRecord(history,time,upper,lower,nodes)
history.time(end+1,1)=time; history.upper(end+1,1)=upper;
history.lower(end+1,1)=lower; history.nodes(end+1,1)=nodes;
for k=1:numel(history.targets)
    if isnan(history.targetTimes(k)) && lower>=(1-history.targets(k))*upper-history.epsilonAbsolute
        history.targetTimes(k)=time; history.targetNodes(k)=nodes;
    end
end
end

function qpOptions=localRemainingQPOptions(baseOptions,startTimer,timeLimit)
qpOptions=baseOptions; remaining=timeLimit-toc(startTimer);
if isfinite(timeLimit), remaining=max(0.001,remaining); else, remaining=Inf; end
if isfield(baseOptions,"TimeLimit") && isfinite(baseOptions.TimeLimit)
    remaining=min(remaining,double(baseOptions.TimeLimit));
end
qpOptions.TimeLimit=remaining;
end

function heap=localHeapCreate()
heap=struct("nodes",{{}},"priorities",zeros(0,1),"count",0);
end
function heap=localHeapPush(heap,node,priority)
heap.count=heap.count+1; i=heap.count; heap.nodes{i,1}=node; heap.priorities(i,1)=priority;
while i>1
    parent=floor(i/2); if heap.priorities(parent)<=heap.priorities(i), break, end
    [heap.nodes{parent},heap.nodes{i}]=deal(heap.nodes{i},heap.nodes{parent});
    [heap.priorities(parent),heap.priorities(i)]=deal(heap.priorities(i),heap.priorities(parent)); i=parent;
end
end
function [heap,node]=localHeapPop(heap)
if heap.count==0, error("Cannot pop an empty heap."); end
node=heap.nodes{1};
if heap.count==1, heap.nodes=cell(0,1); heap.priorities=zeros(0,1); heap.count=0; return, end
heap.nodes{1}=heap.nodes{heap.count}; heap.priorities(1)=heap.priorities(heap.count);
heap.nodes(heap.count)=[]; heap.priorities(heap.count)=[]; heap.count=heap.count-1; i=1;
while true
    left=2*i; right=left+1; if left>heap.count, break, end
    smallest=left; if right<=heap.count && heap.priorities(right)<heap.priorities(left), smallest=right; end
    if heap.priorities(i)<=heap.priorities(smallest), break, end
    [heap.nodes{i},heap.nodes{smallest}]=deal(heap.nodes{smallest},heap.nodes{i});
    [heap.priorities(i),heap.priorities(smallest)]=deal(heap.priorities(smallest),heap.priorities(i)); i=smallest;
end
end
function lower=localHeapMinimumValidLowerBound(heap)
if heap.count==0, lower=Inf; return, end
values=zeros(heap.count,1); for k=1:heap.count, values(k)=heap.nodes{k}.validLowerBound; end
lower=min(values);
end

function [node,stats]=localMarkM1(node,stats)
if ~node.didRunM1
    node.didRunM1=true;
    stats.nodesWithM1=stats.nodesWithM1+1;
end
end
