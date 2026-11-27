function row = nist_tolerance_result_row(result,problem,caseId,relative,absolute,timeLimit,errorMessage)
%NIST_TOLERANCE_RESULT_ROW Scaled objectives and unambiguous node counters.
arguments
    result (1,1) struct
    problem (1,1) struct
    caseId (1,1) string
    relative (1,1) double
    absolute (1,1) double
    timeLimit (1,1) double
    errorMessage (1,1) string = ""
end
row=nist_result_to_row(result,timeLimit,errorMessage);
row.problem_id=string(problem.id); row.instance_id=string(problem.instance_id);
row.n_linear=problem.nLinear; row.n_nonlinear=problem.nNonlinear;
row.n_residual=problem.nResidual;
row.nist_certified_rss=problem.certifiedRSS;
row.case_id=caseId; row.epsilon_relative=relative; row.epsilon_absolute=absolute;
row.objective_scale=problem.objectiveScale;
row.nist_reference_scaled=problem.certifiedRSS/problem.objectiveScale;
row.absolute_tolerance_rss=absolute*problem.objectiveScale;
row.solver_certified=row.certified;
row.pruning_threshold_scaled=(1-relative)*row.objective_scaled-absolute;
row.gap_margin_scaled=row.dual_bound_scaled-row.pruning_threshold_scaled;
row.bounds_consistent=isfinite(row.objective_scaled)&&isfinite(row.dual_bound_scaled)&& ...
    row.objective_scaled>=0&&row.dual_bound_scaled<=row.objective_scaled;
% Use the same configured test for both methods. Keep the solver's original
% flag separately: the legacy SCIP bridge permits a larger numerical margin.
row.certified=row.bounds_consistent&&row.gap_margin_scaled>=0;
row.par2=row.wall_time;
if ~row.certified||~isfinite(row.wall_time), row.par2=2*timeLimit; end
row.nodes_visited=row.nodes;
row.queue_extractions=NaN; row.nodes_created=NaN;
row.scip_nodes_all_runs=NaN;
row.relative_gap_scaled=NaN;
if isfinite(row.objective_scaled)&&isfinite(row.dual_bound_scaled)
    if row.objective_scaled>0
        row.relative_gap_scaled=max(0,(row.objective_scaled-row.dual_bound_scaled)/row.objective_scaled);
    elseif row.objective_scaled==0&&row.dual_bound_scaled==0
        row.relative_gap_scaled=0;
    end
end
names=["pruned_tolerance_only","pruned_at_incumbent","pruned_total", ...
    "pruned_tolerance_only_explicit","pruned_tolerance_only_global", ...
    "pruned_at_incumbent_explicit","pruned_at_incumbent_global", ...
    "nodes_pruned_explicit","nodes_closed_global"];
for name=names, row.(name)=NaN; end
row.prune_count_available=false;
row.max_node_depth=NaN;
row.m1_qp_retries=NaN; row.m2_qp_retries=NaN;
row.m1_unavailable_bounds=NaN; row.m2_unavailable_evaluations=NaN;
row.precision_fallback_m1_evaluations=NaN;
row.precision_limited_nodes=NaN;
stats=localField(result,'stats',struct());
if row.solver=="BnB"&&isfield(stats,'prunedToleranceOnly')
    row.nodes_visited=stats.nodesVisited;
    row.queue_extractions=stats.nodesProcessed;
    row.nodes_created=stats.nodesCreated;
    row.max_node_depth=stats.maxNodeDepth;
    row.m1_qp_retries=localField(stats,'m1QpRetries',0);
    row.m2_qp_retries=localField(stats,'m2QpRetries',0);
    row.m1_unavailable_bounds=localField(stats,'m1UnavailableBounds',0);
    row.m2_unavailable_evaluations=localField(stats,'m2UnavailableEvaluations',0);
    row.precision_fallback_m1_evaluations=localField(stats,'precisionFallbackM1Evaluations',0);
    row.precision_limited_nodes=localField(stats,'precisionLimitedNodes',0);
    fields=["prunedToleranceOnly","prunedAtIncumbent","prunedTotal", ...
        "prunedToleranceOnlyExplicit","prunedToleranceOnlyGlobal", ...
        "prunedAtIncumbentExplicit","prunedAtIncumbentGlobal", ...
        "nodesPruned","closedByGlobalGap"];
    for j=1:numel(names), row.(names(j))=stats.(fields(j)); end
    row.prune_count_available=true;
elseif row.solver=="SCIP"
    row.scip_nodes_all_runs=double(localField(result,'total_nodes',NaN));
end
row.tolerance_only_fraction_of_discarded=NaN;
if row.prune_count_available&&row.pruned_total>0
    row.tolerance_only_fraction_of_discarded=row.pruned_tolerance_only/row.pruned_total;
end
end

function value=localField(s,name,defaultValue)
if isstruct(s)&&isfield(s,name), value=s.(name); else, value=defaultValue; end
end
