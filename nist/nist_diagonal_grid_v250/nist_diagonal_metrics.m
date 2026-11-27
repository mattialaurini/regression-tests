function x=nist_diagonal_metrics(raw,st,options)
%NIST_DIAGONAL_METRICS Final gaps and unique-region cost classes.
U=raw.objective_scaled; L=raw.dual_bound_scaled;
consistent=isfinite(U)&&isfinite(L)&&U>=0&&L<=U;
absolute=NaN; relative=NaN;
if consistent
    absolute=U-L;
    if U>0, relative=absolute/U;
    elseif U==0&&L==0, relative=0;
    else, relative=Inf; end
end
x=struct('bounds_consistent',consistent,'gap_absolute_scaled',absolute, ...
    'gap_relative_scaled',relative,'gap_relative_percent',100*relative, ...
    'gap_absolute_rss',absolute*raw.objective_scale, ...
    'gap_margin_scaled',L-((1-raw.epsilon_relative)*U-raw.epsilon_absolute), ...
    'nodes_created',field(st,'nodesCreated',NaN), ...
    'queue_extractions',field(st,'nodesProcessed',NaN), ...
    'nodes_cheap',NaN,'nodes_expensive',NaN,'node_cost_counts_available',false, ...
    'node_cost_count_source',"unavailable", ...
    'screening_branches',field(st,'m1SkippedByScreening',NaN));
n=raw.nodes_visited;
if isfield(st,'nodesWithM1')&&isfield(st,'nodesWithoutM1')
    x.nodes_expensive=st.nodesWithM1; x.nodes_cheap=st.nodesWithoutM1;
    x.node_cost_counts_available=true; x.node_cost_count_source="instrumented_unique_regions";
elseif options.InitialApproxGridPoints==options.MaximumFineGridPoints && ...
        field(st,'nodesReinserted',NaN)==0 && ...
        field(st,'precisionFallbackM1Evaluations',NaN)==0 && ...
        isfield(st,'m1Evaluations') && st.m1Evaluations>=0 && st.m1Evaluations<=n
    % With M=Nhat and neither grid reinsertion nor precision fallback,
    % there is at most one m1 evaluation per visited region. This is
    % an exact reconstruction, not an assumption that evaluations=nodes.
    x.nodes_expensive=st.m1Evaluations; x.nodes_cheap=n-st.m1Evaluations;
    x.node_cost_counts_available=true;
    x.node_cost_count_source="derived_fixed_grid_no_fallback";
end
if x.node_cost_counts_available
    assert(x.nodes_cheap>=0 && x.nodes_expensive>=0 && ...
        x.nodes_cheap==floor(x.nodes_cheap) && x.nodes_expensive==floor(x.nodes_expensive) && ...
        x.nodes_cheap+x.nodes_expensive==n,'Invalid unique-node cost partition.');
end
end
function v=field(s,n,d)
if isstruct(s)&&isfield(s,n), v=s.(n); else, v=d; end
end
