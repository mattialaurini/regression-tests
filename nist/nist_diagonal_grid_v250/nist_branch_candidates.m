function [candidates,midpoints] = nist_branch_candidates(region,problem,options)
%NIST_BRANCH_CANDIDATES Longest sides with a strictly interior midpoint.
% TieTolerance is relative to the CURRENT maximum length, without a unit floor.
lo=region.zLower(:); hi=region.zUpper(:);
lengths=hi-lo;
assert(all(isfinite(lo))&&all(isfinite(hi))&&all(lengths>=0), ...
    'NISTBranch:InvalidBox','Branching requires a finite, nonempty box.');
midpoints=lo+lengths/2;
% Avoid overflow in hi-lo for finite endpoints of opposite signs.
overflow=~isfinite(midpoints);
midpoints(overflow)=lo(overflow)/2+hi(overflow)/2;
bisectable=isfinite(midpoints)&lo<midpoints&midpoints<hi;
if ~any(bisectable), candidates=zeros(0,1); return, end
if string(options.BranchingRule)=="normalized_longest_edge"
    rootLengths=problem.zUpper(:)-problem.zLower(:);
    scale=max(rootLengths,eps(max(1,max(abs(rootLengths)))));
    lengths=lengths./scale;
end
scores=lengths; scores(~bisectable)=-Inf;
maximum=max(scores);
if maximum>0 && isfinite(maximum)
    relativeGap=(maximum-scores)/maximum;
    candidates=find(bisectable & relativeGap<=options.TieTolerance);
else
    % Underflow in normalization must not suppress representable progress.
    candidates=find(bisectable & scores==maximum);
end
assert(~isempty(candidates),'NISTBranch:NoCandidate','No longest bisectable side selected.');
end
