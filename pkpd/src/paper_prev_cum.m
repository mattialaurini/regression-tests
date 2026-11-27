function c = paper_prev_cum(T1, k, ids, current, nids, tmo, allowed)
%PAPER_PREV_CUM Cumulative PAR2 on the patients ids of the settings of order k that
%   were completed (all nids tuning patients) in earlier rounds, i.e. not in
%   current, restricted to the settings in allowed (optional; e.g. the settings of the
%   same group); [] if there is none. Used by the racing of PAPER_TESTS_LIN.
c = zeros(0, 1);
if ~isfield(T1, 'order'), return; end
S = setdiff(unique(T1.setting(T1.order == k)), current);
if nargin >= 7, S = intersect(S, allowed); end
for s = S(:)'
    sel = T1.order == k & T1.setting == s;
    if sum(sel) < nids, continue; end
    t = T1.time(sel); t(~T1.converged(sel)) = 2 * tmo; idv = T1.id(sel);
    c(end+1, 1) = sum(t(ismember(idv, ids)));
end
end
