function n = paper_nrows(T)
%PAPER_NROWS Number of rows of a results table (struct of column vectors); 0 if empty.
if isstruct(T) && isfield(T, 'id'), n = numel(T.id); else, n = 0; end
end
