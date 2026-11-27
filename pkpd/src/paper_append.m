function R = paper_append(R, row)
%PAPER_APPEND Append a row (struct of scalars) to a struct of column vectors.
f = fieldnames(row);
for i = 1:numel(f)
    if ~isfield(R, f{i}), R.(f{i}) = zeros(0, 1); end
    R.(f{i})(end+1, 1) = row.(f{i});
end
end
