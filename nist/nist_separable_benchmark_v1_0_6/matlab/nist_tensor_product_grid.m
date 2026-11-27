function grid = nist_tensor_product_grid(lowerBound,upperBound,nPoints)
%NIST_TENSOR_PRODUCT_GRID Equally spaced Cartesian grid, one point per column.
lowerBound = double(lowerBound(:)); upperBound = double(upperBound(:));
if isscalar(nPoints), nPoints = repmat(nPoints,numel(lowerBound),1); end
nPoints = double(nPoints(:));
if numel(nPoints)~=numel(lowerBound) || any(nPoints<2) || any(nPoints~=floor(nPoints))
    error("Invalid grid specification.");
end
axes = cell(numel(lowerBound),1);
for k = 1:numel(lowerBound)
    axes{k} = linspace(lowerBound(k),upperBound(k),nPoints(k));
end
mesh = cell(size(axes));
[mesh{:}] = ndgrid(axes{:});
grid = zeros(numel(lowerBound),numel(mesh{1}));
for k = 1:numel(lowerBound)
    grid(k,:) = mesh{k}(:).';
end
end
