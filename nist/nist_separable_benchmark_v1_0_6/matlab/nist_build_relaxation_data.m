function data = nist_build_relaxation_data(problem,region,nPoints)
%NIST_BUILD_RELAXATION_DATA Build sampled residual columns and error bounds.
problem = validate_nist_problem(problem);
grid = nist_tensor_product_grid(region.zLower,region.zUpper,nPoints);
V = problem.linearVertices;
q = size(V,2); nGrid = size(grid,2); m = problem.nResidual;
G = zeros(m,q*nGrid);
for g = 1:nGrid
    columns = (g-1)*q+(1:q);
    G(:,columns) = nist_residual_vertex_matrix(problem,grid(:,g),V);
end
[cVertexMinus,cVertexPlus] = nist_interpolation_error_coefficients( ...
    problem,V,region.zLower,region.zUpper,nPoints);
data.G = G;
data.cMinus = repmat(cVertexMinus,1,nGrid);
data.cPlus = repmat(cVertexPlus,1,nGrid);
data.grid = grid;
data.nColumns = size(G,2);
data.nPoints = nPoints;
end
