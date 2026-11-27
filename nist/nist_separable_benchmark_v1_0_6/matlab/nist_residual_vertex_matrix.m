function H = nist_residual_vertex_matrix(problem,z,vertices)
%NIST_RESIDUAL_VERTEX_MATRIX Residual vectors at all linear-domain vertices.
if nargin < 3 || isempty(vertices)
    vertices = problem.linearVertices;
end
Phi = nist_basis_matrix(problem,z);
H = Phi*double(vertices)-problem.response*ones(1,size(vertices,2));
end
