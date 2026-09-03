function V1 = Get_LeadingEigenvector_RankTwo(theta)
% GET_LEADINGEIGENVECTOR_RANKTWO
%
% Computes the leading eigenvector of the instantaneous phase coherence
% matrix A(i,j) = cos(theta_i - theta_j) WITHOUT ever building the
% n_voxels x n_voxels matrix A.
%
% Rationale:
%   cos(theta_i - theta_j) = cos(theta_i)*cos(theta_j) + sin(theta_i)*sin(theta_j)
% so A = X*X' for X = [cos(theta), sin(theta)] (V x 2). A therefore has
% rank <= 2 regardless of how many voxels V there are. If w is the
% leading eigenvector of the 2x2 matrix X'*X (with eigenvalue lambda),
% then
%   A*(X*w) = X*(X'*X)*w = X*(lambda*w) = lambda*(X*w)
% so X*w, once normalized, IS the leading eigenvector of A. This
% collapses an O(V^2) eigendecomposition (EIGS on an n_voxels x n_voxels
% matrix, as in the original voxel-level implementation) into a 2x2
% eigendecomposition, which is what makes whole-brain voxel-resolution
% LEiDA (tens of thousands of voxels per TR) computationally tractable.
%
% This mirrors the rank-2 trick used in the custom out-of-core Python
% LEiDA pipeline (get_leading_eigenvectors), ported here for the MATLAB
% voxel-level pipeline.
%
% INPUT
%   theta   n_voxels x 1 (or 1 x n_voxels) vector of instantaneous voxel
%           phases (radians) at a single TR
%
% OUTPUT
%   V1      n_voxels x 1 leading eigenvector of cos(theta_i - theta_j),
%           unit norm. Sign is arbitrary (as with any eigenvector) and is
%           NOT fixed here; it is resolved once, vectorized across all
%           time points, by the existing sign convention in
%           Get_EigenVectors_VoxelSpace_Server.m. This keeps the function
%           a drop-in replacement for the previous eigs(...) call.

    theta = theta(:); % ensure column vector, n_voxels x 1

    X = [cos(theta), sin(theta)]; % n_voxels x 2
    M = X' * X;                   % 2 x 2 - cheap regardless of n_voxels

    [evecs, evals] = eig(M);      % eig returns eigenvalues in ascending order
    [~, idx] = max(diag(evals));
    w = evecs(:, idx);

    V1 = X * w;           % n_voxels x 1, eigenvector of the n_voxels x n_voxels matrix
    V1 = V1 / norm(V1);   % normalize to unit length, matching EIGS output
end
