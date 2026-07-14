function [Kmeans_results, rangeK] = LEiDA_cluster_VoxelMNI10mm(data_dir, file_V1, mink, maxk, replicates, results_dir, cluster_file)
% LEiDA_cluster_VoxelMNI10mm clusters leading eigenvectors into K groups 
% using the K-means algorithm.
%
% The algorithm is run for a range of clusters defined by mink and maxk. 
% For each value of K, the clustering is repeated a specified number of times 
% (replicates) to improve robustness.
%
% INPUT:
%   data_dir   - Directory containing the file with leading eigenvectors.
%   file_V1    - Filename of the .mat file containing variable V1_all (the leading eigenvectors).
%   mink       - Minimum number of clusters to explore (K_min).
%   maxk       - Maximum number of clusters to explore (K_max).
%   replicates - Number of independent replicates to run for each value of K.
%   results_dir- Directory where the results (cluster_file) will be saved.
%   cluster_file - Name of the file to save clustering outputs.
%
% OUTPUT:
%   Kmeans_results - Cell array where each cell contains the clustering solution 
%                    for a specific number of clusters. Each cell includes:
%                       .IDX - Cluster assignment for each sample (time point)
%                       .C   - Cluster centroids (functional connectivity patterns)
%   rangeK         - Vector containing the range of cluster numbers (K values) used.
%
% Authors: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

% Display a header message to indicate the start of clustering
disp('');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLUSTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

% -------------------------------------------------------------------------
% Load Leading Eigenvectors
% -------------------------------------------------------------------------
% Construct the full path to the file containing leading eigenvectors.
% The variable V1_all is assumed to be stored in the .mat file, possibly
% as 'single' to save disk space (e.g. Select_Demo_Subsample.m) - kmeans
% below expects double precision, so convert back if needed.
load([data_dir file_V1], 'V1_all');
if ~isa(V1_all, 'double')
    V1_all = double(V1_all);
end

% -------------------------------------------------------------------------
% Define Cluster Range and Prepare Result Storage
% -------------------------------------------------------------------------
% Create a vector of K values over which clustering will be performed.
rangeK = mink:maxk;

% Preallocate a cell array to store K-means clustering results for each K.
Kmeans_results = cell(size(rangeK));

% Display progress to the user.
disp(' ');
disp('Clustering eigenvectors into:');

% -------------------------------------------------------------------------
% Perform K-means Clustering for Each K
% -------------------------------------------------------------------------
for K = 1:length(rangeK)
    % Display the current number of clusters being used.
    disp(['   - ' num2str(rangeK(K)) ' FC states']);
    
    % Run K-means clustering on the eigenvector data (V1_all)
    % using cosine distance, a fixed number of replicates, and a maximum of 
    % 1000 iterations. The 'OnlinePhase' is turned off and parallel processing 
    % is enabled for speed.
    [IDX, C] = kmeans(V1_all, rangeK(K), 'Distance', 'Cosine', ...
                        'Replicates', replicates, 'MaxIter', 1000, ...
                        'OnlinePhase', 'off', 'Display', 'final', ...
                        'Options', statset('UseParallel', 1));
    
    % ---------------------------------------------------------------------
    % Sort Clusters by Occupancy
    % ---------------------------------------------------------------------
    % Calculate the histogram (count) for each cluster.
    % Sort clusters in descending order based on their occupancy (number of 
    % members in each cluster).
    [~, ind_sort] = sort(hist(IDX, 1:rangeK(K)), 'descend');
    
    % Rearrange the indexes so that clustering labels reflect the order.
    % This step ensures that cluster labels are consistently ordered.
    [~, idx_sort] = sort(ind_sort, 'ascend');
    
    % Save the sorted cluster assignments and centroids in the results cell array.
    % IDX is remapped to idx_sort so that the labels are ordered by occupancy.
    Kmeans_results{K}.IDX = idx_sort(IDX);
    % Reorder centroids to follow the sorted order of cluster occupancy.
    Kmeans_results{K}.C = C(ind_sort, :);
end

% -------------------------------------------------------------------------
% Load Additional Data and Save Results
% -------------------------------------------------------------------------
% Reload extra information stored in the same file; these may include:
%   - MNI_lowres_Mask: Low-resolution MNI brain mask.
%   - ind_voxels: Indices of voxels used for analysis.
%   - Time_sessions: Timing information of sessions.
%   - Data_info: Additional data parameters.
load([data_dir file_V1], 'MNI_lowres_Mask', 'ind_voxels', 'Scan_num', 'data_info');

% Save the clustering results, used parameters, and additional data to the specified results directory.
save([results_dir '/' cluster_file], 'Kmeans_results', 'rangeK', 'MNI_lowres_Mask', 'ind_voxels', 'Scan_num', 'data_info', 'replicates');

% Display a completion message with the name of the saved file.
disp(['K-means clustering completed and results saved as ' cluster_file]);
disp(' ');
