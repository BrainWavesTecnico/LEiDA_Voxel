%% run_LEiDA_Vox.m

%==========================================================================
% README and Pipeline Setup for Functional Coupling Analysis using LEiDA 
%==========================================================================
%
% OVERVIEW:
% Pipeline to analyze eigenvector dynamics in fMRI data using LEiDA:
%   1. EIGENVECTORS: Read the fMRI data and get the leading eigenvectors from each volume. 
%   2. CLUSTER: Cluster all the leading eigenvectors extracted from fMRI scans into a range
%      of K clusters (coupling modes).
%   3. PATTERNS: Visualize the spatial pattern of each cluster centroid.
%   4. OCCUPANCY: Compute the fractional occupancy of each mode each coupling for every scan.
%   5.1 COMPARE CONDITIONS: Perform hypothesis testing between conditions, sample permutation and bootstraping.
%   5.2 COMPARE SCORES: Correlate with Scores, both continuous (using correlations) 
%      and binary scores (permutation testing) and correcting for multiple
%      testing
%   6. RESULTS: Visualize the spatial patterns of relevant modes
%               - compare with typical Resting-State Networks
%               - List Brain areas involved
%
% Figures are saved at each step in the results folder in both .fig and .png, 
%      showing the clustering results, statistical outcomes, and 
%      3D renderings of cluster centroids on brain space.
%
% Used in: Campo et al., Cognitive reserve linked to network-specific
% brain-ventricle coupling modes, 2025
%
%
% FUNCTIONS AND THEIR USAGE:
%
%  Get_EigenVectors_VoxelSpace
%
% 1. LEiDA_cluster_VoxelMNI10mm
%    - Purpose: Clusters the leading eigenvectors into a range of K clusters.
%    Note: the leading eigenvectors were previously extracted from the 599 fMRI scans realigned to MNI space
%    using function EigenVectors_VoxelSpace_v1.m
%
%    Inputs:
%      data_dir    : Directory with the eigenvector file.
%      file_V1     : Name of the eigenvector file (e.g., 'LEiDA_V1_all_MNI10mm.mat').
%      mink        : Minimum number of clusters (e.g., 2).
%      maxk        : Maximum number of clusters (e.g., 20).
%      replicates  : Number of replicates (e.g., 100).
%      results_dir : Directory where results will be saved.
%      cluster_file: Output file name for clustering results.
%
%    Example:
%      LEiDA_cluster_VoxelMNI10mm(data_dir, file_V1, mink, maxk, replicates, results_dir, cluster_file);
%
% 2. LEiDA_stats_Voxel_FracOccup
%    - Purpose: Computes the fractional occupancy of each coupling mode for each subject,
%               and performs statistical tests (permutation & bootstrap) between conditions.
%
%    Inputs:
%      results_dir  : Directory where clustering results are stored.
%      file_cluster : Clustering results file name.
%      stats_file   : Output file name for statistical results.
%      cond         : Cell array of condition labels (e.g., {'CN','EMCI','LMCI','AD'}).
%      pair         : 0 for independent subjects, 1 for paired tests.
%      n_permutations: Number of permutation samples (e.g., 1000 for quick tests).
%      n_bootstraps : Number of bootstrap samples (e.g., 2 for large samples).
%
%    Example:
%      LEiDA_stats_Voxel_FracOccup(results_dir, cluster_file, stats_file, Conditions_tag, Paired_tests, n_permutations, n_bootstraps);
%
% 3. Plot_KModes_TransparentBrain
%    - Purpose: Visualize the coupling modes (cluster centroids) and mean +/- standard error bars (fractional occupancy)
%               for a selected K on a transparent brain.
%
%    Inputs:
%      k            : Selected number of FC states to display (e.g., 5).
%      results_dir  : Directory with result files.
%      file_clusters: Clustering results file name.
%      file_stats   : Statistical results file name.
%
%    Example:
%      Plot_KModes_TransparentBrain(5, results_dir, cluster_file, stats_file);
%
% 4. Plot_ClustVoxelCentroid_Pyramid_RSNs
%    - Purpose: Renders a pyramid of centroids, each rendered on a transparent brain
%               with optional RSN overlay and significance markers.
%
%    Inputs:
%      results_dir : Directory where results are stored.
%      file_clusters: Clustering results file name.
%      stats_file  : Statistics file name.
%      save_name   : Base name for saving the output figure.
%      overlap_Yeo : Flag to use Yeo RSN colors (1=yes, 0=no).
%      cortex_dir  : View for rendering ('SideView' or 'TopView').
%      cond_pair   : Index for condition pair (e.g., 3 might indicate CN vs. AD).
%                    Choose which pairs of conditions compared for asterisks:
%                    Choose which pairs of conditions compared for asterisks:
                            %1 : CN - MCI
                            %2 : CN - DEM
                            %3 : MCI - DEM
%      Add_asterisks: Flag to overlay significance markers (1=yes, 0=no).
%
%    Example:
%      Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, 'Centroid_Pyramid_Magenta_CN_AD', 0, 'SideView', 3, 1);
%
% 5. Plot_FracOccup_stats
%    - Purpose: Generates plots summarizing the statistical tests on fractional 
%               occupancy (p-values, barplots of means, and Hedge's effect sizes).
%
%    Inputs:
%      results_dir : Directory where results are stored.
%      stats_file  : Statistics file name.
%
%    Example:
%      Plot_FracOccup_stats(results_dir, stats_file);
%
% 6. Plot_Modes_TransparentBrain
%    - Purpose: Provides detailed 3D rendering for an individual coupling mode
%               considering both the entire brain or only cortical voxels 
%               including its overlap with RSNs.
%
%    Inputs:
%      results_dir : Directory with result files.
%      file_clusters: Clustering results file name.
%      file_stats  : Statistics file name.
%      k           : Selected number of coupling modes (e.g., 4).
%      c           : Index of the mode to visualize (e.g., 2).
%
%    Example:
%      Plot_Modes_TransparentBrain(results_dir, cluster_file, stats_file, 4, 2);
%
% Complete Pipeline Example:
% ----------------------------------------------
%   Set the directories and file names:
%     data_dir    = '/path/to/data/';
%     results_dir = '/path/to/results/';
%     file_V1     = 'LEiDA_V1_all_MNI10mm.mat';
%     cluster_file= 'LEiDA_Clusters_VoxelMNI10mm.mat';
%     stats_file  = 'LEiDA_Stats.mat';
%
%   Add paths if necessary:
%     addpath(genpath(cd));
%
%   1. Cluster the eigenvectors:
%     mink = 2; maxk = 20; replicates = 100;
%     LEiDA_cluster_VoxelMNI10mm(data_dir, file_V1, mink, maxk, replicates, results_dir, cluster_file);
%
%   2. Statistical Analysis:
%     Conditions_tag = {'CN','EMCI','LMCI','AD'};
%     Paired_tests = 0; n_permutations = 1000; n_bootstraps = 2;
%     LEiDA_stats_Voxel_FracOccup(results_dir, cluster_file, stats_file, Conditions_tag, Paired_tests, n_permutations, n_bootstraps);
%
%   3. Generate Figures:
%     a. Repertoire of coupling modes for a selected K:
%        Plot_KModes_TransparentBrain(5, results_dir, cluster_file, stats_file);
%
%     b. Centroid pyramid with RSN overlay:
%        overlap_RSNs = 0; cortex_dir = 'SideView'; Add_asterisks = 1; cond_pair = 3;
%        save_name = 'Centroid_Pyramid_Magenta_CN_AD';
%        Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, save_name, overlap_RSNs, cortex_dir, cond_pair, Add_asterisks);
%
%     c. Plot statistical results, including all p-values and effect sizes
%        Plot_FracOccup_stats(results_dir, stats_file);
%
%     d. Detailed mode visualization:
%        Plot_Modes_TransparentBrain(results_dir, cluster_file, stats_file, 4, 2);
%
%     e. Alternatively, a centroid pyramid with Yeo RSN color overlay:
%        save_name = 'Centroid_Pyramid_RSNoverlap_CN_AD';
%        overlap_RSNs = 1; Add_asterisks = 0; cortex_dir = 'TopView';
%        Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, save_name, overlap_RSNs, cortex_dir, cond_pair, Add_asterisks);
%
%     
%
%==========================================================================
%% Setup Library of Directories and File Names 

% Main directory where the  fMRI data preprocessed in MNI space is stored as Nifti files (.nii or .nii.gz are stored)
fMRI_dir = '/Users/user/Documents/Research/CognitiveDecline/pipeline_cpac-adni-ants-mean-best-for-csf_0p01_0p1/';

% Part of the file name common to all fMRI scans that will be loaded
extension_name= 'rest_space-MNI152NLin6ASym_desc-preproc_bold.nii.gz';

%Directory where the LEiDA scripts are:
leida_dir = '/Users/user/Documents/Research/CognitiveDecline/LEiDA_V3';

% Name of the file containing a binary mask in MNI space, containing 1 in all the voxels 
% to be considered as 'regions of interest' or 'parcels' in the analysis. 
% The total number of voxels with 1 corresponds to the number of elements
% in the eigenvectors. With a 10mm3 size, the entire brain has 1500-3000
% voxels.

Mask_file='Mask_10mm_FullBrain.mat';
Scores_Table='Scores_ADNI_2177scans.mat';

% Directory where the results from this run will be stored
results_dir    = '/Users/user/Documents/Research/CognitiveDecline/FINAL_2177scans_CPAC/';
file_V1     = 'LEiDA_V1_all_MNI10mm_FullMask_2177scans.mat';    % File with leading eigenvectors
cluster_file= 'LEiDA_Clusters_VoxelMNI10mm_2177scans.mat';  % File to save clustering results
stats_file  = 'LEiDA_Stats_FracOccup_2177scans.mat';                  % File to save statistical results

% Add all subfolders of results directory to path
addpath(genpath(results_dir));
addpath(genpath(leida_dir));

%% 1. Get eigenvectors in voxel space
% 
% Get_EigenVectors_VoxelSpace(fMRI_dir,extension_name,results_dir,file_V1,Mask_file)
%
%% 2. Clustering the Eigenvectors into K Clusters
% 
% Define clustering parameters.
% mink = 2;
% maxk = 20;
% replicates = 1;
% % 
% Cluster the eigenvectors.
% LEiDA_cluster_VoxelMNI10mm(results_dir, file_V1, mink, maxk, replicates, results_dir, cluster_file);

%% 2b. Extract Mode Occupancies and correct for site differences

% Get Mode occupancies:
load(fullfile(results_dir, cluster_file));

n_scans = length(data_info);
unique_scans=unique(Scan_num);

% --- Extract features as the mode occupancies for all k and c ---
All_Occupancy = cell(length(rangeK), 1);
for ki = 1:length(rangeK)
    k = rangeK(ki);
    IDX = Kmeans_results{ki}.IDX;
    Occupancy = zeros(n_scans, k);
    for s = 1:n_scans
        Ctime = IDX(Scan_num == unique_scans(s));
        T  = length(Ctime);
        for c = 1:k
            Occupancy(s, c) = sum(Ctime == c) / T;
        end
    end
    All_Occupancy{ki} = Occupancy;
fprintf('Occupancies for k=%d modes in %d scans\n', k, size(All_Occupancy{ki}, 1));
end

% Harmonize Occupancies to account for site differences using COMBAT 

load(Scores_Table, 'Scores_ADNI');

% Covariates to remain unchanged 
age    = double(Scores_ADNI.AGE_AT_SCAN);
sex    = double(Scores_ADNI.PTGENDER == "Male");
edu    = double(Scores_ADNI.PTEDUCAT);
diagnose = Scores_ADNI.DX_num;
mod = [diagnose, age, sex, edu]; 

% Variable to harmonize
site   = double(Scores_ADNI.SITE)';

% --- ComBat ---
All_Occupancies_harmonized = cell(length(rangeK), 1);
for ki = 2:length(rangeK)
    fprintf('ComBat for k=%d\n', rangeK(ki))
    data_to_harmonize = All_Occupancy{ki}';   % k x n_scans
    All_Occupancies_harmonized{ki} = combat(data_to_harmonize, site, mod, 1)';

end

% --- Build harmonized P in the format for the rest of the analysis
P_harmonized = zeros(n_scans, length(rangeK), rangeK(end));
P_original = zeros(n_scans, length(rangeK), rangeK(end));

for ki = 1:length(rangeK)
    k = rangeK(ki);    
    if ki==1
        P_original(:,ki,1:ki)=All_Occupancy{ki};  
        P_harmonized(:,ki,1:ki)=All_Occupancy{ki}; 
    else
        P_original(:,ki,1:ki)=All_Occupancy{ki};  
        P_harmonized(:,ki,1:ki)=All_Occupancies_harmonized{ki}; 
    end
end

save(fullfile(results_dir, 'LEiDA_Occupancies_harmonized.mat'), 'P_original', 'P_harmonized', 'Scores_ADNI', 'rangeK', '-v7.3')

disp('Occupancies harmonized for site differences')

%% 3. Statistical Analysis of Mode Occupancies between conditions

% --- Definir condições ---
Index_Conditions = Scores_ADNI.DX_num+1; % (+ 1 so conditions are 1,2,3)

Condition_values=sort(unique(Index_Conditions));

disp('Number of scans in each condition:')
for cond = 1:length(Condition_values)
    Condition_tags{cond} = Scores_ADNI.DX{find(Scores_ADNI.DX_num+1==Condition_values(cond),1)};
    disp([num2str(cond) '- ' Condition_tags{cond} ' = ' num2str(sum(Index_Conditions==Condition_values(cond)))])
end

Paired_tests   = 0;
n_permutations = 100;
n_bootstraps   = 0;

LEiDA_stats_Voxel_FracOccup_ComBat(results_dir, cluster_file, stats_file, ...
    Condition_tags, Index_Conditions, Paired_tests, n_permutations, n_bootstraps, P_harmonized);

% Statistical Report of Fractional Occupancy across conditions
Plot_FracOccup_stats(results_dir, stats_file);
% Saves figures with main statistical results, including p-values, effect sizes and Barplot pyramid

%% Plot Centroid Pyramid with Statistical Results Overlay.
% Options:
%   - overlap_RSNs : 1 to color the voxels using Yeo RSN colors; 0 otherwise.
%   - cortex_dir : 'SideView' or 'TopView' to chose the camera angle of 3D brains
%   - Add_asterisks: 1 to add significance markers; 0 otherwise (SideView only)
%   - cond_pair    : Specify condition pair to show asterisks; e.g., 2 represents CN vs. AD.
overlap_RSNs = 1;
cortex_dir = 'SideView';
Add_asterisks = 0;

cond_pair = 2;
% Change name if you change the condition pair
save_name = 'Fig2_Centroid_Pyramid_RSNs_CN_vs_DEM';

Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, save_name, overlap_RSNs, cortex_dir, cond_pair, Add_asterisks);



%% 4. Figures with Results

% Select a number of Key Modes to analyze

% Automatically select the Key Modes differing mostly between conditions 
Key_Modes_KC = Choose_Relevant_Modes(results_dir, cluster_file, stats_file);

% OR

% Select all Modes for a given K
% k=8;
% Key_Modes_KC=[ones(1,k)*k;1:k]';

% OR

% Make your own selection in pairs [k c]
% Key_Modes_KC=[[3 4];[5 6];];

%% FIGURE 1. Plot the key modes differing between conditions 

save_name = 'Fig1_Key_modes_Slice_Occupancy_bars_';
Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file,save_name,Key_Modes_KC)

%% FIGURE 3. Plot detailed visualization of Key Modes and get the list of brain areas involved.

Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC);

%% Figure 6: Compare with scores

save_name='Scores_Mode_Stats.mat';
Scores_vs_Mode_Occupancy(Scores_Table,Key_Modes_KC,results_dir,stats_file,save_name)