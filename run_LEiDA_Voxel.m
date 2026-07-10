%% run_LEiDA_Voxel.m

%==========================================================================
% README and Pipeline Setup for Functional Coupling Analysis using LEiDA
%==========================================================================
%
% OVERVIEW:
% Pipeline to analyze eigenvector dynamics in fMRI data using LEiDA, working
% directly in voxel space (a full-brain or custom voxel mask) rather than a
% pre-defined parcellation:
%   0. MASK (optional): Define the set of voxels of interest (e.g. full brain
%      or including CSF compartments) if not using the bundled full-brain mask.
%   1. EIGENVECTORS: Read the fMRI data and get the leading eigenvector of the
%      phase coherence matrix from each volume.
%   2. CLUSTER: Cluster all the leading eigenvectors extracted from fMRI scans
%      into a range of K clusters (coupling modes).
%   2b. HARMONIZE: Compute the fractional occupancy of each mode for every scan,
%      and harmonize it across acquisition sites using ComBat.
%   3. COMPARE CONDITIONS: Perform hypothesis testing between conditions
%      (permutation testing, with optional bootstrap and Hedge's effect sizes).
%   4. RESULTS: Automatically identify the modes that differ most between
%      conditions, and visualize their spatial patterns
%               - as slices with occupancy bar plots
%               - as 3D renders, including overlap with Yeo Resting-State Networks
%   5. SCORES: Correlate the key modes' occupancy with clinical/cognitive scores.
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
% 0. Mask_Voxels_of_Interest  (script, optional)
%    - Purpose: Builds a custom binary voxel mask in MNI space from a set of
%      per-scan brain masks (e.g. the CPAC preprocessing brain-extraction masks),
%      keeping voxels present in at least a chosen proportion of scans, and
%      resizing to the desired voxel size (e.g. 10mm3).
%    - This is a script: edit its "USER INPUT" variables (data_path,
%      extension_name, prop_masks, resize_scale) directly and run it.
%    - Not needed if using the bundled 'utilities/MNI_10mm3_FullBrain.mat' mask.
%
% 1. Get_EigenVectors_VoxelSpace_Server  (script)
%    - Purpose: Loads the preprocessed fMRI NIfTI files, resizes them to the
%      mask's voxel space, computes the signal phase (Hilbert transform), and
%      extracts the leading eigenvector of the phase coherence matrix at every
%      TR, for every scan.
%    - This is a script, not a function: it cannot be called with arguments.
%      Edit its "USER INPUT" section directly (fMRI_dir, extension_name,
%      results_dir, Mask_file, file_V1, TimeMax), then run it.
%    - Inputs (edited inside the script):
%      fMRI_dir      : Directory with the preprocessed fMRI NIfTI files.
%      extension_name: Common suffix of the fMRI filenames to load.
%      Mask_file     : .mat file with the voxel mask (variable MNI_lowres_Mask).
%      file_V1       : Output filename for the leading eigenvectors.
%      TimeMax       : Maximum number of volumes expected per scan.
%    - Output (saved to file_V1): V1_all, ind_voxels, MNI_lowres_Mask, data_info,
%      Scan_num, Scan_length.
%
% 2. LEiDA_cluster_VoxelMNI10mm
%    - Purpose: Clusters the leading eigenvectors into a range of K clusters.
%
%    Inputs:
%      data_dir    : Directory with the eigenvector file.
%      file_V1     : Name of the eigenvector file (output of step 1).
%      mink        : Minimum number of clusters (e.g., 2).
%      maxk        : Maximum number of clusters (e.g., 20).
%      replicates  : Number of replicates (e.g., 100).
%      results_dir : Directory where results will be saved.
%      cluster_file: Output file name for clustering results.
%
%    Output: [Kmeans_results, rangeK]
%
%    Example:
%      [Kmeans_results, rangeK] = LEiDA_cluster_VoxelMNI10mm(data_dir, file_V1, mink, maxk, replicates, results_dir, cluster_file);
%
% 2b. Harmonize mode occupancies with ComBat (inline in this script, uses combat/combat.m)
%    - Purpose: Computes the fractional occupancy of each mode for every scan
%      (for every K), then removes site effects from the occupancies with
%      ComBat (Johnson et al.), keeping diagnosis/age/sex/education as covariates
%      of interest. See the "2b. Extract Mode Occupancies..." section below.
%    - Requires a Scores table (Scores_Table) with per-scan SITE, AGE_AT_SCAN,
%      PTGENDER, PTEDUCAT and DX_num columns. This table is study-specific data
%      and is not included in this repository.
%
% 3. LEiDA_stats_Voxel_FracOccup_ComBat
%    - Purpose: Performs statistical tests (permutation, with optional bootstrap)
%      on the fractional occupancy of each mode between conditions. Uses Welch's
%      t-test for independent samples or a paired permutation test for paired
%      samples, and reports Hedge's effect sizes.
%
%    Inputs:
%      results_dir      : Directory where clustering results are stored.
%      file_cluster      : Clustering results file name.
%      file_stats        : Output file name for statistical results.
%      cond              : Cell array of condition labels (e.g., {'CN','MCI','DEM'}).
%      Index_Conditions   : Vector assigning each scan to a condition.
%      pair              : 0 for independent subjects, 1 for paired tests.
%      n_permutations    : Number of permutation samples (e.g., 1000).
%      n_bootstraps      : Number of bootstrap samples (e.g., 0-50).
%      P                 : Fractional occupancy matrix (raw or ComBat-harmonized).
%
%    Example:
%      LEiDA_stats_Voxel_FracOccup_ComBat(results_dir, cluster_file, stats_file, cond, Index_Conditions, pair, n_permutations, n_bootstraps, P);
%
% 4. Plot_FracOccup_stats
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
% 5. Plot_ClustVoxelCentroid_Pyramid_RSNs
%    - Purpose: Renders a pyramid of all centroids (across all K), each on a
%               transparent brain, with optional RSN overlay and significance markers.
%
%    Inputs:
%      results_dir  : Directory where results are stored.
%      file_clusters: Clustering results file name.
%      stats_file   : Statistics file name.
%      save_name    : Base name for saving the output figure.
%      overlap_Yeo  : Flag to use Yeo RSN colors (1=yes, 0=no).
%      cortex_dir   : View for rendering ('SideView' or 'TopView').
%      cond_pair    : Index for condition pair to report asterisks for, e.g.:
                            %1 : CN - MCI
                            %2 : CN - DEM
                            %3 : MCI - DEM
%      Add_asterisks: Flag to overlay significance markers (1=yes, 0=no).
%
%    Example:
%      Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, 'Centroid_Pyramid_Magenta_CN_AD', 0, 'SideView', 3, 1);
%
% 6. Choose_Relevant_Modes
%    - Purpose: Automatically selects the [k c] modes that differ most between
%               conditions (significant after multiple-testing correction, with
%               a minimum effect size), and groups correlated modes together.
%
%    Inputs:
%      results_dir : Directory where results are stored.
%      cluster_file: Clustering results file name.
%      stats_file  : Statistics file name.
%
%    Output: [Key_Modes_KC, Key_Centroids]
%
%    Example:
%      Key_Modes_KC = Choose_Relevant_Modes(results_dir, cluster_file, stats_file);
%
% 7. Plot_KeyModes_Slices_Stats
%    - Purpose: For each selected key mode, renders the centroid on anatomical
%               slices and plots mean +/- SE fractional occupancy bars per condition.
%
%    Inputs:
%      results_dir  : Directory where result files are stored.
%      cluster_file : Clustering results file name.
%      stats_file   : Statistical results file name.
%      save_name    : Base name for saving the output figure.
%      Key_Modes_KC : Nx2+ matrix with one row per key mode, [k c ...].
%
%    Example:
%      Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file, save_name, Key_Modes_KC);
%
% 8. Plot_Mode_TransparentBrain
%    - Purpose: Provides detailed 3D rendering for each selected coupling mode,
%               including its overlap with Yeo Resting-State Networks.
%
%    Inputs:
%      results_dir  : Directory with result files.
%      cluster_file : Clustering results file name.
%      Key_Modes_KC : Nx2+ matrix with one row per mode, [k c ...].
%
%    Example:
%      Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC);
%
% 9. Scores_vs_Mode_Occupancy
%    - Purpose: Correlates (partial correlation, controlling for age) the
%               occupancy of each key mode with a set of clinical/cognitive
%               scores, plots the results, and exports them to a CSV.
%    - NOTE: The set of score columns used (Genetics/Biomarkers/Cognitive_functions
%      indices) is hardcoded for the ADNI Scores_ADNI table used in Campo et al.;
%      adapt these indices for a different scores table.
%
%    Inputs:
%      Scores_Table : .mat file with the Scores_ADNI table.
%      Key_Modes_KC : Nx2+ matrix with one row per key mode, [k c ...].
%      results_dir  : Directory where results are stored.
%      stats_file   : Statistics file name (for the occupancy matrix P).
%      save_name    : Output .mat filename for the correlation results.
%
%    Example:
%      Scores_vs_Mode_Occupancy(Scores_Table, Key_Modes_KC, results_dir, stats_file, save_name);
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
%   1. Get the leading eigenvectors:
%     Edit and run Get_EigenVectors_VoxelSpace_Server.m (see step 1 above).
%
%   2. Cluster the eigenvectors:
%     mink = 2; maxk = 20; replicates = 100;
%     LEiDA_cluster_VoxelMNI10mm(data_dir, file_V1, mink, maxk, replicates, results_dir, cluster_file);
%
%   3. Harmonize occupancies and run the statistical analysis:
%     cond = {'CN','MCI','DEM'};
%     Paired_tests = 0; n_permutations = 1000; n_bootstraps = 0;
%     LEiDA_stats_Voxel_FracOccup_ComBat(results_dir, cluster_file, stats_file, cond, Index_Conditions, Paired_tests, n_permutations, n_bootstraps, P_harmonized);
%
%   4. Generate Figures:
%     a. Plot statistical results, including all p-values and effect sizes:
%        Plot_FracOccup_stats(results_dir, stats_file);
%
%     b. Centroid pyramid with RSN overlay:
%        overlap_RSNs = 0; cortex_dir = 'SideView'; Add_asterisks = 1; cond_pair = 3;
%        save_name = 'Centroid_Pyramid_Magenta_CN_AD';
%        Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, save_name, overlap_RSNs, cortex_dir, cond_pair, Add_asterisks);
%
%     c. Automatically select and plot the key modes that differ between conditions:
%        Key_Modes_KC = Choose_Relevant_Modes(results_dir, cluster_file, stats_file);
%        Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file, 'Fig1_Key_modes', Key_Modes_KC);
%
%     d. Detailed 3D visualization of the key modes, with RSN overlap:
%        Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC);
%
%     e. Correlate key mode occupancy with clinical/cognitive scores:
%        Scores_vs_Mode_Occupancy(Scores_Table, Key_Modes_KC, results_dir, stats_file, 'Scores_Mode_Stats.mat');
%
%==========================================================================
%% Setup Library of Directories and File Names 

% Main directory where the  fMRI data preprocessed in MNI space is stored as Nifti files (.nii or .nii.gz are stored)
fMRI_dir = '/Users/user/Documents/Research/CognitiveDecline/pipeline_cpac-adni-ants-mean-best-for-csf_0p01_0p1/';

% Part of the file name common to all fMRI scans that will be loaded
extension_name= 'rest_space-MNI152NLin6ASym_desc-preproc_bold.nii.gz';

%Directory where the LEiDA scripts are:
leida_dir = '/Users/user/Documents/Research/CognitiveDecline/LEiDA_V3';

% NOTE on the voxel mask: the binary mask in MNI space (1 in every voxel to be
% used as a 'region of interest', 0 elsewhere) is NOT set here. It is set as
% the Mask_file variable inside Get_EigenVectors_VoxelSpace_Server.m (step 1
% below), since that script cannot be called with arguments. The total number
% of voxels with 1 corresponds to the number of elements in the eigenvectors;
% at 10mm3 resolution the entire brain has ~1500-3000 voxels. See
% utilities/MNI_10mm3_FullBrain.mat for the bundled full-brain mask, or build
% a custom one with Mask_Voxels_of_Interest.m (step 0).

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
% Get_EigenVectors_VoxelSpace_Server.m is a script, not a function, so it
% cannot be called with arguments here. Instead, open that file, set its
% "USER INPUT" variables (fMRI_dir, extension_name, results_dir, Mask_file,
% file_V1, TimeMax) to match this section, and run it directly.
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