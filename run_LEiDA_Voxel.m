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
%      Confirm RAM availibility creating a matrix N_voxels x (T_scanxN_scans).
%   1. EIGENVECTORS: Read the fMRI data realigned to a common template,
%      Resize fMRI volumes to the Mask size, extract the time series from the N voxels of interest,
%      Compute the Hilbert transform and extract the signal phase,
%      get the 1xN leading eigenvector of the phase coherence matrix,
%      concatenate the leading eigenvetors from all the volumes and all scans, with size N_voxels x(TxN_scans)
%   2. CLUSTER: Cluster all the leading eigenvectors extracted from fMRI scans
%      into a range of K clusters (coupling modes).
%   2b. HARMONIZE: Compute the fractional occupancy of each mode for every scan,
%      and, optionally, harmonize it across acquisition sites using ComBat.
%      Both the original and harmonized occupancies are saved, and downstream
%      steps can be run on either one.
%   3a. COMPARE CONDITIONS (optional): Perform hypothesis testing between
%      conditions (permutation testing, with optional bootstrap and Hedge's
%      effect sizes), on either the original or the harmonized occupancies.
%      Skip this step for studies with no discrete conditions to compare.
%   3b. COMPARE SCORES: Correlate every clinical/cognitive score against every
%      mode in the entire pyramid (all K), report/save the significant
%      results. Only needs the occupancy matrix from step 2b - does not
%      depend on step 3a, or on any pre-selected set of key modes.
%   4. RESULTS: Identify a set of key modes - automatically, from step 3a's
%      significant/large-effect modes, or manually - and visualize their
%      spatial patterns
%               - as slices with occupancy bar plots
%               - as 3D renders, including overlap with Yeo Resting-State Networks
%               - against the scores from step 3b
%
% Figures are saved at each step in the results folder in both .fig and .png,
%      showing the clustering results, statistical outcomes, and
%      3D renderings of cluster centroids on brain space.
%
% Used in: Cognitive function linked to temporal occupancy of Brain-Ventricle
% (BraVe) modes. Campo, Miguel, Brattico, Nigro, Tafuri, Logroscino, Cabral
% and the Alzheimer's Disease Neuroimaging Initiative (ADNI).
% bioRxiv 2025.01.04.631289; doi: https://doi.org/10.1101/2025.01.04.631289
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
% 2b. Save_Occupancies_Harmonize (uses combat/combat.m)
%    - Purpose: Computes the fractional occupancy of each mode for every scan
%      (for every K), and, if apply_combat is 1, removes site effects from the
%      occupancies with ComBat (Johnson et al.), keeping diagnosis/age/sex/
%      education as covariates of interest. Saves both P_original and
%      P_harmonized (P_harmonized = P_original when apply_combat is 0), so a
%      later step can choose which one to run statistics on.
%    - Requires a Scores table (Scores_Table) with per-scan SITE, AGE_AT_SCAN,
%      PTGENDER, PTEDUCAT and DX_num columns. This table is study-specific data
%      and is not included in this repository.
%
%    Inputs:
%      results_dir  : Directory with the cluster file.
%      cluster_file : Clustering results file name (output of step 2).
%      Scores_Table : .mat file with the Scores_ADNI table.
%      apply_combat : 1 to harmonize across sites with ComBat; 0 to skip.
%      occup_file   : Output file name for the saved occupancies.
%
%    Output: [P_original, P_harmonized, rangeK, Scores_ADNI]
%
%    Example:
%      [P_original, P_harmonized, rangeK, Scores_ADNI] = Save_Occupancies_Harmonize(results_dir, cluster_file, Scores_Table, apply_combat, occup_file);
%
% 3a. LEiDA_stats_Voxel_FracOccup_ComBat
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
%      P                 : Fractional occupancy matrix; pass P_original or
%                          P_harmonized from step 2b, whichever you want to
%                          test. Pass the same P to Scores_vs_Mode_Occupancy
%                          (step 3b) to keep the two analyses consistent.
%
%    Example:
%      LEiDA_stats_Voxel_FracOccup_ComBat(results_dir, cluster_file, stats_file, cond, Index_Conditions, pair, n_permutations, n_bootstraps, P);
%
% 3b. Scores_vs_Mode_Occupancy
%    - Purpose: Correlates every clinical/cognitive score against every mode
%      in the entire pyramid (every K in rangeK), partial correlation
%      controlling for age. Does NOT depend on any pre-selected set of modes
%      and does NOT generate a figure - it only computes and saves the
%      p-values, and reports the most significant results to the command
%      line (any mode surviving Bonferroni across the whole pyramid AND all
%      scores tested: 0.05/sum(rangeK)/N_scores). Use Plot_KeyModes_vs_Scores
%      (step 9 below) afterwards to plot a selected set of modes against
%      these scores. The saved p-values can also be loaded by
%      Plot_ClustVoxelCentroid_Pyramid_RSNs.m (step 5) in place of the
%      permutation-test p-values (pass pyramid_stats_file as its stats_file;
%      stat_of_interest then selects which score to display).
%    - Takes the occupancy matrix P directly (e.g. P_original or P_harmonized
%      from step 2b), NOT from stats_file, so it does not depend on step 3a
%      (LEiDA_stats_Voxel_FracOccup_ComBat) having been run. This makes it
%      usable on its own for studies with no discrete conditions to compare,
%      only continuous scores.
%    - NOTE: The set of score columns used (Genetics/Biomarkers/Cognitive_functions
%      indices) is hardcoded for the ADNI Scores_ADNI table used in Campo et al.;
%      adapt these indices for a different scores table.
%
%    Inputs:
%      P                  : Fractional occupancy matrix (P_original or P_harmonized).
%      Scores_Table       : .mat file with the Scores_ADNI table.
%      results_dir        : Directory where the mat output is saved.
%      cluster_file       : Clustering results file name (used to load rangeK).
%      pyramid_stats_file : Output .mat filename for the entire-pyramid p-values.
%
%    Example:
%      Scores_vs_Mode_Occupancy(P, Scores_Table, results_dir, cluster_file, pyramid_stats_file);
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
%               transparent brain, with optional RSN overlay. Add_asterisks
%               overlays a significance marker on top of each cluster
%               centroid, showing the chosen statistical result
%               (stat_of_interest) for that mode.
%
%    Inputs:
%      results_dir      : Directory where results are stored.
%      file_clusters     : Clustering results file name.
%      stats_file        : Statistics file name - either the output of step
%                          3a (P_pval indexed by condition pair) or the
%                          pyramid_stats_file from step 3b (P_pval indexed by
%                          score instead).
%      save_name         : Base name for saving the output figure. The name
%                          of the selected comparison is appended
%                          automatically (see stat_of_interest below).
%      overlap_Yeo       : Flag to use Yeo RSN colors (1=yes, 0=no).
%      cortex_dir        : View for rendering ('SideView' or 'TopView').
%      stat_of_interest  : Index into P_pval's first dimension, selecting
%                          WHICH statistical result is displayed on top of
%                          each centroid: a condition pair (e.g. CN-AD) when
%                          stats_file comes from step 3a, e.g.:
                            %1 : CN - MCI
                            %2 : CN - DEM
                            %3 : MCI - DEM
%                          or a score (e.g. MoCA) when stats_file comes from
%                          step 3b.
%      Add_asterisks     : Flag to overlay significance markers (1=yes, 0=no).
%
%    Example:
%      Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, 'Centroid_Pyramid_Magenta', 0, 'SideView', 3, 1);
%
% 6. Choose_Relevant_Modes
%    - Purpose: Automatically selects the [ki c] modes that differ most between
%               conditions (significant after multiple-testing correction, with
%               a minimum effect size), and groups correlated modes together.
%
%    Inputs:
%      results_dir : Directory where results are stored.
%      cluster_file: Clustering results file name.
%      stats_file  : Statistics file name.
%
%    Output: [Key_Modes_KC, Key_Centroids]
%      Key_Modes_KC is an Nx12 matrix, one row per key mode: [ki, c, slope, ...].
%      IMPORTANT: ki is the POSITION of the clustering solution in rangeK
%      (i.e. Kmeans_results{ki}, P(:,ki,:)), NOT the literal number of
%      clusters - e.g. if rangeK = 2:20, ki=1 means K=2 clusters, ki=2 means
%      K=3, etc. If you build Key_Modes_KC manually instead of via this
%      function (see step 4 below), use the same ki convention.
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
%      Key_Modes_KC : Nx3+ matrix, [ki c slope ...] - see step 6 above for the
%                     ki (position, not literal K) convention.
%      Scores_Table : .mat file with the Scores_ADNI table (used to split by sex).
%
%    Example:
%      Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file, save_name, Key_Modes_KC, Scores_Table);
%
% 8. Plot_Mode_TransparentBrain
%    - Purpose: Provides detailed 3D rendering for each selected coupling mode,
%               including its overlap with Yeo Resting-State Networks.
%
%    Inputs:
%      results_dir  : Directory with result files.
%      cluster_file : Clustering results file name.
%      Key_Modes_KC : Nx2+ matrix, [ki c ...] - see step 6 above for the ki
%                     (position, not literal K) convention.
%
%    Example:
%      Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC);
%
% 9. Plot_KeyModes_vs_Scores
%    - Purpose: Figure-generation counterpart to Scores_vs_Mode_Occupancy
%      (step 3b): plots one bar (of correlation coefficients across scores)
%      per selected key mode, and exports a CSV. Unlike step 3b, this DOES
%      depend on a pre-selected set of modes (Key_Modes_KC).
%    - NOTE: The set of score columns used (Genetics/Biomarkers/Cognitive_functions
%      indices) is hardcoded for the ADNI Scores_ADNI table used in Campo et al.;
%      adapt these indices for a different scores table.
%
%    Inputs:
%      P            : Fractional occupancy matrix (P_original or P_harmonized).
%      Scores_Table : .mat file with the Scores_ADNI table.
%      Key_Modes_KC : Nx2+ matrix, [ki c ...] - see step 6 above for the ki
%                     (position, not literal K) convention.
%      results_dir  : Directory where the figure/CSV/mat outputs are saved.
%      save_name    : Output .mat filename for the correlation results.
%
%    Example:
%      Plot_KeyModes_vs_Scores(P, Scores_Table, Key_Modes_KC, results_dir, save_name);
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
%   3. Extract occupancies, harmonize across sites (optional), and choose which
%      occupancy matrix (raw or harmonized) the rest of the pipeline uses:
%     apply_combat = 1; occup_file = 'LEiDA_Occupancies_harmonized.mat';
%     [P_original, P_harmonized, rangeK, Scores_ADNI] = Save_Occupancies_Harmonize(results_dir, cluster_file, Scores_Table, apply_combat, occup_file);
%     use_harmonized_occupancies = 1;   % 1: use P_harmonized, 0: use P_original
%     P = P_harmonized;   % or P_original
%
%   4a. (Optional) Statistical analysis between discrete conditions - skip
%      this step entirely if your study only has continuous scores to correlate:
%     cond = {'CN','MCI','DEM'};
%     Paired_tests = 0; n_permutations = 1000; n_bootstraps = 0;
%     LEiDA_stats_Voxel_FracOccup_ComBat(results_dir, cluster_file, stats_file, cond, Index_Conditions, Paired_tests, n_permutations, n_bootstraps, P);
%
%   4b. Correlate every score against every mode in the entire pyramid (does
%      NOT depend on step 4a, or on any pre-selected set of modes):
%     pyramid_stats_file = 'Scores_Pyramid_Pval.mat';
%     Scores_vs_Mode_Occupancy(P, Scores_Table, results_dir, cluster_file, pyramid_stats_file);
%
%   5. Generate Figures:
%     a. Plot statistical results, including all p-values and effect sizes
%        (requires step 4a):
%        Plot_FracOccup_stats(results_dir, stats_file);
%
%     b. Centroid pyramid with RSN overlay, marked by a condition-pair result
%        (requires step 4a) - or pass pyramid_stats_file instead of stats_file
%        and a score index instead of a condition-pair index to mark by score
%        significance instead (requires step 4b):
%        overlap_RSNs = 0; cortex_dir = 'SideView'; Add_asterisks = 1; stat_of_interest = 3;
%        save_name = 'Centroid_Pyramid_Magenta';
%        Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, save_name, overlap_RSNs, cortex_dir, stat_of_interest, Add_asterisks);
%
%     c. Automatically select the key modes that differ between conditions
%        (requires step 4a). Without step 4a, specify Key_Modes_KC manually
%        instead, e.g. Key_Modes_KC = [2 4; 4 6] (1st column is the POSITION
%        in rangeK, not the literal K - see step 6's docs above):
%        Key_Modes_KC = Choose_Relevant_Modes(results_dir, cluster_file, stats_file);
%
%     d. Plot the key modes on anatomical slices with occupancy bars
%        (requires step 4a, for the Key_Modes_KC computed in c):
%        Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file, 'Fig1_Key_modes', Key_Modes_KC, Scores_Table);
%
%     e. Detailed 3D visualization of the key modes, with RSN overlap
%        (only requires cluster_file and Key_Modes_KC, not step 4a):
%        Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC);
%
%     f. Plot the key modes' occupancy against the scores from step 4b (only
%        requires P and Key_Modes_KC, not step 4a):
%        Plot_KeyModes_vs_Scores(P, Scores_Table, Key_Modes_KC, results_dir, 'Scores_Mode_Stats.mat');
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

%% 2b. Extract Mode Occupancies and Harmonize across sites (optional)

occup_file   = 'LEiDA_Occupancies_harmonized.mat';   % File to save original + harmonized occupancies
apply_combat = 1;   % 1 to harmonize occupancies across sites with ComBat; 0 to skip (P_harmonized = P_original)

[P_original, P_harmonized, rangeK, Scores_ADNI] = ...
    Save_Occupancies_Harmonize(results_dir, cluster_file, Scores_Table, apply_combat, occup_file);

% Choose, once, whether the rest of the pipeline (condition statistics in step
% 3a AND the score correlations in step 3b) runs on the raw or the harmonized
% occupancies. This is independent of whether step 3a is run at all, so
% step 3b can be used on its own for studies with no discrete conditions to
% test, only continuous scores to correlate with mode occupancy.
use_harmonized_occupancies = 1;   % 1: use P_harmonized; 0: use P_original
if use_harmonized_occupancies
    P = P_harmonized;
else
    P = P_original;
end

%% 3a. Statistical Analysis of Mode Occupancies between conditions
% (Skip this section entirely if your study has no discrete conditions to
% compare; step 3b below only needs P, not this section's output.)

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
    Condition_tags, Index_Conditions, Paired_tests, n_permutations, n_bootstraps, P);

%% 3b. Correlate every score with every mode in the entire pyramid
% Does NOT depend on step 3a, or on any pre-selected set of modes. Reports
% the most significant results to the command line and saves the p-values
% for use by Plot_ClustVoxelCentroid_Pyramid_RSNs.m (step 5b below) and by
% Plot_KeyModes_vs_Scores.m (step 5f below, once Key_Modes_KC is chosen).

pyramid_stats_file='Scores_Pyramid_Pval.mat';
Scores_vs_Mode_Occupancy(P,Scores_Table,results_dir,cluster_file,pyramid_stats_file)

%% 5. Figures with Results

% Statistical Report of Fractional Occupancy across conditions (requires 3a)
Plot_FracOccup_stats(results_dir, stats_file);
% Saves figures with main statistical results, including p-values, effect sizes and Barplot pyramid

%% 5b. Plot Centroid Pyramid with Statistical Results Overlay.
% Options:
%   - overlap_RSNs   : 1 to color the voxels using Yeo RSN colors; 0 otherwise.
%   - cortex_dir     : 'SideView' or 'TopView' to chose the camera angle of 3D brains
%   - Add_asterisks  : 1 to add significance markers; 0 otherwise (SideView only)
%   - stat_of_interest : Which statistical result to mark on each centroid.
%       With stats_file (requires 3a): index into a condition pair, e.g. 2 = CN vs. AD.
%       With pyramid_stats_file instead (requires 3b): index into a score,
%       e.g. see the command-line report from 3b for which scores had
%       significant modes.
overlap_RSNs = 1;
cortex_dir = 'SideView';
Add_asterisks = 0;

stat_of_interest = 2;
save_name = 'Fig2_Centroid_Pyramid_RSNs';

Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, save_name, overlap_RSNs, cortex_dir, stat_of_interest, Add_asterisks);

% OR, to mark by score significance instead of condition-comparison significance:
% Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, pyramid_stats_file, save_name, overlap_RSNs, cortex_dir, stat_of_interest, Add_asterisks);

%% 5c. Select a number of Key Modes to analyze

% Automatically select the Key Modes differing mostly between conditions
% (requires stats_file from step 3a, i.e. discrete conditions to compare):
Key_Modes_KC = Choose_Relevant_Modes(results_dir, cluster_file, stats_file);

% OR

% Select all Modes for a given K (e.g. K=8 clusters). Key_Modes_KC's 1st
% column must be the POSITION of that K in rangeK, not the literal K itself:
% K_wanted = 8;
% ki = find(rangeK == K_wanted);
% Key_Modes_KC=[ones(1,K_wanted)*ki; 1:K_wanted]';

% OR

% Make your own selection in pairs [ki c] (does not require step 3a at all,
% e.g. for studies with only continuous scores and no discrete conditions).
% ki is the POSITION in rangeK, not the literal K - e.g. with rangeK=2:20,
% ki=2 means K=3 clusters. Use find(rangeK==K_wanted) to convert, as above:
% Key_Modes_KC=[[2 4];[4 6];];

%% 5d. Plot the key modes differing between conditions (requires 3a for Key_Modes_KC above)

save_name = 'Fig1_Key_modes_Slice_Occupancy_bars_';
Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file,save_name,Key_Modes_KC,Scores_Table)

%% 5e. Plot detailed visualization of Key Modes and get the list of brain areas involved.

Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC);

%% 5f. Plot the key modes' occupancy against the scores from step 3b
% Only requires P and Key_Modes_KC, not step 3a.

save_name='Scores_Mode_Stats.mat';
Plot_KeyModes_vs_Scores(P,Scores_Table,Key_Modes_KC,results_dir,save_name)
