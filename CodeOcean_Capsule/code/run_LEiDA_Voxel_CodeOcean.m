%% run_LEiDA_Voxel_CodeOcean.m
%
% Code Ocean capsule script for the voxel-space LEiDA pipeline, starting from
% already-extracted leading eigenvectors. Meant to be run directly (not
% called as a function) from a Code Ocean capsule with the standard
% code/data/results layout (this file living in code/).
%
% Clusters the leading eigenvectors for K = 2:20, extracts mode occupancies
% (optionally ComBat-harmonized), tests for differences between conditions
% (Step 3a) and correlates mode occupancy with clinical/cognitive scores
% (Step 3b), and generates all summary figures (Step 4). Eigenvector
% extraction (Step 1: Get_EigenVectors_VoxelSpace_Server.m, which reads raw
% fMRI NIfTI files) is NOT run here - it is a separate script that needs the
% full-size raw fMRI data that isn't part of this capsule. Run it offline
% beforehand; this script only reads its output.
%
% Expected files in data_dir (rename file_V1/Scores_Table below to match):
%   'LEiDA_V1_all_MNI10mm_s90demo.mat' - eigenvectors, as saved by
%       Get_EigenVectors_VoxelSpace_Server.m: V1_all, ind_voxels,
%       MNI_lowres_Mask, data_info, Scan_num, Scan_length.
%   'Scores_ADNI_s90demo.mat' - table Scores_ADNI with SITE, AGE_AT_SCAN,
%       PTGENDER, PTEDUCAT, DX_num, DX columns (plus the score columns used
%       by Scores_vs_Mode_Occupancy.m), subsetted to the same scans as file_V1.
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

%% Directories - edit these if not running inside a Code Ocean capsule
data_dir    = '../data/';      % Directory with the eigenvector file and Scores table
results_dir = '../results/';   % Directory where all outputs are saved

if ~exist(results_dir, 'dir'), mkdir(results_dir); end

addpath(genpath(fileparts(mfilename('fullpath'))));

%% Step 1: Leading eigenvectors (extracted OFFLINE - not run here)
% Get_EigenVectors_VoxelSpace_Server.m reads raw fMRI NIfTI files and is a
% script, not a function, so it cannot be called from this capsule. Run it
% offline on your full dataset (or use Select_Demo_Subsample.m to build a
% small demo subsample from an existing full-cohort eigenvector file), then
% place the resulting eigenvector file and a matching Scores_ADNI table in
% data_dir - edit the two filenames below to match what you name them.
file_V1      = 'LEiDA_V1_all_MNI10mm_s90demo.mat';
Scores_Table = fullfile(data_dir, 'Scores_ADNI_s90demo.mat');

cluster_file = 'LEiDA_Clusters_VoxelMNI10mm_demo.mat';
occup_file   = 'LEiDA_Occupancies_demo.mat';
stats_file   = 'LEiDA_Stats_FracOccup_demo.mat';

%% Step 2: Cluster the eigenvectors for K = 2:20
mink = 2; maxk = 20; replicates = 20;   % fewer replicates than the full study, for capsule runtime
fprintf('\n=== Clustering leading eigenvectors (K = %d:%d) ===\n', mink, maxk);
LEiDA_cluster_VoxelMNI10mm(data_dir, file_V1, mink, maxk, replicates, results_dir, cluster_file);

%% Step 2b: Extract mode occupancies (ComBat harmonization optional)
% Off by default: with a small demo subsample, per-site N is typically too
% low for reliable ComBat estimation. Set to 1 if your demo data spans
% multiple sites with enough scans per site.
apply_combat = 0;
fprintf('\n=== Extracting mode occupancies (apply_combat = %d) ===\n', apply_combat);
[P_original, P_harmonized, rangeK, Scores_ADNI] = ...
    Save_Occupancies_Harmonize(results_dir, cluster_file, Scores_Table, apply_combat, occup_file);

if apply_combat
    P = P_harmonized;
else
    P = P_original;
end

%% Step 3a: Statistical analysis between conditions
Index_Conditions = Scores_ADNI.DX_num + 1; % (+ 1 so conditions are 1,2,3,...)
Condition_values = sort(unique(Index_Conditions));
Condition_tags = cell(1, length(Condition_values));

fprintf('\n=== Step 3a: Statistical analysis across conditions ===\n');
disp('Number of scans in each condition:')
for cnd = 1:length(Condition_values)
    Condition_tags{cnd} = Scores_ADNI.DX{find(Scores_ADNI.DX_num+1==Condition_values(cnd),1)};
    disp([num2str(cnd) '- ' Condition_tags{cnd} ' = ' num2str(sum(Index_Conditions==Condition_values(cnd)))])
end

Paired_tests = 0; n_permutations = 500; n_bootstraps = 0;
LEiDA_stats_Voxel_FracOccup_ComBat(results_dir, cluster_file, stats_file, ...
    Condition_tags, Index_Conditions, Paired_tests, n_permutations, n_bootstraps, P);

%% Step 3b: Correlate every score with every mode in the entire pyramid
% Does NOT depend on step 3a, or on any pre-selected set of modes.
fprintf('\n=== Step 3b: Comparing mode occupancy with scores ===\n');

pyramid_stats_file = 'Scores_Pyramid_Pval.mat';
Scores_vs_Mode_Occupancy(P, Scores_Table, results_dir, cluster_file, pyramid_stats_file);

%% Step 4: Figures
fprintf('\n=== Step 4: Generating figures ===\n');
Plot_FracOccup_stats(results_dir, stats_file);

overlap_RSNs = 1; cortex_dir = 'SideView'; Add_asterisks = 0; stat_of_interest = 1;
Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, cluster_file, stats_file, ...
    'Centroid_Pyramid_RSNs', overlap_RSNs, cortex_dir, stat_of_interest, Add_asterisks);

% Select the key modes that differ most between conditions (falling back to a
% fixed selection if none survive significance on this reduced sample).
Key_Modes_KC = Choose_Relevant_Modes(results_dir, cluster_file, stats_file);
if isempty(Key_Modes_KC)
    % Expected with a reduced demo sample: statistical power may be too low
    % for any mode to survive the significance threshold. Fall back to a
    % fixed selection (mid-range K) so the figure steps still produce output,
    % computing the same [ki c slope] format Plot_KeyModes_Slices_Stats expects
    % (ki is the POSITION in rangeK, not the literal number of clusters).
    warning('run_LEiDA_Voxel_CodeOcean:noSignificantModes', ...
        ['No mode survived the significance threshold on this demo subsample; ' ...
         'falling back to a fixed selection of modes for illustration.']);
    ki_demo = min(4, length(rangeK));
    n_demo = min(3, rangeK(ki_demo));
    Key_Modes_KC = zeros(n_demo, 3);
    for m = 1:n_demo
        mean_P_cond = zeros(1, length(Condition_tags));
        for j = 1:length(Condition_tags)
            mean_P_cond(j) = nanmean(P(Index_Conditions == j, ki_demo, m));
        end
        Key_Modes_KC(m,:) = [ki_demo, m, mean_P_cond(end) - mean_P_cond(1)];
    end
end

Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file, 'Fig1_Key_modes', Key_Modes_KC, Scores_Table);
Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC);
Plot_KeyModes_vs_Scores(P, Scores_Table, Key_Modes_KC, results_dir, 'Scores_Mode_Stats.mat');

fprintf('\n=== LEiDA Voxel pipeline complete. Results saved to %s ===\n', results_dir);
