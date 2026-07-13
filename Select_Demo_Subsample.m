% Select_Demo_Subsample.m
%
% One-off, local data-prep script: selects a small, balanced demo subsample
% (unique participants, equal N per condition) from your full leading
% eigenvectors + Scores_ADNI table, and saves the two files run_LEiDA_Voxel_CodeOcean.m
% expects in its data/ folder for a Code Ocean capsule.
%
% This is NOT part of the analysis pipeline itself - run it once, locally,
% against your full dataset, to produce the small demo files you upload to
% the capsule. It never touches raw fMRI data, only the already-extracted
% leading eigenvectors (output of Get_EigenVectors_VoxelSpace_Server.m).
%
% ASSUMPTIONS (adjust the USER INPUT section below to match your data):
%   - Scores_ADNI has one row per scan, in the same order as data_info /
%     unique(Scan_num) in the eigenvector file (this is required by the rest
%     of the pipeline too, e.g. Save_Occupancies_Harmonize.m).
%   - Scores_ADNI has a column that uniquely identifies the participant
%     (so the same person isn't sampled twice) - set ID_column below.
%   - Scores_ADNI.DX_num encodes condition as 0/1/2 (CN/MCI/DEM), as used
%     elsewhere in this pipeline (e.g. run_LEiDA_Voxel.m: Index_Conditions
%     = Scores_ADNI.DX_num + 1).
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

%% USER INPUT - edit to match your files/columns
full_data_dir   = '/path/to/your/full/results/';   % where the full-cohort files live
file_V1_full    = 'LEiDA_V1_all_MNI10mm_FullMask_2177scans.mat';
Scores_Table_full = fullfile(full_data_dir, 'Scores_ADNI_2177scans.mat');

ID_column   = 'RID';   % CHANGE to your table's unique-participant column, e.g. 'RID' or 'PTID'
n_per_condition = 30;  % scans per condition in the demo (30+30+30 = 90 total)
rng_seed = 42;          % fixed seed for reproducibility

out_dir  = 'data/';                              % local output folder for the capsule
file_V1_demo      = 'LEiDA_V1_all_MNI10mm_demo.mat';
Scores_Table_demo = 'Scores_ADNI_demo.mat';

%% Load full data
load(fullfile(full_data_dir, file_V1_full), ...
    'V1_all', 'ind_voxels', 'MNI_lowres_Mask', 'data_info', 'Scan_num', 'Scan_length');
load(Scores_Table_full, 'Scores_ADNI');

n_scans_full = length(data_info);
assert(height(Scores_ADNI) == n_scans_full, ...
    'Scores_ADNI must have exactly one row per scan (same order as data_info).');

%% Select n_per_condition unique-participant scans per condition
rng(rng_seed);
Index_Conditions = Scores_ADNI.DX_num + 1;   % 1=CN, 2=MCI, 3=DEM (matches run_LEiDA_Voxel.m)
Condition_values = sort(unique(Index_Conditions));

selected_scans = [];
for cnd = 1:length(Condition_values)
    cond_scan_idx = find(Index_Conditions == Condition_values(cnd));
    cond_ids = Scores_ADNI.(ID_column)(cond_scan_idx);

    [unique_ids, first_occurrence] = unique(cond_ids, 'stable');
    n_available = numel(unique_ids);
    if n_available < n_per_condition
        error(['Only %d unique participants available for condition %d, ' ...
               'need %d. Lower n_per_condition or check ID_column.'], ...
               n_available, Condition_values(cnd), n_per_condition);
    end

    pick = randperm(n_available, n_per_condition);
    selected_scans = [selected_scans; cond_scan_idx(first_occurrence(pick))]; %#ok<AGROW>

    disp([Scores_ADNI.DX{cond_scan_idx(first_occurrence(pick(1)))} ...
          ': selected ' num2str(n_per_condition) ' of ' num2str(n_available) ' unique participants']);
end

selected_scans = sort(selected_scans);
n_scans_demo = numel(selected_scans);
fprintf('Selected %d scans total (%d unique participants, %d conditions x %d each).\n', ...
    n_scans_demo, n_scans_demo, length(Condition_values), n_per_condition);

%% Subset Scores_ADNI, data_info, Scan_length (row-aligned to selected_scans)
Scores_ADNI = Scores_ADNI(selected_scans, :);
data_info   = data_info(selected_scans);
Scan_length = Scan_length(selected_scans);

%% Subset V1_all / Scan_num (row-per-TR), remapping scan numbers to 1:n_scans_demo
keep_TR = ismember(Scan_num, selected_scans);
V1_all_demo  = V1_all(keep_TR, :);
Scan_num_old = Scan_num(keep_TR);

Scan_num_demo = zeros(size(Scan_num_old));
for new_id = 1:n_scans_demo
    Scan_num_demo(Scan_num_old == selected_scans(new_id)) = new_id;
end

V1_all = V1_all_demo;
Scan_num = Scan_num_demo;

%% Save the demo files
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

save(fullfile(out_dir, file_V1_demo), ...
    'V1_all', 'ind_voxels', 'MNI_lowres_Mask', 'data_info', 'Scan_num', 'Scan_length', '-v7.3');
save(fullfile(out_dir, Scores_Table_demo), 'Scores_ADNI', '-v7.3');

fprintf('Saved demo eigenvectors to %s and demo scores table to %s\n', ...
    fullfile(out_dir, file_V1_demo), fullfile(out_dir, Scores_Table_demo));
