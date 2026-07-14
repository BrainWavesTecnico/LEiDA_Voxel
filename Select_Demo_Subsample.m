% Select_Demo_Subsample.m
%
% One-off, local data-prep script: selects a small demo subsample of scans,
% balanced by condition and, within each condition, stratified by sex and
% age-tertile (so age and sex distributions are matched as closely as
% possible across conditions) - one scan per unique participant. Reads your
% full leading eigenvectors + Scores_ADNI table, and saves the two files
% CodeOcean_Capsule/code/run_LEiDA_Voxel_CodeOcean.m expects in
% CodeOcean_Capsule/data/.
%
% This is NOT part of the analysis pipeline itself - run it once, locally,
% against your full dataset, to produce the small demo files you upload to
% the capsule. It never touches raw fMRI data, only the already-extracted
% leading eigenvectors (output of Get_EigenVectors_VoxelSpace_Server.m).
%
% Scores_ADNI IS reduced to just the selected scans (not kept at full size
% and indexed at runtime): the rest of the pipeline assumes Scores_ADNI has
% exactly one row per scan, in the same order as data_info/Scan_num, so a
% full-size table would misalign against the demo eigenvector file. Reducing
% it also avoids shipping clinical/genetic data for participants who aren't
% actually part of the demo.
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
%   - Scores_ADNI.AGE_AT_SCAN and Scores_ADNI.PTGENDER ('Male'/'Female')
%     exist, as used elsewhere in this pipeline (e.g. Save_Occupancies_Harmonize.m).
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

%% USER INPUT - edit to match your files/columns
full_data_dir   = '/path/to/your/full/results/';   % where the full-cohort files live
file_V1_full    = 'LEiDA_V1_all_MNI10mm_FullMask_2177scans.mat';
Scores_Table_full = fullfile(full_data_dir, 'Scores_ADNI_2177scans.mat');

ID_column   = 'RID';   % CHANGE to your table's unique-participant column, e.g. 'RID' or 'PTID'
n_per_condition = 30;  % scans per condition in the demo (30+30+30 = 90 total)
n_age_bins  = 3;       % age tertiles (young/mid/old) to stratify by, within each condition
rng_seed = 42;          % fixed seed for reproducibility

out_dir  = 'CodeOcean_Capsule/data/';            % local output folder for the capsule
file_V1_demo      = 'LEiDA_V1_all_MNI10mm_s90demo.mat';
Scores_Table_demo = 'Scores_ADNI_s90demo.mat';

%% Load full data
load(fullfile(full_data_dir, file_V1_full), ...
    'V1_all', 'ind_voxels', 'MNI_lowres_Mask', 'data_info', 'Scan_num', 'Scan_length');
load(Scores_Table_full, 'Scores_ADNI');

n_scans_full = length(data_info);
assert(height(Scores_ADNI) == n_scans_full, ...
    'Scores_ADNI must have exactly one row per scan (same order as data_info).');

%% Select n_per_condition unique-participant scans per condition,
%  stratified by sex x age-tertile so age/sex are balanced across conditions
rng(rng_seed);
Index_Conditions = Scores_ADNI.DX_num + 1;   % 1=CN, 2=MCI, 3=DEM (matches run_LEiDA_Voxel.m)
Condition_values = sort(unique(Index_Conditions));

Age = Scores_ADNI.AGE_AT_SCAN;
Sex = string(Scores_ADNI.PTGENDER);
Sex_values = unique(Sex(~ismissing(Sex)));

valid_demo = ~isnan(Age) & ~ismissing(Sex);
if any(~valid_demo)
    fprintf('Excluding %d scan(s) with missing age/sex from eligibility.\n', sum(~valid_demo));
end

% Age-tertile edges from the full (valid) cohort, so bin edges are the same
% across conditions - that's what makes the strata comparable between groups.
age_edges = quantile(Age(valid_demo), linspace(0, 1, n_age_bins + 1));
age_edges([1 end]) = [-Inf, Inf];
Age_bin = discretize(Age, age_edges);

n_strata = numel(Sex_values) * n_age_bins;
base_per_stratum = floor(n_per_condition / n_strata);
remainder = n_per_condition - base_per_stratum * n_strata;

selected_scans = [];
for cnd = 1:length(Condition_values)
    cond_scan_idx = find(Index_Conditions == Condition_values(cnd) & valid_demo);
    cond_ids = Scores_ADNI.(ID_column)(cond_scan_idx);
    [~, first_occurrence] = unique(cond_ids, 'stable');
    cond_scan_idx = cond_scan_idx(first_occurrence);   % one scan per participant, this condition

    cond_selected = [];
    strat_i = 0;
    for si = 1:numel(Sex_values)
        for ai = 1:n_age_bins
            strat_i = strat_i + 1;
            target_n = base_per_stratum + (strat_i <= remainder);   % spread the remainder over the first strata
            in_stratum = cond_scan_idx(Sex(cond_scan_idx) == Sex_values(si) & Age_bin(cond_scan_idx) == ai);
            n_take = min(target_n, numel(in_stratum));
            if n_take > 0
                pick = in_stratum(randperm(numel(in_stratum), n_take));
                cond_selected = [cond_selected; pick]; %#ok<AGROW>
            end
        end
    end

    % Backfill from the remaining eligible participants in this condition if
    % some strata didn't have enough people (relaxes strict age/sex matching
    % only for the shortfall, so the total N per condition is always met).
    shortfall = n_per_condition - numel(cond_selected);
    if shortfall > 0
        remaining_pool = setdiff(cond_scan_idx, cond_selected);
        if numel(remaining_pool) < shortfall
            error(['Only %d eligible participants available for condition %d, ' ...
                   'need %d. Lower n_per_condition, n_age_bins, or check ID_column.'], ...
                   numel(cond_scan_idx), Condition_values(cnd), n_per_condition);
        end
        backfill = remaining_pool(randperm(numel(remaining_pool), shortfall));
        cond_selected = [cond_selected; backfill];
        warning('Select_Demo_Subsample:strataShortfall', ...
            '%s: %d participant(s) backfilled outside strict age/sex stratification (some strata had too few eligible participants).', ...
            Scores_ADNI.DX{cond_scan_idx(1)}, shortfall);
    end

    selected_scans = [selected_scans; cond_selected]; %#ok<AGROW>

    sel_age = Age(cond_selected);
    sel_sex = Sex(cond_selected);
    fprintf('%s: %d participants selected, age %.1f +/- %.1f, %d male / %d female\n', ...
        Scores_ADNI.DX{cond_scan_idx(1)}, numel(cond_selected), mean(sel_age), std(sel_age), ...
        sum(sel_sex == "Male"), sum(sel_sex == "Female"));
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
% V1_all is saved as 'single' (halves the file size); LEiDA_cluster_VoxelMNI10mm.m
% converts it back to 'double' automatically before clustering.
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

V1_all = single(V1_all);
save(fullfile(out_dir, file_V1_demo), ...
    'V1_all', 'ind_voxels', 'MNI_lowres_Mask', 'data_info', 'Scan_num', 'Scan_length', '-v7.3');
save(fullfile(out_dir, Scores_Table_demo), 'Scores_ADNI', '-v7.3');

fprintf('Saved demo eigenvectors to %s and demo scores table to %s\n', ...
    fullfile(out_dir, file_V1_demo), fullfile(out_dir, Scores_Table_demo));
