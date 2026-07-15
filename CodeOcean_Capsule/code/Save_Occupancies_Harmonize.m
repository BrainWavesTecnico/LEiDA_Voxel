function [P_original, P_harmonized, rangeK, Scores_ADNI] = Save_Occupancies_Harmonize(results_dir, cluster_file, Scores_Table, apply_combat, occup_file)
% Save_Occupancies_Harmonize extracts the fractional occupancy of each mode
% (for every K) from a LEiDA_cluster_VoxelMNI10mm clustering solution, and
% optionally harmonizes it across acquisition sites with ComBat.
%
% Both the original (raw) and harmonized occupancies are saved together, so
% downstream steps (LEiDA_stats_Voxel_FracOccup_ComBat, Scores_vs_Mode_Occupancy)
% can be run on either one.
%
% INPUT:
%   results_dir  - Directory containing the cluster file.
%   cluster_file - Clustering results filename (output of LEiDA_cluster_VoxelMNI10mm;
%                  must contain Kmeans_results, rangeK, Scan_num, data_info).
%   Scores_Table - .mat file with the Scores_ADNI table (must contain SITE,
%                  AGE_AT_SCAN, PTGENDER, PTEDUCAT and DX_num columns when
%                  apply_combat is 1).
%   apply_combat - 1 to harmonize occupancies across SITE with ComBat, keeping
%                  diagnosis/age/sex/education as covariates of interest.
%                  0 to skip harmonization (P_harmonized is set equal to P_original).
%   occup_file   - Output filename (saved in results_dir) with the occupancies.
%
% OUTPUT (also saved to results_dir/occup_file):
%   P_original   - Raw fractional occupancy (N_scans x length(rangeK) x rangeK(end)).
%   P_harmonized - Occupancy after ComBat harmonization (same size), or equal
%                  to P_original when apply_combat is 0.
%   rangeK       - Vector containing the range of cluster numbers (K values) used.
%   Scores_ADNI  - The scores table loaded from Scores_Table (returned for
%                  convenience, e.g. to derive Index_Conditions downstream).
%
% Example:
%   [P_original, P_harmonized, rangeK, Scores_ADNI] = ...
%       Save_Occupancies_Harmonize(results_dir, cluster_file, Scores_Table, 1, 'LEiDA_Occupancies_harmonized.mat');
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

%% Load clustering results
load([results_dir cluster_file], 'Kmeans_results', 'rangeK', 'Scan_num', 'data_info');

n_scans = length(data_info);
unique_scans = unique(Scan_num);

%% --- Extract fractional occupancy for all k and c ---
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

%% Load scores table (needed for ComBat covariates, and returned for downstream use)
load(Scores_Table, 'Scores_ADNI');

%% --- Build P_original in the format used by the rest of the analysis ---
P_original = zeros(n_scans, length(rangeK), rangeK(end));
for ki = 1:length(rangeK)
    P_original(:, ki, 1:rangeK(ki)) = All_Occupancy{ki};
end

%% --- Harmonize occupancies across sites with ComBat (optional) ---
if apply_combat
    disp('Harmonizing occupancies across sites using ComBat.');

    % Covariates to remain unchanged
    age      = double(Scores_ADNI.AGE_AT_SCAN);
    sex      = double(Scores_ADNI.PTGENDER == "Male");
    edu      = double(Scores_ADNI.PTEDUCAT);
    diagnose = Scores_ADNI.DX_num;
    mod = [diagnose, age, sex, edu];

    % Variable to harmonize
    site = double(Scores_ADNI.SITE)';

    % ComBat needs at least 2 scans per site to estimate a within-batch
    % variance. With a single scan, combat.m's `var(s_data(:,indices)')`
    % sees a 1xK slice, which MATLAB's var() treats as a plain vector and
    % collapses to one scalar instead of one variance per feature/mode -
    % this crashes the vertcat of delta_hat across batches ("Dimensions of
    % arrays being concatenated are not consistent"). Common with a small
    % demo subsample spread across many sites. Merge any site with fewer
    % than 2 scans into one pooled "other sites" batch so every batch
    % ComBat sees has >= 2 members.
    site_levels = unique(site);
    site_counts = arrayfun(@(v) sum(site == v), site_levels);
    small_sites = site_levels(site_counts < 2);
    if ~isempty(small_sites)
        n_small_scans = sum(ismember(site, small_sites));
        if n_small_scans < 2
            error('Save_Occupancies_Harmonize:comBatTooFewScans', ...
                ['Only %d scan(s) belong to a site with fewer than 2 scans, and ' ...
                 'merging them still leaves fewer than 2 - ComBat cannot estimate a ' ...
                 'batch effect from a single scan. Increase the sample size, or set ' ...
                 'apply_combat=0.'], n_small_scans);
        end
        warning('Save_Occupancies_Harmonize:mergedSmallSites', ...
            ['%d site(s) with fewer than 2 scans (%d scans total) were merged into ' ...
             'one pooled "other sites" batch for ComBat, since a batch effect cannot ' ...
             'be estimated from a single scan.'], numel(small_sites), n_small_scans);
        site(ismember(site, small_sites)) = max(site) + 1;
    end

    All_Occupancies_harmonized = cell(length(rangeK), 1);
    for ki = 2:length(rangeK)
        fprintf('ComBat for k=%d\n', rangeK(ki))
        data_to_harmonize = All_Occupancy{ki}';   % k x n_scans
        All_Occupancies_harmonized{ki} = combat(data_to_harmonize, site, mod, 1)';
    end

    P_harmonized = zeros(n_scans, length(rangeK), rangeK(end));
    for ki = 1:length(rangeK)
        if ki == 1
            % First K in rangeK is left unharmonized (e.g. K=1 is trivially all-ones)
            P_harmonized(:, ki, 1:rangeK(ki)) = All_Occupancy{ki};
        else
            P_harmonized(:, ki, 1:rangeK(ki)) = All_Occupancies_harmonized{ki};
        end
    end

    disp('Occupancies harmonized for site differences');
else
    disp('Skipping ComBat harmonization (apply_combat=0): P_harmonized = P_original.');
    P_harmonized = P_original;
end

%% Save both original and harmonized occupancies
save(fullfile(results_dir, occup_file), 'P_original', 'P_harmonized', 'Scores_ADNI', 'rangeK', '-v7.3')
disp(['Occupancies saved to ' occup_file]);
