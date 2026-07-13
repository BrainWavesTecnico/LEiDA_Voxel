function Scores_vs_Mode_Occupancy(P,Scores_Table,results_dir,cluster_file,pyramid_stats_file)
% Scores_vs_Mode_Occupancy correlates every clinical/cognitive score against
% every mode in the entire pyramid (every K in rangeK), partial correlation
% controlling for age. It does not depend on any pre-selected set of modes,
% and does not generate a figure - it only computes and saves the p-values,
% and reports the most significant results to the command line. Use
% Plot_KeyModes_vs_Scores.m afterwards to plot a selected set of modes (e.g.
% the output of Choose_Relevant_Modes.m) against these scores.
%
% P is taken directly as an input (e.g. P_original or P_harmonized from
% Save_Occupancies_Harmonize) rather than loaded from a LEiDA stats file, so
% this function does not depend on LEiDA_stats_Voxel_FracOccup_ComBat having
% been run. This makes it usable for studies with no discrete conditions to
% compare, only continuous scores to correlate with mode occupancy.
%
% For each score, prints to the command line whether any mode survives
% significance at 0.05/sum(rangeK)/N_scores (Bonferroni across the whole
% pyramid AND across all scores tested), and the most significant result.
% The p-values are saved in a format that Plot_ClustVoxelCentroid_Pyramid_RSNs.m
% can load in place of the permutation-test p-values (pass pyramid_stats_file
% as its stats_file, and use stat_of_interest to select which score to display).
%
% INPUT:
%   P                  - Fractional occupancy matrix (N_scans x length(rangeK) x rangeK(end)),
%                        e.g. P_original or P_harmonized from Save_Occupancies_Harmonize.
%   Scores_Table       : .mat file with the Scores_ADNI table.
%   results_dir        : Directory where the mat output is saved.
%   cluster_file       : Clustering results file name (used to load rangeK).
%   pyramid_stats_file : Output .mat filename for the entire-pyramid p-values,
%                        loadable by Plot_ClustVoxelCentroid_Pyramid_RSNs.m.
%
% NOTE: The set of score columns used (Genetics/Biomarkers/Cognitive_functions
% indices below) is hardcoded for the ADNI Scores_ADNI table used in Campo et
% al.; adapt these indices for a different scores table.
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

load(Scores_Table,'Scores_ADNI')
load([results_dir cluster_file], 'rangeK')

Genetics=[20];
Biomarkers=[21 23 24 25 26 27];
Cognitive_functions=[28 29  38 56 57 34 35 36 37 30 31 32 33  40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 ];

Selected_scores=flip([Genetics Biomarkers Cognitive_functions]);

Age=Scores_ADNI.AGE_AT_SCAN;

%% --- Entire pyramid: correlate every score with every mode, for every K ---
% Bonferroni threshold across the whole pyramid (sum(rangeK) modes) AND across
% all scores tested, matching the convention used elsewhere in this pipeline
% (e.g. Choose_Relevant_Modes.m uses 0.05/sum(rangeK) for a single comparison).
N_scores = length(Selected_scores);
N_modes_pyramid = sum(rangeK);
alpha_pyramid = 0.05 / N_modes_pyramid / N_scores;

Rho_pyramid  = zeros(N_scores, length(rangeK), rangeK(end));
Pval_pyramid = zeros(N_scores, length(rangeK), rangeK(end));

fprintf('\n=== Correlating %d scores with all %d modes across K=%d:%d (Bonferroni alpha = %.2e) ===\n', ...
    N_scores, N_modes_pyramid, rangeK(1), rangeK(end), alpha_pyramid);

for Score = 1:N_scores
    col = Scores_ADNI{:, Selected_scores(Score)};
    valid_values = find(~isnan(col));

    for ki = 1:length(rangeK)
        for c = 1:rangeK(ki)
            P_mode_vec = squeeze(P(:, ki, c));
            [Rho_pyramid(Score, ki, c), Pval_pyramid(Score, ki, c)] = ...
                partialcorr(col(valid_values), P_mode_vec(valid_values), Age(valid_values), 'rows', 'complete');
        end
    end

    score_name = Scores_ADNI.Properties.VariableNames{Selected_scores(Score)};
    Pval_score = squeeze(Pval_pyramid(Score, :, :));
    sig_mask = Pval_score > 0 & Pval_score < alpha_pyramid;

    if any(sig_mask(:))
        [ki_sig, c_sig] = find(sig_mask);
        [~, best] = min(Pval_score(sig_mask));
        fprintf('%s: SIGNIFICANT - %d mode(s) survive Bonferroni (min p = %.2e at K=%d, mode %d)\n', ...
            score_name, numel(ki_sig), Pval_score(ki_sig(best), c_sig(best)), rangeK(ki_sig(best)), c_sig(best));
    else
        fprintf('%s: not significant (min p = %.2e)\n', score_name, min(Pval_score(Pval_score > 0)));
    end
end

% --- Save in a format Plot_ClustVoxelCentroid_Pyramid_RSNs.m can load directly
% as stats_file: P_pval here has one "row" per score (instead of per condition
% pair), so pass stat_of_interest as a score index when calling that function.
P_pval = Pval_pyramid; %#ok<NASGU>
Score_labels = Scores_ADNI.Properties.VariableNames(Selected_scores); %#ok<NASGU>
cond = Score_labels; %#ok<NASGU>
condRow = 1:N_scores; %#ok<NASGU>
condCol = 1:N_scores; %#ok<NASGU>
Index_Conditions = Scores_ADNI.DX_num + 1; %#ok<NASGU>

save(fullfile(results_dir, pyramid_stats_file), ...
    'P', 'P_pval', 'Rho_pyramid', 'cond', 'condRow', 'condCol', 'Index_Conditions', ...
    'Score_labels', 'Selected_scores', 'rangeK', '-v7.3');
fprintf('\nPyramid-wide score p-values saved to %s\n', pyramid_stats_file);
