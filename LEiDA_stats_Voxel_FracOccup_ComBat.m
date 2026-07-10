function LEiDA_stats_Voxel_FracOccup_ComBat(results_dir, file_cluster, file_stats, cond, Index_Conditions, pair, n_permutations, n_bootstraps, P)
% LEiDA_stats_Voxel_FracOccup performs statistical tests on fractional
% occupancy of LEiDA coupling modes, comparing experimental conditions.
%
% For independent samples, Welch's t-test is used throughout (robust under
% both equal and unequal variances). For paired samples, a paired
% permutation test is used.
%
% INPUT:
%   results_dir      - Directory where LEiDA results are saved.
%   file_cluster     - Filename of the .mat file containing K-means results
%                      (must contain rangeK).
%   file_stats       - Filename to save statistical results.
%   P                - Fractional occupancy matrix (N_scans x length(rangeK) x rangeK(end).
%                      Pass either raw or ComBat/age-deconfounded values.
%   cond             - Cell array with condition tags.
%   Index_Conditions - Vector assigning each scan to a condition (same length as N_scans).
%   pair             - 0 for independent samples, 1 for paired samples.
%   n_permutations   - Number of permutation iterations (e.g. 10000).
%   n_bootstraps     - Number of bootstrap samples per permutation (e.g. 50).
%
% OUTPUT (saved to file_stats):
%   P                - Fractional occupancy matrix (as passed in).
%   P_pval           - Two-sided permutation p-values
%                      (n_pairs x length(rangeK) x rangeK(end)).
%   effectsize       - Hedge's effect size for each comparison.
%   cond, rangeK, file_cluster, Index_Conditions, pair,
%   condRow, condCol, n_bootstraps, n_permutations.
%
% Authors: Joana Cabral, Universidty of Lisbon,
% joanabcabral@tecnico.ulisboa.pt
%          Miguel Farinha, ICVS/2CA-Braga, miguel.farinha@ccabraga.pt

%% Load clustering parameters
load([results_dir file_cluster], 'rangeK');

N_scans = size(P, 1);
n_Cond  = size(cond, 2);

%% Define condition pairs for comparison
n_pairs  = n_Cond * (n_Cond - 1) / 2;
condRow  = zeros(1, n_pairs);
condCol  = zeros(1, n_pairs);

disp('Pairs of conditions compared:')
cond_pair = 1;
for cond1 = 1:n_Cond - 1
    for cond2 = cond1 + 1:n_Cond
        condRow(cond_pair) = cond1;
        condCol(cond_pair) = cond2;
        disp([num2str(cond_pair) ' : ' num2str(cond{cond1}) ' vs ' num2str(cond{cond2})])
        cond_pair = cond_pair + 1;
    end
end

%% Permutation tests
disp(' ');
disp(['Running permutation tests (' num2str(n_permutations) ' permutations, ' num2str(n_bootstraps) ' bootstraps each).']);
if pair == 0
    disp('Independent samples: using Welch''s t-test throughout.');
else
    disp('Paired samples: using paired permutation test.');
end
disp(' ');

P_pval     = zeros(n_pairs, length(rangeK), rangeK(end));
effectsize = zeros(n_pairs, length(rangeK), rangeK(end));

for k = 1:length(rangeK)
    if rangeK(k)>1
    disp(['K = ' num2str(rangeK(k))])
    for c = 1:rangeK(k)
        cond_pair = 1;
        for cond1 = 1:n_Cond - 1
            for cond2 = cond1 + 1:n_Cond

                a = squeeze(P(Index_Conditions == cond1, k, c))';
                b = squeeze(P(Index_Conditions == cond2, k, c))';
                design = [ones(1, numel(a)), 2 * ones(1, numel(b))];

                if pair == 1
                    stats = bootstrap_within_permutation_paired_samples( ...
                        [a, b], design, n_permutations, n_bootstraps, 0.05);
                else
                    stats = bootstrap_within_permutation_ttest2( ...
                        [a, b], design, n_permutations, n_bootstraps, 0.05, 'welchtest');
                end

                P_pval(cond_pair, k, c)     = stats.pvals_2sided;
                effectsize(cond_pair, k, c) = stats.eff;
                cond_pair = cond_pair + 1;
            end
        end
    end
    end
end

%% Save results
save([results_dir '/' file_stats], ...
    'P', 'P_pval', 'effectsize', ...
    'cond', 'rangeK', 'file_cluster', 'Index_Conditions', 'pair', ...
    'condCol', 'condRow', 'n_bootstraps', 'n_permutations');

disp(' ');
disp(['Results saved to ' file_stats]);
disp(' ');