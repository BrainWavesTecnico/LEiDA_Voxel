function Plot_FracOccup_stats(results_dir, stats_file)
% Plot_FracOccup_stats visualizes the results from hypothesis tests that compare
% the mean fractional occupancy of FC states between conditions.
%
% This function generates several figures:
%   Fig1: Plot of two-sided p-values obtained from permutation tests across K,
%         for each pair of conditions.
%   Fig2: Barplot of the mean fractional occupancy of each FC state for all conditions.
%   Fig3: Barplots of the mean fractional occupancy for each FC state for each
%         pair of conditions.
%   Fig4: Plot of Hedge's effect sizes from the permutation tests across K.
%
% INPUT:
%   results_dir - Directory with the results from the hypothesis tests on 
%                 the fractional occupancy of FC states.
%   stats_file  - Filename (MAT-file) containing statistics (e.g., 'LEiDA_Stats_Voxel_FracOccup.mat').
%
% OUTPUT:
%   The function creates and saves figures to disk (both PNG and FIG formats).
%
% Author: Joana Cabral, Tecnico, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

%% Load Required Data
% Load condition labels, fractional occupancy matrix (P), p-values (P_pval), index map,
% effect sizes, and the range of K values (number of FC states) from the statistics file.
load([results_dir stats_file], 'cond', 'P', 'P_pval', 'Index_Conditions', 'effectsize', 'rangeK');

% Determine the number of conditions in the experiment.
n_Cond = size(cond, 2);

%% Create Figure Layout Maps for Subplots
% Create a mapping to organize subplots for pairwise comparisons.
% The total number of condition pairs is (n_Cond*(n_Cond-1)/2).
pos = power(n_Cond-1, 2);
index_fig = reshape(1:pos, n_Cond-1, n_Cond-1).';
% Use the upper triangular part of a matrix to select subplot positions.
subplot_map = triu(ones(n_Cond-1));
% Find the indices in the subplot map where figures should be placed.
subplot_indices = find(subplot_map);

%% List Condition Pair Comparisons
% Create arrays that define each pair of conditions.
condRow = zeros(1, n_Cond*(n_Cond-1)/2);
condCol = zeros(1, n_Cond*(n_Cond-1)/2);
cond_pair = 1;
disp('Pairs of conditions compared:')
for cond1 = 1:n_Cond-1
    for cond2 = cond1+1:n_Cond
        condRow(cond_pair) = cond1;
        condCol(cond_pair) = cond2;
        disp([num2str(cond{cond1}) ' - ' num2str(cond{cond2})])
        cond_pair = cond_pair + 1;
    end
end

%% Define Color Tones for Plots
darkBlue = [0.07, 0.62, 1.00];
darkGreen = [0.3, 0.5, 0];
darkRed = [0.6, 0, 0];

%% Figure 1: Two-sided p-values from Hypothesis (Permutation) Tests
Fig1 = figure;
disp('Plotting two-sided p-values from hypothesis (permutation) tests:');
for s_ind = 1:length(subplot_indices)
    % Determine the subplot index to use.
    subplot_ind = subplot_indices(s_ind);
    subplot(size(subplot_map, 1), size(subplot_map, 2), index_fig(subplot_ind))
    
    % Plot reference lines:
    % Red dashed line for p = 0.05,
    % green dashed line for p = 0.05/K,
    % blue dashed line for p = 0.05/sum(K).
    semilogy(rangeK(1)-1:rangeK(end)+1, 0.05 * ones(1, length(rangeK)+2), '--','color', darkRed, 'LineWidth', 1.5)
    hold on
    semilogy(rangeK(1)-1:rangeK(end)+1, 0.05 ./ (rangeK(1)-1:rangeK(end)+1) .* ones(1, length(rangeK)+2), '--', 'Color', darkGreen, 'LineWidth', 1.5)
    semilogy(rangeK(1)-1:rangeK(end)+1, 0.05 ./ sum(rangeK) * ones(1, length(rangeK)+2), '--', 'color', darkBlue, 'LineWidth', 1.5)
    
    % Loop over each clustering solution (K) and state (c) to plot the p-values.
    for k = 1:length(rangeK)
        for c = 1:rangeK(k)
            % Plot markers with different symbols and colors based on significance thresholds.
            if P_pval(s_ind, k, c) > 0.05
                semilogy(rangeK(k), P_pval(s_ind, k, c), '.k', 'Markersize', 6);
            elseif P_pval(s_ind, k, c) < 0.05 && P_pval(s_ind, k, c) > (0.05 / rangeK(k))
                semilogy(rangeK(k), P_pval(s_ind, k, c), '+r', 'Markersize', 6);
            elseif P_pval(s_ind, k, c) < (0.05 / rangeK(k)) && P_pval(s_ind, k, c) > (0.05 / sum(rangeK))
                semilogy(rangeK(k), P_pval(s_ind, k, c), 'o', 'Markersize', 6, 'Color', darkGreen);
            elseif P_pval(s_ind, k, c) <= (0.05 / sum(rangeK))
                semilogy(rangeK(k), P_pval(s_ind, k, c), '*', 'color', darkBlue, 'Markersize', 6);
            end
        end
    end
     
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}], 'interpreter', 'none')
    ylabel({'Two-sided {\itp}-value', 'Fractional Occupancy'})
    xlabel('Number of FC States K')
    xticks([2 4 6 8 10 12 14 16 18 20 22 24])
    xlim([rangeK(1)-1 rangeK(end)+1])
    box off
    if min(P_pval(P_pval>0)) < 0.002
        ylim([min(P_pval(P_pval>0))/2 1])
    else
        ylim([0.001 1])
    end
    set(gca, 'YTick', 10.^(-10:1:0))
end
% Save figure 1
saveas(Fig1, fullfile(results_dir, 'Figure_FracOccup_pvalues.png'), 'png');
saveas(Fig1, fullfile(results_dir, 'Figure_FracOccup_pvalues.fig'), 'fig');
disp('- Plot successfully saved as FracOccup_pvalues');

%% Figure 2: Barplot of Mean Fractional Occupancy Across All Conditions
if n_Cond > 2
    Fig2 = figure;
    disp('Plotting barplot of the mean fractional occupancy across conditions and K:');
    for k_ind = 1:length(rangeK)
        for c = 1:rangeK(k_ind)
            % Create a subplot for the current K and state.
            subplot_tight(length(rangeK), rangeK(end), c + (k_ind-1) * rangeK(end),0.05)
            hold on
            
            % Compute mean occupancy and standard error (ste) for each condition.
            P_cond = cell(1, n_Cond);
            mean_P_cond = zeros(1, n_Cond);
            ste = zeros(1, n_Cond);
            for j = 1:n_Cond
                P_cond{j} = P(Index_Conditions == j, k_ind, c);
                mean_P_cond(j) = mean(P(Index_Conditions == j, k_ind, c));
                ste(j) = std(P(Index_Conditions == j, k_ind, c)) / sqrt(numel(P(Index_Conditions == j, k_ind, c)));
            end
            
            % Create a barplot with error bars.
            bar(1:n_Cond, mean_P_cond, 'EdgeColor', 'k', 'LineWidth', 0.8, 'FaceColor', 'none')
            hold on
            errorbar(mean_P_cond, ste, 'LineStyle', 'none', 'Color', 'k', 'CapSize', 5, 'LineWidth', 0.4);
            if k_ind == length(rangeK)
                set(gca, 'XTick', 1:n_Cond, 'XTickLabel', cond, 'Fontsize', 6, 'TickLabelInterpreter', 'none')
            else
                set(gca, 'XTick', 1:n_Cond, 'XTickLabel', [], 'Fontsize', 6, 'TickLabelInterpreter', 'none')
            end
            set(gca, 'color', 'none')
            hold off
            box off
        end
    end
    % Save figure 2
    saveas(Fig2, fullfile(results_dir, 'Figure_FracOccup_Barplot_Allconditions.png'), 'png');
    saveas(Fig2, fullfile(results_dir, 'Figure_FracOccup_Barplot_Allconditions.fig'), 'fig');
    disp('- Plot successfully saved as FracOccup_Barplot_Allconditions');
end

%% Figure 3: Barplots for Each Pair of Conditions
disp('Plotting barplots of mean fractional occupancy for each pair of conditions:');
n_compare = n_Cond * (n_Cond - 1) / 2;
for i = 1:n_compare
    Fig3 = figure('Name', ['Occupancy ' cond{condRow(i)} ' vs ' cond{condCol(i)}]);
    for k_ind = 1:length(rangeK)
        for c = 1:rangeK(k_ind)
            subplot_tight(length(rangeK), rangeK(end), c + (k_ind-1) * rangeK(end), 0.05)
            hold on
            
            % Extract occupancy for the two conditions from the current pair.
            P_cond1 = squeeze(P(Index_Conditions == condRow(i), k_ind, c))';
            P_cond2 = squeeze(P(Index_Conditions == condCol(i), k_ind, c))';
            
            % Choose bar edge color based on significance threshold.
            if P_pval(i, k_ind, c) > 0.05
                bar(1:2, [mean(P_cond1) mean(P_cond2)], 'EdgeColor', 'k', 'LineWidth', 0.8, 'FaceColor', 'none')
            elseif P_pval(i, k_ind, c) < 0.05 && P_pval(i, k_ind, c) > (0.05 / rangeK(k_ind))
                bar(1:2, [mean(P_cond1) mean(P_cond2)], 'EdgeColor', 'r', 'LineWidth', 0.8, 'FaceColor', 'none')
            elseif P_pval(i, k_ind, c) < (0.05 / rangeK(k_ind)) && P_pval(i, k_ind, c) > (0.05 / sum(rangeK))
                bar(1:2, [mean(P_cond1) mean(P_cond2)], 'EdgeColor', darkGreen, 'LineWidth', 0.8, 'FaceColor', 'none')
            elseif P_pval(i, k_ind, c) <= (0.05 / sum(rangeK))
                bar(1:2, [mean(P_cond1) mean(P_cond2)], 'EdgeColor', darkBlue, 'LineWidth', 0.8, 'FaceColor', 'none')
            end
            
            hold on
            % Plot error bars.
            errorbar([mean(P_cond1) mean(P_cond2)], [std(P_cond1)/sqrt(numel(P_cond1)) std(P_cond2)/sqrt(numel(P_cond2))],...
                     'LineStyle', 'none', 'Color', 'k', 'CapSize', 5, 'LineWidth', 0.4);
            if k_ind == length(rangeK)
                set(gca, 'XTick', 1:2, 'XTickLabel', {cond{condRow(i)} cond{condCol(i)}}, 'Fontsize', 6, 'TickLabelInterpreter', 'none')
            else
                set(gca, 'XTick', 1:2, 'XTickLabel', [])
            end
            set(gca, 'color', 'none')
            hold off
            box off
        end
    end
    % Save figure 3 for current pair of conditions.
    saveas(Fig3, fullfile(results_dir, ['Figure_FracOccup_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)} '.png']), 'png');
    saveas(Fig3, fullfile(results_dir, ['Figure_FracOccup_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)} '.fig']), 'fig');
    disp(['- Plot successfully saved as FracOccup_Barplot_' cond{condRow(i)} '_vs_' cond{condCol(i)}]);
end

%% Figure 4: Plot of Hedge's Effect Size from Hypothesis Tests
disp('Plotting Hedge''s effect size from hypothesis (permutation) tests:');
Fig4 = figure;
for s_ind = 1:length(subplot_indices)
    subplot_ind = subplot_indices(s_ind);
    subplot(size(subplot_map, 1), size(subplot_map, 2), index_fig(subplot_ind))
    
    % Plot reference lines for effect size thresholds.
    plot([rangeK(1)-1 rangeK(end)+1], [0.8 0.8], '--', 'color', darkBlue, 'LineWidth', 1.5);   
    hold on;
    plot([rangeK(1)-1 rangeK(end)+1], [0.5 0.5], '--', 'color', darkGreen, 'LineWidth', 1.5);
    plot([rangeK(1)-1 rangeK(end)+1], [0.2 0.2], 'r--', 'LineWidth', 1.5);
    
    % Loop over clustering solutions and states to plot effect sizes.
    for k_ind = 1:length(rangeK)
        for c = 1:rangeK(k_ind)
            if effectsize(s_ind, k_ind, c) < 0.2
                plot(rangeK(k_ind), effectsize(s_ind, k_ind, c), '.k', 'Markersize', 6);
            elseif effectsize(s_ind, k_ind, c) > 0.2 && effectsize(s_ind, k_ind, c) < 0.5
                plot(rangeK(k_ind), effectsize(s_ind, k_ind, c), '+r', 'Markersize', 6);
            elseif effectsize(s_ind, k_ind, c) > 0.5 && effectsize(s_ind, k_ind, c) < 0.8
                plot(rangeK(k_ind), effectsize(s_ind, k_ind, c), 'o', 'color', darkGreen, 'Markersize', 6);
            elseif effectsize(s_ind, k_ind, c) > 0.8
                plot(rangeK(k_ind), effectsize(s_ind, k_ind, c), '*', 'color', darkBlue, 'Markersize', 6);
            end
        end
    end
    
    title([cond{condRow(s_ind)} ' vs ' cond{condCol(s_ind)}], 'interpreter', 'none')
    ylabel('Hedge''s effect size')
    xlabel('Number of FC States K')
    xticks([2 4 6 8 10 12 14 16 18 20 22 24])
    xlim([rangeK(1)-1 rangeK(end)+1])
    box off
    ylim([0 1])
    set(gca, 'YTick', [0, 0.2, 0.5, 0.8, 1])
end
% Save figure 4.
saveas(Fig4, fullfile(results_dir, 'Figure_FracOccup_effetcsize.png'), 'png');
saveas(Fig4, fullfile(results_dir, 'Figure_FracOccup_effetcsize.fig'), 'fig');
disp('- Plot successfully saved as FracOccup_effetcsize');
disp(' ');
