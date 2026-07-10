function Plot_KeyModes_Slices_Stats(results_dir, cluster_file, stats_file,save_name,Key_Modes_KC)
% Plot_KeyModes_Slices_Stats renders each selected [k c] key mode and displays
% mean +/- SE fractional-occupancy bars per condition, on anatomical slices.
%
% This function:
%   - Loads clustering and statistical data.
%   - Sorts the centroids by their occupancy in the control condition.
%   - Plots a barplot with error bars (for each condition) for each centroid.
%   - Renders the corresponding centroid as a 3D image using transparent brain masks
%     and anatomical slices in axial, coronal, and sagittal views.
%
% INPUT:
%   results_dir  - Directory where the cluster and stats files are stored.
%   cluster_file - Filename containing clustering results (e.g., centroids, mask).
%   stats_file   - Filename containing statistical results (occupancy, p-values, etc.).
%   save_name    - Base name used when saving the output figure.
%   Key_Modes_KC - Nx2+ matrix with one row per key mode, [k c ...], as returned
%                  by Choose_Relevant_Modes (or built manually).
%
% NOTE: this function also loads a hardcoded 'Scores_ADNI_2177scans.mat' from
% results_dir to split subjects by sex; that file is study-specific and not
% included in this repository.
%
% Author: Joana Cabral, University of Lisbon,
% joanabcabral@tecnico.ulisboapt
 
%% Load Required Data
% Load clustering centroids, rangeK, brain mask, and voxel indices.
load([results_dir cluster_file], 'Kmeans_results', 'rangeK', 'MNI_lowres_Mask', 'ind_voxels');
% Load statistical data (fractional occupancy etc.).
load([results_dir stats_file], 'cond','P','P_pval', 'Index_Conditions', 'effectsize', 'condRow', 'condCol');

load([results_dir 'Scores_ADNI_2177scans.mat'],'Scores_ADNI')
Index_Sex=Scores_ADNI.PTGENDER=='Male';

%% Setup Mask and Color Map
% Define expected output volume dimensions.
size_MNI = [91 109 91];
lim_mask = 0.45;

% Load a predefined colorbar for visualization.
load ColormapBlueRed cmap_blue_red

% Get number of conditions.
n_Cond = size(cond, 2);

% Create and smooth the brain mask, then resize it to match target dimensions.
Mask_Brain = zeros(size(MNI_lowres_Mask));
Mask_Brain(ind_voxels) = 1;
Mask_Brain = flip(smooth3(Mask_Brain, 'gaussian', 3, 0.8));
Mask_Brain = imresize3(Mask_Brain, size_MNI, 'Method', 'linear');

%% Create Figure of Modes Decreasing Occupancy
Key_Modes_Decrease=Key_Modes_KC(Key_Modes_KC(:,3)<0,:);
N_Modes=size(Key_Modes_Decrease,1);

Fig = figure('Color'); %, 'w', 'Position',[14          -1         765         340*N_Modes]);
colormap(cmap_blue_red);

n_slices=9;
n_columns = 3 + n_slices;

% Loop over each centroid (sorted) for the selected K.
for Mode = 1:N_Modes

    k=Key_Modes_Decrease(Mode,1);
    c=Key_Modes_Decrease(Mode,2);

    % Create subplot for barplot (error bars)
    subplot_tight(N_Modes, n_columns, (1:3) + (Mode-1)*n_columns,0.1)
    
    % Initialize containers for occupancy data for each condition.
    P_cond = cell(1, n_Cond);
    mean_P_cond = zeros(1, n_Cond);
    ste = zeros(1, n_Cond);
    for j = 1:n_Cond
        P_cond{j} = P(Index_Conditions == j, rangeK == k, c);
        mean_P_cond(j) = nanmean(P(Index_Conditions == j, rangeK == k, c));
        ste(j) = std(P(Index_Conditions == j, rangeK == k, c)) / sqrt(numel(P(Index_Conditions == j, rangeK == k, c)));
        % FEMALE
        P_condF{j} = P(~Index_Sex & Index_Conditions == j, rangeK == k, c);
        mean_P_condF(j) = nanmean(P(~Index_Sex & Index_Conditions == j, rangeK == k, c));
        steF(j) = std(P(~Index_Sex & Index_Conditions == j, rangeK == k, c)) / sqrt(numel(P(~Index_Sex & Index_Conditions == j, rangeK == k, c)));
        % MALE
        P_condM{j} = P(Index_Sex & Index_Conditions == j, rangeK == k, c);
        mean_P_condM(j) = nanmean(P(Index_Sex & Index_Conditions == j, rangeK == k, c));
        steM(j) = std(P(Index_Sex & Index_Conditions == j, rangeK == k, c)) / sqrt(numel(P(Index_Sex & Index_Conditions == j, rangeK == k, c)));
    end
    hold on
    x = 1:n_Cond;
    offset = 0.2;  % adjust this to taste
    errorbar(x, mean_P_cond, ste,"_",'MarkerSize',8,'LineStyle', 'none', 'Color', 'k', 'CapSize', 6, 'LineWidth', 1);
%    errorbar(x - offset, mean_P_condF, steF,"_",'MarkerSize',6,'LineStyle', 'none', 'Color', [0.95,0.40,0.77], 'CapSize', 4, 'LineWidth', .5);
%    errorbar(x + offset, mean_P_condM, steM,"_",'MarkerSize',6,'LineStyle', 'none', 'Color', [0.18,0.75,0.94], 'CapSize', 4, 'LineWidth', .5);
    set(gca, 'XTick', 1:n_Cond, 'XTickLabel', cond, 'Fontsize', 8, 'TickLabelInterpreter', 'none')
    set(gca, 'color', 'none')
%    legend({'All', 'F', 'M'}, 'Location', 'northwest', 'FontSize', 8)
    xlim([.4 n_Cond+.6])
    ylabel('Occupancy')
    hold off, box off
    
    
    % Retrieve the centroid vector and reshape it to a 3D volume.
    Vc = Kmeans_results{rangeK==k}.C(c, :);
    Vc_3D = zeros(size(MNI_lowres_Mask));
    Vc_3D(ind_voxels) = Vc;
    Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.4);
    Vc_3D = imresize3(Vc_3D, size_MNI, 'Method', 'linear');
    
    % Zero out low-intensity values using the brain mask.
    Vc_3D(Mask_Brain < lim_mask) = 0;
    Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.8);
    
    % Define display limits based on centroid variability.
    lim = 2 * std(Vc);
    
    %%% Plot Slices: Axial, Coronal, and Sagittal %%%
    sagital_start = 22;
    sagital_end = 72;
    sagital_slices=round(sagital_start:( sagital_end-sagital_start)/(n_slices-1): sagital_end);
   
    for slice = 1:length(sagital_slices)
        subplot_tight(N_Modes*3, n_columns, 3+ slice + (Mode-1)*3*(n_columns), 0.0005)
        imagesc(squeeze(Vc_3D(sagital_slices(slice),:,:))', 'AlphaData',squeeze(Mask_Brain(sagital_slices(slice),:,:)>lim_mask)',[-lim lim])
        axis xy, axis image, xlim([5 105]), ylim([0 120]), xticks(''), yticks(''), axis off
    end
    
    axial_start=73;
    axial_end=13;
    axial_slices = round(axial_start:( axial_end-axial_start)/(n_slices-1): axial_end);
    for slice = 1:length(axial_slices)
        subplot_tight(N_Modes*3, n_columns, 3+slice + (Mode-1)*3*(n_columns)+ n_columns, 0.0005)
        imagesc(squeeze(Vc_3D(:,:,axial_slices(slice)))', 'AlphaData',squeeze(Mask_Brain(:,:,axial_slices(slice))>lim_mask)',[-lim lim])
        axis xy, axis image, xlim([0 91]), ylim([8 105]), xticks(''), yticks(''), axis off
    end
    
    coronal_start=84;
    coronal_end=22;
    coronal_slices = round(coronal_start:(coronal_end-coronal_start)/(n_slices-1): coronal_end);
    for slice = 1:length(coronal_slices)
        subplot_tight(N_Modes*3, n_columns, 3+slice + (Mode-1)*3*(n_columns)+ 2*n_columns, 0.0005)
        imagesc(squeeze(Vc_3D(:,coronal_slices(slice),:))', 'AlphaData',squeeze(Mask_Brain(:,coronal_slices(slice),:)>lim_mask)',[-lim lim])
        axis xy, axis image, xlim([0 90]), ylim([-10 92]), xticks(''), yticks(''), axis off
    end
    
    % Display results for the current centroid.
    disp(['Results for k = ' num2str(k) ' c = ' num2str(c)]);
    for condpair = 1:size(P_pval, 1)
        disp([num2str(condpair) ' : ' num2str(cond{condRow(condpair)}) ' - ' num2str(cond{condCol(condpair)})]);
        disp(['     permutation p-value = ' num2str(P_pval(condpair, rangeK == k, c), 3)]);
        disp(['     effect size = ' num2str(effectsize(condpair, rangeK == k, c), 3)]);
    end
end

%% Save the Figure
% Save the figure as both PNG and MATLAB FIG file.
saveas(Fig, fullfile(results_dir, ['_' save_name '_Decrease.png']), 'png');
saveas(Fig, fullfile(results_dir, [save_name '_Decrease.fig']), 'fig');
disp(['- Plot successfully saved as ' save_name]);
disp(' ');

%% Create Figure of Modes Increasing Occupancy
Key_Modes_Increase=Key_Modes_KC(Key_Modes_KC(:,3)>0,:);
N_Modes=size(Key_Modes_Increase,1);

Fig = figure('Color', 'w'); % , 'Position',[14          -1         765        340*N_Modes]);
colormap(cmap_blue_red);

% Loop over each centroid (sorted) for the selected K.
for Mode = 1:N_Modes

    k=Key_Modes_Increase(N_Modes-Mode+1,1);   
    c=Key_Modes_Increase(N_Modes-Mode+1,2);

    % Create subplot for barplot (error bars)
    subplot_tight(N_Modes, n_columns, (1:3) + (Mode-1)*n_columns,0.1)
    
    
    % Initialize containers for occupancy data for each condition.
    P_cond = cell(1, n_Cond);
    mean_P_cond = zeros(1, n_Cond);
    ste = zeros(1, n_Cond);
    for j = 1:n_Cond
        P_cond{j} = P(Index_Conditions == j, rangeK == k, c);
        mean_P_cond(j) = nanmean(P(Index_Conditions == j, rangeK == k, c));
        ste(j) = std(P(Index_Conditions == j, rangeK == k, c)) / sqrt(numel(P(Index_Conditions == j, rangeK == k, c)));
        % FEMALE
        P_condF{j} = P(~Index_Sex & Index_Conditions == j, rangeK == k, c);
        mean_P_condF(j) = nanmean(P(~Index_Sex & Index_Conditions == j, rangeK == k, c));
        steF(j) = std(P(~Index_Sex & Index_Conditions == j, rangeK == k, c)) / sqrt(numel(P(~Index_Sex & Index_Conditions == j, rangeK == k, c)));
        % MALE
        P_condM{j} = P(Index_Sex & Index_Conditions == j, rangeK == k, c);
        mean_P_condM(j) = nanmean(P(Index_Sex & Index_Conditions == j, rangeK == k, c));
        steM(j) = std(P(Index_Sex & Index_Conditions == j, rangeK == k, c)) / sqrt(numel(P(Index_Sex & Index_Conditions == j, rangeK == k, c)));
    end
    hold on
        x = 1:n_Cond;
    offset = 0.2;  % adjust this to taste
    errorbar(x, mean_P_cond, ste,"_",'MarkerSize',8,'LineStyle', 'none', 'Color', 'k', 'CapSize', 6, 'LineWidth', 1);
%    errorbar(x - offset, mean_P_condF, steF,"_",'MarkerSize',6,'LineStyle', 'none', 'Color', [0.95,0.40,0.77], 'CapSize', 4, 'LineWidth', .5);
%    errorbar(x + offset, mean_P_condM, steM,"_",'MarkerSize',6,'LineStyle', 'none', 'Color', [0.18,0.75,0.94], 'CapSize', 4, 'LineWidth', .5);
     set(gca, 'XTick', 1:n_Cond, 'XTickLabel', cond, 'Fontsize', 8, 'TickLabelInterpreter', 'none')
    set(gca, 'color', 'none')
    xlim([.4 n_Cond+.6])
    ylabel('Occupancy')
%    legend({'All', 'F', 'M'}, 'Location', 'northwest', 'FontSize', 8)
    hold off, box off
    
    % Retrieve the centroid vector and reshape it to a 3D volume.
    Vc = Kmeans_results{k}.C(c, :);
    Vc_3D = zeros(size(MNI_lowres_Mask));
    Vc_3D(ind_voxels) = Vc;
    Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.4);
    Vc_3D = imresize3(Vc_3D, size_MNI, 'Method', 'linear');
    
    % Zero out low-intensity values using the brain mask.
    Vc_3D(Mask_Brain < lim_mask) = 0;
    Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.8);
    
    % Define display limits based on centroid variability.
    lim = 2 * std(Vc);
    
    %%% Plot Slices: Axial, Coronal, and Sagittal %%%
    %sagital_slices = 22:4:72;
    for slice = 1:length(sagital_slices)
        subplot_tight(N_Modes*3, n_columns, 3+ slice + (Mode-1)*3*(n_columns), 0.0005)
        imagesc(squeeze(Vc_3D(sagital_slices(slice),:,:))', 'AlphaData',squeeze(Mask_Brain(sagital_slices(slice),:,:)>lim_mask)',[-lim lim])
        axis xy, axis image, xlim([5 105]), ylim([0 120]), xticks(''), yticks(''), axis off
    end
    
    
    %axial_slices = 73:-5:10;
    for slice = 1:length(axial_slices)
        subplot_tight(N_Modes*3, n_columns, 3+slice + (Mode-1)*3*(n_columns)+ n_columns, 0.0005)
        imagesc(squeeze(Vc_3D(:,:,axial_slices(slice)))', 'AlphaData',squeeze(Mask_Brain(:,:,axial_slices(slice))>lim_mask)',[-lim lim])
        axis xy, axis image, xlim([0 91]), ylim([8 105]), xticks(''), yticks(''), axis off
    end
    
    for slice = 1:length(coronal_slices)
        subplot_tight(N_Modes*3, n_columns, 3+slice + (Mode-1)*3*(n_columns)+ 2*n_columns, 0.0005)
        imagesc(squeeze(Vc_3D(:,coronal_slices(slice),:))', 'AlphaData',squeeze(Mask_Brain(:,coronal_slices(slice),:)>lim_mask)',[-lim lim])
        axis xy, axis image, xlim([0 90]), ylim([-10 92]), xticks(''), yticks(''), axis off
    end
    
    % Display results for the current centroid.
    disp(['Results for k = ' num2str(k) ' c = ' num2str(c)]);
    for condpair = 1:size(P_pval, 1)
        disp([num2str(condpair) ' : ' num2str(cond{condRow(condpair)}) ' - ' num2str(cond{condCol(condpair)})]);
        disp(['     permutation p-value = ' num2str(P_pval(condpair, rangeK == k, c), 3)]);
        disp(['     effect size = ' num2str(effectsize(condpair, rangeK == k, c), 3)]);
    end
end

%% Save the Figure
% Save the figure as both PNG and MATLAB FIG file.
saveas(Fig, fullfile(results_dir, ['_' save_name '_Increase.png']), 'png');
saveas(Fig, fullfile(results_dir, [save_name '_Increase.fig']), 'fig');
disp(['- Plot successfully saved as ' save_name '_Increase']);
disp(' ');
