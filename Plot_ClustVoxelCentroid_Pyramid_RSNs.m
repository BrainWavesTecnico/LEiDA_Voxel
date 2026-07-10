function Plot_ClustVoxelCentroid_Pyramid_RSNs(results_dir, file_clusters, stats_file, save_name, overlap_Yeo, cortex_dir, cond_pair, Add_asterisks)
% Plot_ClustVoxelCentroid_Pyramid_RSNs renders all the clustering centroids in cortical
% space.
%
% The function loads the clustering results and related statistical
% results,  and creates a rendered image displaying the centroids on a transparent brain.
% Optionally, it colors the modes based on correlation with RSNs and/or marks centroids based on significance levels.
%
% INPUT:
%   results_dir   - Directory where the clustering centroids and statistical results are saved.
%   file_clusters - Filename (.mat) containing K-means clustering results (centroids, rangeK, mask).
%   stats_file    - Filename (.mat) with statistical results and occupancy values (P, P_pval, cond, etc.).
%   save_name     - Base name for saving the output figure.
%   overlap_Yeo   - Flag (1/0) to indicate whether to color voxels using Yeo network colors.
%   cortex_dir    - View direction for rendering (e.g., 'TopView' or 'SideView').
%   cond_pair     - Index of the condition pair for which p-values will be reported.
%   Add_asterisks - Flag to indicate whether to overlay significance markers on the plot.
%
% OUTPUT:
%   Fig           - A figure rendering all the centroids, saved in .fig and .jpg
%
% Author: Joana Cabral, Tecnico, University of Lisbon, joanabcabral@tecnico.ulisboa.pt
% Version from 9 April 2025

%% Load Required Data

% Load clustering centroids and masks from the clusters file.
load([results_dir file_clusters], 'Kmeans_results', 'rangeK', 'MNI_lowres_Mask', 'ind_voxels');

% Load statistical occupancy data and condition info.
load([results_dir stats_file], 'P', 'P_pval', 'cond', 'condCol', 'condRow', 'Index_Conditions');
n_Cond=size(P_pval,1);

% Display condition pair information (e.g., "Condition A vs Condition B")
disp(['P-values reported for ' cond{condRow(cond_pair)} ' vs ' cond{condCol(cond_pair)}]);

% Extract p-values for the specified condition pair for each clustering solution.
P_pval = squeeze(P_pval(cond_pair, :, :));

%% Load RSN Parcellation and Prepare Mask
% Load the mask of the Yeo parcellation in MNI 2mm space.
V_Yeo = struct2array(load('ParcelsMNI2mm', 'V_Yeo7'));

if overlap_Yeo
    % Identify voxels belonging to any of the 7 Yeo networks.
    ind_voxels_Yeo = find(V_Yeo > 0 & V_Yeo < 8);
    % Define Yeo network colors (scaled to [0,1]) based on original paper
    YeoColor = [180 115 208; 76 112 210; 100 174 117; 230 150 200; 207 225 185; 250 221 107;   227 120 91] ./ 256;
end

% Define additional color tones for significance markers.
darkBlue = [0.07, 0.62, 1.00];
darkGreen = [0.3, 0.5, 0];
darkRed = [0.6, 0, 0];

% Determine the size of the MNI volume to later resize masks appropriately.
size_MNI = size(V_Yeo);

% Create a smoothed, resampled brain mask for the transparent brain
lim_mask = 0.5;
Mask_Brain = zeros(size(MNI_lowres_Mask));
Mask_Brain(ind_voxels) = 1;
Mask_Brain = flip(smooth3(Mask_Brain, 'gaussian', 3, 0.8));
Mask_Brain = imresize3(Mask_Brain, size_MNI, 'Method', 'linear');

%% Render and Plot Clusters
% Create a full-screen figure for rendering.
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT CENTROIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(' ');

disp(['Rendering the centroids (using ' cortex_dir '):'])
Fig = figure('Position', get(0, 'Screensize'));

% Loop over each clustering solution based on different K values.
for k = 1:length(rangeK)
    disp(['- K = ' num2str(rangeK(k))])
    % Loop over each state within the current clustering solution.
    for Centroid = 1:rangeK(k)
        % Get the centroid order from the reordering performed earlier.
        %Centroid = c_reorder{k}(c);
        
        % Retrieve the centroid's vector and reshape it into 3D volume using voxel indices.
        Vc = Kmeans_results{k}.C(Centroid, :);
        Vc_3D = zeros(size(MNI_lowres_Mask));
        Vc_3D(ind_voxels) = Vc;
        % Smooth and resize the centroid volume to match the MNI mask.
        Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.4);
        Vc_3D = imresize3(Vc_3D, size_MNI, 'Method', 'linear');

        % Compute overlap with Yeo networks (if selected)
        if overlap_Yeo
            for net = 1:7
                dice_net(net) = dice(V_Yeo(ind_voxels_Yeo) == net, Vc_3D(ind_voxels_Yeo) > 0);
            end
        end

        % Threshold and further smooth the centroid volume by masking out
        % values outside the transparent brain.
        Vc_3D(Mask_Brain < lim_mask) = 0;
        Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.8);

        % Create a subplot for the current cluster state.
        subplot_tight(length(rangeK), rangeK(end)*2, Centroid*2 - 1 + (k - 1) * rangeK(end)*2, 0.01);
        hold on
        % Plot a transparent brain surface for context.
        Brainpatch = patch(isosurface(Mask_Brain, 0.3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

        reducepatch(Brainpatch, 0.4);
        isonormals(Mask_Brain, Brainpatch);

        % Plot the centroid patch.
        if overlap_Yeo && sum(dice_net > 0.55)
            % Color based on Yeo network color corresponding to maximum overlap.
            Mode_patch = patch(isosurface(smooth3(Vc_3D > 0), 0.6), 'FaceColor', YeoColor(dice_net == max(dice_net), :), 'EdgeColor', 'none');
        elseif overlap_Yeo
            % Render in grey if no significant overlap with any network.
            Mode_patch = patch(isosurface(smooth3(Vc_3D > 0), 0.6), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
        else
            % Otherwise, use a default magenta color.
            Mode_patch = patch(isosurface(smooth3(Vc_3D > 0), 0.6), 'FaceColor', 'm', 'EdgeColor', 'none');
        end

        reducepatch(Mode_patch, 0.5);
        isonormals(Vc_3D, Mode_patch);

        % Configure material and lighting for the rendered patch.
        material dull;
        lighting gouraud;
        % Set the view according to the requested direction.
        switch cortex_dir
            case 'TopView'
                view(-90, 90);
            case 'SideView'
                view(0, 0);
        end

        daspect([1 1 1]);
        camlight;
        % Set axis limits and hide axes.
        xlim([0 110]); ylim([0 95]); zlim([0 90]);
        axis off;

        % If required, add significance markers based on p-values.
        pval = P_pval(k, Centroid);
        clear t
        if Add_asterisks

            % Different markers and colors are chosen based on p-value thresholds.
            if  pval <= (0.05 / sum(rangeK)) / 10 / 10 / 10 / 10 / 10 / 10
                % Draw multiple markers when p-values are extremely significant.
                t=title('*******', 'Color', darkBlue);
                t.Position(1)=t.Position(1)*0.8;
                t.Position(3)=t.Position(3)*0.8;
            elseif pval <= (0.05 / sum(rangeK)) / 10 / 10 / 10 / 10 /10
                t=title('******', 'Color', darkBlue);
                t.Position(1)=t.Position(1)*0.8;
                t.Position(3)=t.Position(3)*0.8;
            elseif pval <= (0.05 / sum(rangeK)) / 10 / 10 / 10 / 10
                t=title('*****', 'Color', darkBlue);
                t.Position(1)=t.Position(1)*0.8;
                t.Position(3)=t.Position(3)*0.8;
            elseif pval <= (0.05 / sum(rangeK)) / 10 / 10 / 10
                t=title('****', 'Color', darkBlue);
                t.Position(1)=t.Position(1)*0.8;
                t.Position(3)=t.Position(3)*0.8;
            elseif pval <= (0.05 / sum(rangeK)) / 10 / 10
                t=title('***', 'Color', darkBlue);
                t.Position(1)=t.Position(1)*0.8;
                t.Position(3)=t.Position(3)*0.8;
            elseif pval <= (0.05 / sum(rangeK)) / 10
                t=title('**', 'Color', darkBlue);
                t.Position(1)=t.Position(1)*0.8;
                t.Position(3)=t.Position(3)*0.8;
            elseif pval <= (0.05 / sum(rangeK)) % survive Bonferroni
                t=title('*', 'Color', darkBlue);
                t.Position(1)=t.Position(1)*0.8;
                t.Position(3)=t.Position(3)*0.8;
            elseif pval < 0.05 / rangeK(k)
                t=title('o', 'Color', darkGreen);
                t.Position(1)=t.Position(1)*0.9;
                t.Position(3)=t.Position(3)*0.9;
            elseif pval < 0.05
                t=title('+', 'Color', darkRed);
                t.Position(1)=t.Position(1)*0.9;
                t.Position(3)=t.Position(3)*0.9;
            end

        end

        % Create a subplot for the with the errorbars
        subplot_tight(length(rangeK), rangeK(end)*2, Centroid*2 + (k - 1) * rangeK(end)*2, 0.01);
         % Initialize containers for occupancy data for each condition.
         P_cond = cell(1, n_Cond);
         mean_P_cond = zeros(1, n_Cond);
         ste = zeros(1, n_Cond);
         for j = 1:n_Cond
             P_cond{j} = P(Index_Conditions == j, rangeK == k, Centroid);
             mean_P_cond(j) = nanmean(P(Index_Conditions == j, rangeK == k, Centroid));
             ste(j) = std(P(Index_Conditions == j, rangeK == k, Centroid)) / sqrt(numel(P(Index_Conditions == j, rangeK == k, Centroid)));
         end
         hold on
         if pval <= (0.05 / sum(rangeK))
            errorbar(mean_P_cond, ste,"_",'MarkerSize',2,'LineStyle', 'none', 'Color', darkBlue, 'CapSize', 4, 'LineWidth', 1);
         elseif pval < 0.05 / rangeK(k)
            errorbar(mean_P_cond, ste,"_",'MarkerSize',2,'LineStyle', 'none', 'Color', darkGreen, 'CapSize', 4, 'LineWidth', 1); 
         elseif pval < 0.05
            errorbar(mean_P_cond, ste,"_",'MarkerSize',2,'LineStyle', 'none', 'Color', darkRed, 'CapSize', 4, 'LineWidth', 1);
         else
            errorbar(mean_P_cond, ste,"_",'MarkerSize',2,'LineStyle', 'none', 'Color', 'k', 'CapSize', 4, 'LineWidth', 1); 
         end

         set(gca, 'color', 'none')
         set(gca, 'XTick', 1:n_Cond, 'XTickLabel', '', 'Fontsize', 14, 'TickLabelInterpreter', 'none')
         set(gca, 'YTick', [])
         xlim([.4 n_Cond+.6])
         hold off, box off
    end
end

%% Save the Rendered Figure
% Save the figure as both PNG and MATLAB FIG file.
saveas(Fig, fullfile(results_dir, ['_' save_name '_' cortex_dir '.png']), 'png');
saveas(Fig, fullfile(results_dir, [save_name '_' cortex_dir '.fig']), 'fig');
disp(['- Plot successfully saved as ' save_name '_' cortex_dir]);
disp(' ');
