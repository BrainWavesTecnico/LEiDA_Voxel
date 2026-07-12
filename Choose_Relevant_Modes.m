function [Key_Modes_KC, Key_Centroids] = Choose_Relevant_Modes(results_dir, cluster_file, stats_file)
% Choose_Relevant_Modes automatically selects the modes that differ most
% between conditions after LEiDA_stats_Voxel_FracOccup_ComBat, keeping modes
% that are significant (after multiple-testing correction) with an effect
% size > 0.35, then grouping strongly-correlated modes (corr > 0.65) and
% keeping one representative (most significant) mode per group.
%
% INPUT:
%   results_dir  - Directory where the cluster and stats files are stored.
%   cluster_file - Clustering results file name.
%   stats_file   - Statistical results file name (output of
%                  LEiDA_stats_Voxel_FracOccup_ComBat).
%
% OUTPUT:
%   Key_Modes_KC   - Nx12 matrix, one row per selected key mode:
%                    [ki, c, slope, n_sig_pairs, pval(1:3), effectsize(1:3), group, group_order].
%                    ki is the POSITION of the clustering solution in rangeK
%                    (i.e. Kmeans_results{ki}), NOT the literal number of
%                    clusters - e.g. if rangeK = 2:20, ki=1 means K=2
%                    clusters, ki=2 means K=3, etc. c is the mode/centroid
%                    index within that solution. slope (column 3) is the mean
%                    occupancy change between the last and first condition.
%   Key_Centroids  - Corresponding centroid vectors, one row per key mode.
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

% Function to detect the most relevant modes after statistical analysis

generate_pyramid_groups=0;

% Load statistical occupancy data and condition info.
load([results_dir stats_file], 'P', 'P_pval', 'cond', 'condCol', 'condRow', 'Index_Conditions','rangeK','effectsize');

P_pval_sig_sum=squeeze(sum(P_pval< (0.05 / sum(rangeK)/ 3/2) & P_pval>0 & effectsize> 0.35));

% NOTE: ind2sub must use P_pval_sig_sum's actual size (length(rangeK) x rangeK(end)).
% A previous version hardcoded [20, 20], which only matched when rangeK == 1:20;
% for any other mink/maxk it silently produced wrong (K,c) subscripts.
[K, c]=ind2sub(size(P_pval_sig_sum),find(P_pval_sig_sum));

Signif_Modes_KC=[K, c];

cortex_dir='SideView';

% Load clustering centroids from the clusters file.
load([results_dir cluster_file], 'Kmeans_results')

Centroids_Signif=zeros(size(Signif_Modes_KC,1),size(Kmeans_results{1}.C,2));

for Mode=1:size(Signif_Modes_KC,1)
    Centroids_Signif(Mode,:)=Kmeans_results{Signif_Modes_KC(Mode,1)}.C(Signif_Modes_KC(Mode,2),:);
    mean_P_cond = zeros(1, size(cond,2));
    for j = 1:size(cond,2)
        mean_P_cond(j) = nanmean(P(Index_Conditions == j, Signif_Modes_KC(Mode,1), Signif_Modes_KC(Mode,2)));
    end
    Signif_Modes_KC(Mode,3)=mean_P_cond(end)-mean(mean_P_cond(1));
end

Signif_Modes_KC=cat(2,Signif_Modes_KC,zeros(size(Signif_Modes_KC,1),8));

for Mode=1:size(Signif_Modes_KC,1)
    k=Signif_Modes_KC(Mode,1);
    c=Signif_Modes_KC(Mode,2);
    Signif_Modes_KC(Mode,4)=P_pval_sig_sum(k,c);
    Signif_Modes_KC(Mode,5:7)=P_pval(:,k,c);
    Signif_Modes_KC(Mode,8:10)=effectsize(:,k,c);

end

Centroids_Signif_reorder=zeros(size(Centroids_Signif));
Signif_Modes_KC_reorder=zeros(size(Signif_Modes_KC));

G=1;
ind_start=1;
while size(Centroids_Signif,1)>0

    [line_index, ~]=ind2sub(size(Signif_Modes_KC),find(Signif_Modes_KC==min(min(Signif_Modes_KC(:,5:7)))));
    %[line_index, ~]=ind2sub(size(Signif_Modes_KC),find(squeeze(max(effectsize))==max(max(squeeze(effectsize(2,:,:))))));

    CorrMode=corrcoef(Centroids_Signif');
    Group=find(CorrMode(line_index,:)>0.65);
    Centroids_Signif_reorder(ind_start:ind_start+length(Group)-1,:)=Centroids_Signif(Group,:);
    Signif_Modes_KC_reorder(ind_start:ind_start+length(Group)-1,:)=Signif_Modes_KC(Group,:);
    Signif_Modes_KC_reorder(ind_start:ind_start+length(Group)-1,11)=G;

    Centroids_Signif(Group,:)=[];
    Signif_Modes_KC(Group,:)=[];

    ind_start=ind_start+length(Group);
    G=G+1;
end
N_groups=G-1; %




Change_Slope=zeros(1,N_groups);
% Order groups by most decrease to most increase
for g=1:N_groups
    Change_Slope(g)=mean(Signif_Modes_KC_reorder(Signif_Modes_KC_reorder(:,11)==g,3));
end
[~, Group_reorder]=sort(Change_Slope, 'ascend');

for Mode=1:size(Signif_Modes_KC_reorder,1)
    Signif_Modes_KC_reorder(Mode,12)=find(Group_reorder==Signif_Modes_KC_reorder(Mode,11));
end



Key_member_index=zeros(N_groups,1);

for G=1:N_groups

    Group_members=Signif_Modes_KC_reorder(:,12)==G;

        %Group_members=Group_members(Signif_Modes_KC_reorder(Group_members,4) == max(Signif_Modes_KC_reorder(Group_members,4)));
        [Key_member_index(G), b]=ind2sub(size(Signif_Modes_KC_reorder),find(Signif_Modes_KC_reorder==min(min(Signif_Modes_KC_reorder(Group_members,5:7)))));

end

[~, b]=sort(Signif_Modes_KC_reorder(:,12));

Centroids_Signif_reorder=Centroids_Signif_reorder(b,:);

% % Figure correlation between modes
% figure
% colormap(jet)
% imagesc(corrcoef(Centroids_Signif_reorder'), [0 1]);
% 
% labels = arrayfun(@(a,b) sprintf('%d %d', a, b), ...
%     round(Signif_Modes_KC_reorder(b,1)), ...
%     round(Signif_Modes_KC_reorder(b,2)), ...
%     'UniformOutput', false);
% 
% set(gca, "XTick", 1:size(Signif_Modes_KC_reorder, 1),'FontSize',6)
% set(gca, "XTickLabel", labels)
% set(gca, "YTick", 1:size(Signif_Modes_KC_reorder, 1),'FontSize',6)
% set(gca, "YTickLabel", labels)
% colorbar
% axis square




Key_Centroids=Centroids_Signif_reorder(Key_member_index,:);
Key_Modes_KC=Signif_Modes_KC_reorder(Key_member_index,:);

%%

if generate_pyramid_groups
    % Load clustering centroids and masks from the clusters file.
    load([results_dir cluster_file], 'MNI_lowres_Mask', 'ind_voxels');

    V_Yeo = struct2array(load('ParcelsMNI2mm', 'V_Yeo7'));
    size_MNI = size(V_Yeo);

    % Create a smoothed, resampled brain mask for the transparent brain
    lim_mask = 0.5;
    Mask_Brain = zeros(size(MNI_lowres_Mask));
    Mask_Brain(ind_voxels) = 1;
    Mask_Brain = flip(smooth3(Mask_Brain, 'gaussian', 3, 0.8));
    Mask_Brain = imresize3(Mask_Brain, size_MNI, 'Method', 'linear');

    %Render and Plot Clusters
    % Create a full-screen figure for rendering.
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT CENTROIDS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    disp(' ');

    Group_cmap=[137 207 240; 227 120 91; 250 221 107; 207 225 185] ./ 256;

    disp('Rendering the centroids colored according to group assignement:')
    Fig = figure('Position', get(0, 'Screensize'));

    % Loop over each clustering solution based on different K values.
    for k = 1: length(rangeK)
        disp(['- K = ' num2str(rangeK(k))])
        % Loop over each state within the current clustering solution.
        for Cluster = 1:rangeK(k)

            % Get the Cluster order from the reordering performed earlier.
            %Cluster = c_reorder{k}(c);

            % Retrieve the Cluster's vector and reshape it into 3D volume using voxel indices.
            Vc = Kmeans_results{k}.C(Cluster, :);
            Vc_3D = zeros(size(MNI_lowres_Mask));
            Vc_3D(ind_voxels) = Vc;
            % Smooth and resize the Cluster volume to match the MNI mask.
            Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.4);
            Vc_3D = imresize3(Vc_3D, size_MNI, 'Method', 'linear');


            % Threshold and further smooth the Cluster volume by masking out
            % values outside the transparent brain.
            Vc_3D(Mask_Brain < lim_mask) = 0;
            Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.8);

            % Create a subplot for the current cluster state.
            subplot_tight(length(rangeK), rangeK(end), Cluster + (k - 1) * rangeK(end), 0.01);
            hold on
            % Plot a transparent brain surface for context.
            Brainpatch = patch(isosurface(Mask_Brain, 0.3), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.2);

            reducepatch(Brainpatch, 0.4);
            isonormals(Mask_Brain, Brainpatch);

            % Plot the Cluster patch.
            if P_pval_sig_sum(k,Cluster)>0
                group_assignement=Signif_Modes_KC_reorder((Signif_Modes_KC_reorder(:,1)==k & Signif_Modes_KC_reorder(:,2)==Cluster),12);
                % Color based on Yeo network color corresponding to maximum overlap.
                Mode_patch = patch(isosurface(smooth3(Vc_3D > 0), 0.6), 'FaceColor', Group_cmap(group_assignement, :), 'EdgeColor', 'none');
            else
                % Render in grey if not significant
                Mode_patch = patch(isosurface(smooth3(Vc_3D > 0), 0.6), 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
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

            if Key_Modes_KC(:,1)==k & Key_Modes_KC(:,2)==Cluster
                axis on
                title('*','FontColor', Group_cmap(group_assignement, :));
            end


        end
    end

    %% Save the Rendered Figure
% Save the figure as both PNG and MATLAB FIG file.
saveas(Fig, fullfile(results_dir, ['GroupPyramid_' cortex_dir '.png']), 'png');
saveas(Fig, fullfile(results_dir, ['Group_Pyramid_' cortex_dir '.fig']), 'fig');
disp('- Plot successfully saved ');
disp(' ');

end

