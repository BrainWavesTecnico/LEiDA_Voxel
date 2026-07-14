function Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC)
% Plot_Mode_TransparentBrain provides an in-depth analysis of one or more selected
% modes (centroids) from the clustering solution by rendering each in 3D and
% computing its overlap with Yeo's functional networks.
%
% The function performs the following operations:
%   - Loads clustering results.
%   - Selects and reshapes each chosen centroid (mode) into a 3D volume.
%   - Computes the correlation between the mode and each of 7 RSNs.
%   - Renders the mode in 3D, first as a purple patch on a transparent brain,
%     and then colors the mode in cortex alone using RSN colors.
%
% INPUT:
%   results_dir  - Directory containing the cluster file.
%   cluster_file - Filename with clustering results (centroids, mask, etc.).
%   Key_Modes_KC - Nx2+ matrix with one row per mode, [ki c ...], as returned
%                  by Choose_Relevant_Modes (or built manually). ki is the
%                  POSITION of the clustering solution in rangeK (i.e.
%                  Kmeans_results{ki}), NOT the literal number of clusters -
%                  e.g. if rangeK = 2:20, ki=1 means K=2 clusters, ki=2 means
%                  K=3, etc. c is the mode/centroid index within that solution.
%
% OUTPUT:
%   Displays multiple 3D renderings and saves figures.
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt
% Funded by BPI-LaCaixa Foundation
% and Portuguese Foundation for Science and Technology

%% Load Required Data
% Load clustering centroids and associated mask information.
load([results_dir '/' cluster_file], 'MNI_lowres_Mask', 'ind_voxels','Kmeans_results','rangeK');

N_Modes=size(Key_Modes_KC,1);


%% Prepare RSN Mask and Brain Mask
size_MNI = [91 109 91];
lim_mask = 0.5;
% Load Yeo parcellation in MNI space at 2mm resolution.
V_Yeo = struct2array(load('ParcelsMNI2mm', 'V_Yeo7'));
ind_voxels_Yeo = find(V_Yeo > 0 & V_Yeo < 8);
% Define Yeo network colors as in Yeo et al., 2011.
YeoColor = [180 115 208; 76 112 210; 100 174 117; 230 150 200; 207 225 185; 250 221 107;   227 120 91] ./ 256;

% Load AAL parcels to identify brain areas involved.
load ParcelsMNI2mm.mat V_AAL120 label120
load ParcelsMNI2mm.mat V_Desikan labelDesikan70

% Mode_colors=[137 207 240; 227 120 91; 250 221 107; 207 225 185] ./ 256;

% Create and smooth a brain mask.
Mask_Brain = zeros(size(MNI_lowres_Mask));
Mask_Brain(ind_voxels) = 1;
Mask_Brain = flip(smooth3(Mask_Brain, 'gaussian', 3, 0.8));
Mask_Brain = imresize3(Mask_Brain, size_MNI, 'Method', 'linear');


%% Loop over each Mode

% Use figure handles (not hardcoded numbers) so this always opens new
% figures, even if earlier plotting steps already created Figure 1/2/3.
Fig1 = figure('Position', [521         229        1220         699]);
Fig2 = figure('Position', [ 423     2   801   998]);
Fig3 = figure('Position', [ 423     2   801   998]);
  

for Mode = 1:N_Modes
    % Transform Selected Mode into 3D Volume and Compute Overlap with RSNs
    % Retrieve the centroid vector and reshape it.
    ki=Key_Modes_KC(Mode,1);
    c=Key_Modes_KC(Mode,2);

    Vc = Kmeans_results{ki}.C(c, :);
    % Create an empty 3D volume with the size of the MNI brain with 1cm3 voxels
    Vc_3D = zeros(size(MNI_lowres_Mask));
    % Insert the elements of the eigenvector on the corresponding volxels
    Vc_3D(ind_voxels) = Vc;
    % Smooth before Resizing
    Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.4);
    Vc_3D = imresize3(Vc_3D, size_MNI, 'Method', 'linear');

    % Compute correlations with each Yeo network.
    for net = 1:7
        dice_net(net) = dice(V_Yeo(ind_voxels_Yeo) == net, Vc_3D(ind_voxels_Yeo) > 0);
    end

    % Clear activations outside the brain mask and smooth.
    Vc_3D(Mask_Brain < lim_mask) = 0;
    Vc_3D = smooth3(Vc_3D, 'gaussian', 3, 0.8);

    %%%% Left - Render Mode as 3D Red Patch on Transparent Brain
    figure(Fig1)
    subplot_tight(N_Modes, 6, 1 +(Mode-1)*6, 0.02)
    hold on
    % Plot a transparent brain patch.
    Brainpatch = patch(isosurface(Mask_Brain, 0.3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    reducepatch(Brainpatch, 0.4);
    isonormals(Mask_Brain, Brainpatch);
    % Overlay the mode as a red patch.
    Redpatch = patch(isosurface(smooth3(Vc_3D > 0), 0.5), 'FaceColor', [1 .2 .2], 'EdgeColor', 'none');
    reducepatch(Redpatch, 0.5);
    isonormals(Vc_3D, Redpatch);
    material dull; lighting flat; daspect([1 1 1]);
    Az = 180; El = 0;
    view(Az, El) % Sagittal view
    camlight; xlim([10 105]); ylim([5 85]); zlim([2 85]); axis off

    % Additional View of the Mode
    % Azimuth and Elevation of camera


    subplot_tight(N_Modes, 6, 2 +(Mode-1)*6, 0.02)
    hold on
    % Plot a transparent brain patch.
    Brainpatch = patch(isosurface(Mask_Brain, 0.3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    reducepatch(Brainpatch, 0.4);
    isonormals(Mask_Brain, Brainpatch);
    % Overlay the mode as a red patch.
    Redpatch = patch(isosurface(smooth3(Vc_3D > 0), 0.5), 'FaceColor', [1 .2 .2], 'EdgeColor', 'none');
    reducepatch(Redpatch, 0.5);
    isonormals(Vc_3D, Redpatch);
    material dull; lighting flat; daspect([1 1 1]);
    Az = -90; El = 90;
    view(Az, El) % Sagittal view
    camlight; xlim([10 105]); ylim([5 85]); zlim([2 85]); axis off

    if Mode>1
        % Second Rendering: Mode Overlap with RSNs
        subplot_tight(N_Modes, 6, 3 +(Mode-1)*6 , 0.02)
        hold on
        % Plot a transparent brain with a grey tone.
        Brainpatch = patch(isosurface(Mask_Brain, 0.3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        reducepatch(Brainpatch, 0.4);
        isonormals(Mask_Brain, Brainpatch);
        % Create a patch using only voxels overlapping with Yeo parcellation.
        Patch_cortical = zeros(size(Vc_3D));
        Patch_cortical(ind_voxels_Yeo) = Vc_3D(ind_voxels_Yeo);
        if max(dice_net > 0.4)
            Mode_patch = patch(isosurface(smooth3(Patch_cortical > 0), 0.4), 'FaceColor', YeoColor(dice_net == max(dice_net), :), 'EdgeColor', 'none');
        else
            Mode_patch = patch(isosurface(smooth3(Patch_cortical > 0), 0.4), 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
        end
        reducepatch(Mode_patch, 0.5);
        isonormals(Patch_cortical, Mode_patch);
        material dull; lighting flat; daspect([1 1 1]);
        Az = 180; El = 0;
        view(Az, El);     camlight; xlim([10 105]); ylim([5 85]); zlim([2 85]); axis off

        subplot_tight(N_Modes, 6, 4 +(Mode-1)*6 , 0.02)
        hold on
        % Plot a transparent brain with a grey tone.
        Brainpatch = patch(isosurface(Mask_Brain, 0.3), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        reducepatch(Brainpatch, 0.4);
        isonormals(Mask_Brain, Brainpatch);
        % Create a patch using only voxels overlapping with Yeo parcellation.
        Patch_cortical = zeros(size(Vc_3D));
        Patch_cortical(ind_voxels_Yeo) = Vc_3D(ind_voxels_Yeo);
        if max(dice_net > 0.4)
            Mode_patch = patch(isosurface(smooth3(Patch_cortical > 0), 0.4), 'FaceColor', YeoColor(dice_net == max(dice_net), :), 'EdgeColor', 'none');
        else
            Mode_patch = patch(isosurface(smooth3(Patch_cortical > 0), 0.4), 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
        end
        reducepatch(Mode_patch, 0.5);
        isonormals(Patch_cortical, Mode_patch);
        material dull; lighting flat; daspect([1 1 1]);
        Az = -90; El = 90;
        view(Az, El);     camlight; xlim([10 105]); ylim([5 85]); zlim([2 85]); axis off


        % Bar Plot of the Overlap with RSNs
        subplot(4, 6, 6 +(Mode-1)*6)
        b = bar(1, dice_net);
        for i = 1:7
            set(b(i), 'FaceColor', YeoColor(i, :), 'BarWidth', 0.6);
        end
        xticks([]); xlabel('RSNs', 'FontSize', 8);
        ylabel({'Dice index'}, 'FontSize', 8);

        if Mode==2
            legend({'Visual','Somatomotor','Dorsal Attention','Ventral Attention','Limbic','Frontoparietal','Default Mode'},...
                'Location', 'northwest', 'FontSize', 8);
        end
        box off;
        ylim([0 0.75]);
        set(gca,"YGrid","on")


        Vc_AAL120 = zeros(1, max(V_AAL120(:)));
        Vc_3D = Vc_3D(:);
        for v = 1:max(V_AAL120(:))
            ind_area = V_AAL120 == v;
            Vc_AAL120(v) = mean(Vc_3D(ind_area) > 0);
        end
        [Vc_AAL_reoder, order_indices] = sort(Vc_AAL120, 'descend');
        label120_reorder = label120(order_indices, :);

        figure(Fig2)
        subplot(1,N_Modes-1,Mode-1)
        barh(Vc_AAL_reoder(end:-1:1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none')
        yticks(1:max(V_AAL120(:)))
        yticklabels(label120_reorder(end:-1:1, :))

        set(gca, 'Fontsize', 8)
        ylim([0 size(Vc_AAL_reoder, 2)+.5])
        title(['BraVe Mode ' num2str(Mode)], 'FontSize', 10)
        set(gcf, 'Position', [620 65 372 982])
        xlabel({'Proportion of ','Phase-shifted voxels'})


        Vc_Desikan = zeros(1, max(V_Desikan(:)));
        Vc_3D = Vc_3D(:);
        for v = 1:max(V_Desikan(:))
            ind_area = V_Desikan == v;
            Vc_Desikan(v) = mean(Vc_3D(ind_area) > 0);
        end
        [Vc_Desikan_reoder, order_indices] = sort(Vc_Desikan, 'descend');
        labelDesikan70_reorder = labelDesikan70(order_indices, :);

        figure(Fig3)
        subplot(1, N_Modes-1, Mode-1)
        barh(Vc_Desikan_reoder(end:-1:1), 'FaceColor', [0.3 0.3 0.3], 'EdgeColor', 'none')
        yticks(1:max(V_Desikan(:)))
        yticklabels(labelDesikan70_reorder(end:-1:1))
        ax = gca;
        ax.TickLabelInterpreter = 'none';

        set(gca, 'Fontsize', 8)
        ylim([0 size(labelDesikan70_reorder, 1)+.5])
        title(['BraVe Mode ' num2str(Mode)], 'FontSize', 10)
        set(gcf, 'Position', [620 65 372 982])
        xlabel({'Proportion of ','Phase-shifted voxels'})




    end


end

%%

% % Save the figures.
% figure(1)
% saveas(gcf, [results_dir 'ModeTranspOverlap' '_K' num2str(k) '_c' num2str(c)], 'fig');
% saveas(gcf, [results_dir '_Fig3_ModeTranspOverlap_K' num2str(k) '_c' num2str(c)], 'png');
% 
% figure(2)
% saveas(gcf, [results_dir 'ListAreas_Key_Modes_AAL2'], 'fig')
% 
% figure(3)
% saveas(gcf, [results_dir 'ListAreas_Key_Modes_Desikan'], 'fig')
% 
% %% Plot list of the brain areas located on dephased poles (using AAL2 atlas)
% 
% 
% %% Save the Figure
% % Save the figure as both PNG and MATLAB FIG file.
% saveas(gcf, [results_dir '_List_Brain_areas_K' num2str(k) '_c' num2str(c) '.png'], 'png');
% saveas(gcf, [results_dir 'List_Brain_areas_K' num2str(k) '_c' num2str(c) '.fig'], 'fig');
% disp('- Plot successfully saved as List_Brain_areas');
% disp(' ');