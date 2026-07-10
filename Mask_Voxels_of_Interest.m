
% Script to Define a Mask with the Voxels of Interest for use in LEiDA
% (instead of using a parcelation scheme)

% In our case the CPAC preprocessing pipeline outputs a mask used for brain extraction
% where we kept signal in CSF compartments


% Find all files with _mask extension in the folder with preprocessed data
data_path='/Users/user/Documents/Research/CognitiveDecline/bold_masks_200/bold_masks';
extension_name='-bold_mask';
file_names = dir(fullfile(data_path, '**', ['*' extension_name '*.nii.gz']));

% Since not all scans include the same voxels, and since a few scans may
% miss relevant voxels, define minimum proportion of masks with that voxel
% to include that voxel in the final mask.

prop_masks=0.95; % 

% Resize the entire volume to tune the number of voxels to include
resize_scale=0.2; % Here going from 2 mm to 10 mm (2/0.2 =10)

n_scans=length(file_names);

fmri_mask=niftiread([file_names(1).folder '/' file_names(1).name]);

masks_total=zeros(n_scans,size(fmri_mask,1),size(fmri_mask,2),size(fmri_mask,3));
for scan=1:n_scans

    fmri_mask=niftiread([file_names(scan).folder '/' file_names(scan).name]);
    masks_total(scan,:,:,:)=fmri_mask;
end

Sum_Mask=squeeze(sum(masks_total));

figure
for s=1:91
subplot_tight(7,13,s)
        imagesc(squeeze(Sum_Mask(:,:,s)),[0 n_scans])
        axis image
        axis off
end


%% Visualize the mask

% 
figure
for s=1:91
subplot_tight(7,13,s)
        imagesc(squeeze(Sum_Mask(:,:,s))>=round(n_scans*prop_masks),[0 1])
        axis image
        axis off
end

Resize_mask=imresize3(Sum_Mask>=round(n_scans*prop_masks),resize_scale,'Method','linear');

figure
for s=1:size(Resize_mask,3)
        subplot_tight(4,5,s)
        imagesc(squeeze(Resize_mask(:,:,s)),[0 1])
        axis image
        axis off
end

% 

save(['Mask_10mm_n' num2str(n_scans) '_95thr.mat'],'masks_total','Sum_Mask','Resize_mask')
