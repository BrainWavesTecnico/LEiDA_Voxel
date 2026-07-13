% Get_EigenVectors_VoxelSpace_Server
%
% Function to save the leading eigenvectors from the fMRI data.
% Vectors are concatenated acorss all scans ready for clustering.
%
% This function loads the fMRI signals aligned to a brain template.
% Resizes the template to the spatial resolution of the brain mask. 
% To have a reasonable number of voxels of voxels in the brain mask, we use voxels with 10mm3 size, resulting in 1821 voxels of interest.
% Gets the Hilbert transform from signals in the voxels of interest in the brain mask.
% Computes the fMRI phase leading eigenvectors for each TR for all participants. 
%
% INPUT:
% data_dir      directory where the fMRI data in NIFTI format is saved
% leida_dir     directory where the results from running LEiDA are saved
%
% OUTPUT:
% V1_MNI10mm.mat    leading eigenvectors in MNI 10mm space
%
% Author: Joana Cabral, University of Lisbon, joanabcabral@tecnico.ulisboa.pt

% USER INPUT 
% Main directory where the  fMRI data preprocessed in MNI space is stored as Nifti files (.nii or .nii.gz are stored)
fMRI_dir = '/pipeline_cpac-adni-ants-mean-best-for-csf_0p01_0p1/';
% Part of the file name common to all fMRI scans that we will load
extension_name= 'rest_space-MNI152NLin6ASym_desc-preproc_bold.nii.gz';
% Directory where the results from this run will be stored
results_dir    = fMRI_dir;

% File containing a binary mask of brain voxels of interest in MNI space
% Adjust voxel size to tune the number of Voxels according to RAM space
Mask_file='MNI_10mm3_FullBrain.mat';

% File to save with the leading eigenvectors
file_V1     = 'LEiDA_V1_all_MNI10mm_FullMask.mat';    

% Set here the maximum number of volumes in a scan to allocate space for the entire matrix of eigenvectors
TimeMax=200;

%%

% Add all subfolders of results directory to path
addpath(genpath(results_dir));
addpath(genpath(fMRI_dir))

% Get number of files in folder
data_info = dir([fMRI_dir '/**/*' extension_name]);
num_scans = numel(data_info);
disp(['Total scans = ' num2str(num_scans)])

% Load 10mm MNI voxel space
load(Mask_file,'MNI_lowres_Mask');
sz = size(MNI_lowres_Mask); % size of the 10mm MNI mask
ind_voxels = find(MNI_lowres_Mask(:) > 0); % find the non-zero elements in the mask
n_voxels = length(ind_voxels);

disp(['Number of brain voxels N=' num2str(n_voxels)])

% Matrix to store the leading eigenvectors of all subjects at each TR
V1_all = zeros(num_scans*TimeMax,n_voxels);
Scan_num=zeros(num_scans*TimeMax,1);
Scan_length=zeros(num_scans,1);

t_all = 0;
for s = 1:num_scans
        disp(['Scan ' num2str(s) ' ' data_info(s).name ':']);
        % Read the nii file
        fMRI_MNI = niftiread([data_info(s).folder '/' data_info(s).name]);
        T = size(fMRI_MNI,4); % number of volumes
        Scan_length(s)=T;

        disp('- Resizing NIFTI file to Mask size space');
        % Files will be resized in order to be accomodated to the MNI10mm template
        fMRI_MNI_lowres = zeros(sz(1), sz(2), sz(3), T);
        for t = 1:T
            fMRI_volume=squeeze(fMRI_MNI(:,:,:,t));
            fMRI_MNI_lowres(:,:,:,t) = imresize3(fMRI_volume,sz,'method','linear');
        end
        clear fMRI_MNI

        disp('- Computing the signal phases using the Hilbert transform');
        % Store the fMRI signal phase using the Hilbert transform
        fMRI_ts = zeros(n_voxels,T);

        for v = 1:n_voxels
            [I1,I2,I3] = ind2sub(sz,ind_voxels(v));
            fMRI_ts(v,:) = squeeze(fMRI_MNI_lowres(I1,I2,I3,:))';
        end
        clear fMRI_MNI_lowres

        disp('- Computing the leading eigenvectors at each time point');

        % De-mean and get the signal phase using the Hilbert transform
        for v = 1:n_voxels
            fMRI_ts(v,:) = angle(hilbert(detrend(fMRI_ts(v,:) - mean(fMRI_ts(v,:)))));
        end

        % Compute leading eigenvector of each phase coherence matrix
        for t = 2:T-1 % exclude 1st and last TR after hilbert transform

            % Save the leading eigenvector for time t
            [v1,~] = eigs(cos(fMRI_ts(:,t) - fMRI_ts(:,t)'),1);

            t_all = t_all + 1; % time point in V1_all

            % row t_all correponds to the computed eigenvector at time t for scan s
            V1_all(t_all,:) = v1;

            Scan_num(t_all)=s;
        end

end
% Reduce size in case some scans have less TRs than tmax
V1_all(t_all+1:end,:) = [];
Scan_num(t_all+1:end)=[];
V1_all(sum(V1_all,2)>0,:)= -(V1_all(sum(V1_all,2)>0,:));


disp(' ');
disp('- Saving leading eigenvectors in MNI 10mm space');
save([results_dir  file_V1], 'V1_all', 'ind_voxels', 'MNI_lowres_Mask', 'data_info', 'Scan_num','Scan_length', '-v7.3')

