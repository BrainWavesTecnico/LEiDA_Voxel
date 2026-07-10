# LEiDA_Voxel

Leading Eigenvector Dynamics Analysis (LEiDA) applied to voxels defined in 3D volumes rather than pre-defined parcels.

This is a MATLAB pipeline to analyze eigenvector dynamics in fMRI data using LEiDA, working directly in voxel space (a full-brain or custom voxel mask) instead of a fixed brain parcellation. It clusters recurrent BOLD phase-coherence patterns into coupling modes, computes their fractional occupancy per scan, harmonizes occupancy across acquisition sites with ComBat, and tests for differences between conditions (with permutation/bootstrap statistics), visualizing results as 3D brain renders and correlating key modes with clinical/cognitive scores.

Used in: Campo et al., *Cognitive reserve linked to network-specific brain-ventricle coupling modes*, 2025.

## Requirements

- MATLAB with the Statistics and Machine Learning Toolbox (`kmeans`, `statset`) and Image Processing Toolbox (`imresize3`, `niftiread`).
- Preprocessed resting-state fMRI data in NIfTI format, aligned to a common MNI template.

## Pipeline overview

The entry point is [`run_LEiDA_Voxel.m`](run_LEiDA_Voxel.m), which documents and drives the full pipeline:

| Step | Script/Function | Purpose |
|---|---|---|
| 0 (optional) | `Mask_Voxels_of_Interest.m` | Build a custom binary voxel mask in MNI space from a set of per-scan brain masks, keeping voxels present in a chosen proportion of scans, at a chosen voxel size. Not needed if using the bundled full-brain mask. |
| 1 | `Get_EigenVectors_VoxelSpace_Server.m` | Load the preprocessed fMRI NIfTI files, resize them to the mask's voxel space, compute the signal phase (Hilbert transform), and extract the leading eigenvector of the phase coherence matrix at every TR, for every scan. |
| 2 | `LEiDA_cluster_VoxelMNI10mm.m` | Cluster all leading eigenvectors into a range of K clusters (coupling modes) with K-means (cosine distance), for `mink:maxk`. |
| 2b | ComBat (inline in `run_LEiDA_Voxel.m`, via `combat/combat.m`) | Compute the fractional occupancy of each mode per scan, then harmonize it across acquisition sites, keeping diagnosis/age/sex/education as covariates of interest. |
| 3 | `LEiDA_stats_Voxel_FracOccup_ComBat.m` | Test each mode's occupancy between conditions (Welch's t-test for independent samples, paired permutation test for paired samples), with permutation p-values and Hedge's effect sizes. |
| 4 | `Plot_FracOccup_stats.m` | Summary plots of the statistical tests: p-values, mean occupancy barplots, effect sizes. |
| 5 | `Plot_ClustVoxelCentroid_Pyramid_RSNs.m` | Render a pyramid of all centroids across K, with optional Yeo RSN color overlay and significance markers. |
| 6 | `Choose_Relevant_Modes.m` | Automatically select the `[k c]` modes that differ most between conditions (significant after multiple-testing correction, minimum effect size), grouping correlated modes. |
| 7 | `Plot_KeyModes_Slices_Stats.m` | Render each key mode on anatomical slices with mean ± SE occupancy bars per condition. |
| 8 | `Plot_Mode_TransparentBrain.m` | Detailed 3D rendering of each key mode, including overlap with Yeo Resting-State Networks. |
| 9 | `Scores_vs_Mode_Occupancy.m` | Partial correlation (controlling for age) between key mode occupancy and clinical/cognitive scores, plotted and exported to CSV. |

Figures are saved at each step in the results folder in both `.fig` and `.png`/`.jpg`.

## Repository structure

```
run_LEiDA_Voxel.m                    Main pipeline script (documents and runs all steps)
Mask_Voxels_of_Interest.m            Step 0: build a custom voxel mask
Get_EigenVectors_VoxelSpace_Server.m Step 1: extract leading eigenvectors from fMRI data
LEiDA_cluster_VoxelMNI10mm.m         Step 2: K-means clustering
LEiDA_stats_Voxel_FracOccup_ComBat.m Step 3: statistics on mode occupancy
Plot_FracOccup_stats.m               Step 4: statistics summary plots
Plot_ClustVoxelCentroid_Pyramid_RSNs.m  Step 5: centroid pyramid render
Choose_Relevant_Modes.m              Step 6: automatic key-mode selection
Plot_KeyModes_Slices_Stats.m         Step 7: key-mode slice + occupancy plots
Plot_Mode_TransparentBrain.m         Step 8: detailed 3D key-mode render
Scores_vs_Mode_Occupancy.m           Step 9: mode-occupancy vs. scores correlation
combat/                              ComBat site-harmonization toolbox (Johnson et al.)
utilities/                           Colormaps, MNI masks, Yeo RSN parcellation, plotting/stats helpers
```

## Usage

1. Edit the "USER INPUT" section inside `Get_EigenVectors_VoxelSpace_Server.m` (data directory, filename pattern, voxel mask, output filename) and run it to extract leading eigenvectors from your fMRI data.
2. Edit the directories and filenames at the top of `run_LEiDA_Voxel.m` (`fMRI_dir`, `leida_dir`, `Scores_Table`, `results_dir`, `file_V1`, `cluster_file`, `stats_file`) and run the script section by section to cluster, harmonize, test, and visualize the coupling modes.

See the header comments in `run_LEiDA_Voxel.m` for the full function reference, input/output descriptions, and a worked example of each step.

## Notes

- `Scores_vs_Mode_Occupancy.m` uses a fixed set of score-table column indices tailored to the ADNI `Scores_ADNI` table used in Campo et al.; adapt these indices when using a different scores table.
- `Plot_KeyModes_Slices_Stats.m` additionally expects a study-specific `Scores_ADNI_2177scans.mat` file (for sex-based grouping) in the results directory; this file is not included in this repository.

## Author

Joana Cabral, Tecnico, University of Lisbon — joanabcabral@tecnico.ulisboa.pt
