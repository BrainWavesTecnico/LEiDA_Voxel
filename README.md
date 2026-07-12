# LEiDA_Voxel

Leading Eigenvector Dynamics Analysis (LEiDA) applied to voxels defined in 3D volumes rather than pre-defined parcels.

This is a MATLAB pipeline to analyze eigenvector dynamics in fMRI data using LEiDA, working directly in voxel space (a full-brain or custom voxel mask) instead of a fixed brain parcellation. It clusters recurrent BOLD phase-coherence patterns into coupling modes, computes their fractional occupancy per scan, harmonizes occupancy across acquisition sites with ComBat, and tests for differences between conditions (with permutation/bootstrap statistics), visualizing results as 3D brain renders and correlating key modes with clinical/cognitive scores.

Used in: *Cognitive function linked to temporal occupancy of Brain-Ventricle (BraVe) modes*. Campo, Miguel, Brattico, Nigro, Tafuri, Logroscino, Cabral and the Alzheimer's Disease Neuroimaging Initiative (ADNI). bioRxiv 2025.01.04.631289; doi: [https://doi.org/10.1101/2025.01.04.631289](https://doi.org/10.1101/2025.01.04.631289).

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
| 2b | `Save_Occupancies_Harmonize.m` (via `combat/combat.m`) | Compute the fractional occupancy of each mode per scan, and, optionally (`apply_combat`), harmonize it across acquisition sites with ComBat, keeping diagnosis/age/sex/education as covariates of interest. Saves both `P_original` and `P_harmonized`. |
| 3 (optional) | `LEiDA_stats_Voxel_FracOccup_ComBat.m` | Test each mode's occupancy between conditions (Welch's t-test for independent samples, paired permutation test for paired samples), with permutation p-values and Hedge's effect sizes. Runs on either `P_original` or `P_harmonized`, chosen when calling it. Skip for studies with no discrete conditions to compare. |
| 4 | `Plot_FracOccup_stats.m` | Summary plots of the statistical tests: p-values, mean occupancy barplots, effect sizes. Requires step 3. |
| 5 | `Plot_ClustVoxelCentroid_Pyramid_RSNs.m` | Render a pyramid of all centroids across K, with optional Yeo RSN color overlay and significance markers. Requires step 3. |
| 6 | `Choose_Relevant_Modes.m` | Automatically select the `[ki c]` modes that differ most between conditions (significant after multiple-testing correction, minimum effect size), grouping correlated modes. Requires step 3; without it, specify `Key_Modes_KC` manually instead. |
| 7 | `Plot_KeyModes_Slices_Stats.m` | Render each key mode on anatomical slices with mean ± SE occupancy bars per condition. Requires step 3. |
| 8 | `Plot_Mode_TransparentBrain.m` | Detailed 3D rendering of each key mode, including overlap with Yeo Resting-State Networks. Only needs `cluster_file` and `Key_Modes_KC` — does not require step 3. |
| 9 | `Scores_vs_Mode_Occupancy.m` | Partial correlation (controlling for age) between key mode occupancy and clinical/cognitive scores, plotted and exported to CSV. Takes the occupancy matrix `P` directly (from step 2b) and `Key_Modes_KC` — does not require step 3, so it works standalone for studies with only continuous scores. |

Figures are saved at each step in the results folder in both `.fig` and `.png`/`.jpg`.

## Repository structure

```
run_LEiDA_Voxel.m                    Main pipeline script (documents and runs all steps)
run_LEiDA_vox.m                      Code Ocean capsule entry point (steps 2-4, from pre-extracted eigenvectors)
Select_Demo_Subsample.m              One-off local script: build a balanced demo subsample for the capsule's data/
Mask_Voxels_of_Interest.m            Step 0: build a custom voxel mask
Get_EigenVectors_VoxelSpace_Server.m Step 1: extract leading eigenvectors from fMRI data
LEiDA_cluster_VoxelMNI10mm.m         Step 2: K-means clustering
Save_Occupancies_Harmonize.m         Step 2b: mode occupancy extraction + optional ComBat harmonization
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

## Code Ocean capsule

[`run_LEiDA_vox.m`](run_LEiDA_vox.m) is a batch entry point (a function, not a script with hardcoded personal paths) for reproducing the pipeline as a Code Ocean capsule, starting from already-extracted leading eigenvectors (i.e. skipping step 0/1, which need the raw fMRI NIfTI files):

1. On your own machine, with the full cohort's already-extracted eigenvectors and `Scores_ADNI` table, run [`Select_Demo_Subsample.m`](Select_Demo_Subsample.m) to pick a small, balanced, unique-participant demo subsample (e.g. 30 scans per condition) and save the two demo files.
2. In the capsule, put those two files in `data/` (named to match `file_V1`/`Scores_Table` at the top of `run_LEiDA_vox.m`, or edit those two lines to match your filenames).
3. Put all the `.m` files (including `combat/` and `utilities/`) in `code/`.
4. `results/` is written to automatically; call `run_LEiDA_vox()` with no arguments and it clusters K = 2:20, extracts occupancies, runs the condition statistics, and saves every figure there.

ComBat harmonization is off by default in the capsule (`apply_combat = 0` inside `run_LEiDA_vox.m`) since a small demo subsample is unlikely to have enough scans per site for reliable harmonization — turn it on if your demo data spans multiple well-populated sites. If no mode survives the significance threshold on the reduced sample (likely with far fewer scans than the full study), the driver falls back to a fixed mid-K mode selection so the figure-generation steps still produce output, and logs a warning explaining why.

## Notes

- **`Key_Modes_KC`'s first column (`ki`) is a POSITION into `rangeK`, not the literal number of clusters.** `Kmeans_results{ki}` and `P(:,ki,:)` are indexed by position (e.g. if `rangeK = 2:20`, `ki=1` means K=2 clusters, `ki=2` means K=3). `Choose_Relevant_Modes.m` returns `ki` correctly; when building `Key_Modes_KC` manually, convert with `ki = find(rangeK == K_wanted)` rather than using the literal K value directly — this was a real bug in earlier versions of `Plot_KeyModes_Slices_Stats.m`/`Plot_Mode_TransparentBrain.m`/`Choose_Relevant_Modes.m` (which only happened to work when `mink=1`), now fixed to use `ki` consistently everywhere.
- **Studies with only continuous scores, no discrete conditions**: skip step 3 (`LEiDA_stats_Voxel_FracOccup_ComBat`) and steps 4-7 that depend on it. Run steps 0-2b as usual, pick `P = P_original` or `P = P_harmonized`, choose `Key_Modes_KC` manually (e.g. `Key_Modes_KC = [2 4; 4 6]`, using the `ki` convention above), then call `Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC)` and `Scores_vs_Mode_Occupancy(P, Scores_Table, Key_Modes_KC, results_dir, save_name)` directly — neither depends on step 3's output.
- `Save_Occupancies_Harmonize.m` saves both the raw (`P_original`) and ComBat-harmonized (`P_harmonized`) occupancies; `run_LEiDA_Voxel.m` exposes a `use_harmonized_occupancies` toggle to pick which one `P` refers to for the rest of the pipeline (both step 3 and `Scores_vs_Mode_Occupancy.m`).
- `Scores_vs_Mode_Occupancy.m` uses a fixed set of score-table column indices tailored to the ADNI `Scores_ADNI` table used in Campo et al.; adapt these indices when using a different scores table.
- `Plot_KeyModes_Slices_Stats.m` additionally takes `Scores_Table` (used for sex-based grouping via `PTGENDER`).

## Author

Joana Cabral, University of Lisbon — joanabcabral@tecnico.ulisboa.pt
