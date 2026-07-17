# LEiDA_Voxel

Leading Eigenvector Dynamics Analysis (LEiDA) applied to voxels defined in 3D volumes rather than pre-defined parcels.

This is a MATLAB pipeline to analyze eigenvector dynamics in fMRI data using LEiDA, working directly in voxel space (a full-brain or custom voxel mask) instead of a fixed brain parcellation. It clusters recurrent phase-coherence patterns into coupling modes, computes their fractional occupancy per scan, harmonizes occupancy across acquisition sites with ComBat, and tests for differences between conditions (with permutation/bootstrap statistics), visualizing results as 3D brain renders and correlating key modes with clinical/cognitive scores.

Used in: *Cognitive function linked to temporal occupancy of Brain-Ventricle (BraVe) modes*. Campo, Miguel, Brattico, Nigro, Tafuri, Logroscino, Cabral and the Alzheimer's Disease Neuroimaging Initiative (ADNI). bioRxiv 2025.01.04.631289; doi: [https://doi.org/10.1101/2025.01.04.631289](https://doi.org/10.1101/2025.01.04.631289).

## Requirements

- MATLAB, tested on R2025b (should also work on other recent releases). Required toolboxes: Statistics and Machine Learning Toolbox (`kmeans`, `statset`, `partialcorr`) and Image Processing Toolbox (`imresize3`, `niftiread`). Parallel Computing Toolbox is recommended (not required) — clustering requests parallel K-means and falls back to serial execution without it.
- No non-standard hardware required; a standard multi-core desktop/laptop CPU is sufficient.
- Preprocessed resting-state fMRI data in NIfTI format, aligned to a common MNI template. See [`Preprocessing/`](Preprocessing/) for the C-PAC pipeline used to preprocess the ADNI data for this study.
- Expected runtime: ~2 hours for the [Code Ocean capsule](CodeOcean_Capsule/)'s 300-scan demo sample; ~24 hours for the full 2177-scan manuscript dataset (K=2:20, all pipeline steps), on a normal desktop.

## Pipeline overview

The entry point is [`run_LEiDA_Voxel.m`](run_LEiDA_Voxel.m), which documents and drives the full pipeline:

| Step | Script/Function | Purpose |
|---|---|---|
| 0 (optional) | `Mask_Voxels_of_Interest.m` | Build a custom binary voxel mask in MNI space from a set of per-scan brain masks, keeping voxels present in a chosen proportion of scans, at a chosen voxel size. Not needed if using the bundled full-brain mask. |
| 1 | `Get_EigenVectors_VoxelSpace_Server.m` | Load the preprocessed fMRI NIfTI files, resize them to the mask's voxel space, compute the signal phase (Hilbert transform), and extract the leading eigenvector of the phase coherence matrix at every TR, for every scan. |
| 2 | `LEiDA_cluster_VoxelMNI10mm.m` | Cluster all leading eigenvectors into a range of K clusters (coupling modes) with K-means (cosine distance), for `mink:maxk`. |
| 2b | `Save_Occupancies_Harmonize.m` (via `combat/combat.m`) | Compute the fractional occupancy of each mode per scan, and, optionally (`apply_combat`), harmonize it across acquisition sites with ComBat, keeping diagnosis/age/sex/education as covariates of interest. Saves both `P_original` and `P_harmonized`. |
| 3a (optional) | `LEiDA_stats_Voxel_FracOccup_ComBat.m` | Test each mode's occupancy between conditions (Welch's t-test for independent samples, paired permutation test for paired samples), with permutation p-values and Hedge's effect sizes. Runs on either `P_original` or `P_harmonized`, chosen when calling it. Skip for studies with no discrete conditions to compare. |
| 3b | `Scores_vs_Mode_Occupancy.m` | Correlates every clinical/cognitive score against every mode in the entire pyramid (all K), partial correlation controlling for age. Reports per-score Bonferroni significance (`0.05/sum(rangeK)/N_scores`) to the command line and saves the p-values — no figure, no dependency on any pre-selected set of modes or on step 3a. The saved p-values can be rendered by `Plot_ClustVoxelCentroid_Pyramid_RSNs.m` in place of the permutation-test p-values. |
| 4 | `Plot_FracOccup_stats.m` | Summary plots of the statistical tests: p-values, mean occupancy barplots, effect sizes. Requires step 3a. |
| 5 | `Plot_ClustVoxelCentroid_Pyramid_RSNs.m` | Render a pyramid of all centroids across K, with optional Yeo RSN color overlay and significance markers for a chosen `stat_of_interest` (a condition pair from step 3a, or a score from step 3b). |
| 6 | `Choose_Relevant_Modes.m` | Automatically select the `[ki c]` modes that differ most between conditions (significant after multiple-testing correction, minimum effect size), grouping correlated modes. Requires step 3a; without it, specify `Key_Modes_KC` manually instead. |
| 7 | `Plot_KeyModes_Slices_Stats.m` | Render each key mode on anatomical slices with mean ± SE occupancy bars per condition. Requires step 3a. |
| 8 | `Plot_Mode_TransparentBrain.m` | Detailed 3D rendering of each key mode, including overlap with Yeo Resting-State Networks. Only needs `cluster_file` and `Key_Modes_KC` — does not require step 3a. |
| 9 | `Plot_KeyModes_vs_Scores.m` | Figure-generation counterpart to step 3b: plots a bar of score correlations per selected key mode and exports a CSV. Only needs `P` and `Key_Modes_KC` — does not require step 3a. |

Figures are saved at each step in the results folder in both `.fig` and `.png`/`.jpg`.

## Repository structure

```
run_LEiDA_Voxel.m                    Main pipeline script (documents and runs all steps)
Preprocessing/                       C-PAC pipeline used to preprocess the raw ADNI fMRI data (upstream of step 1)
CodeOcean_Capsule/                   Standalone Code Ocean capsule (self-contained copy of the pipeline)
  code/run_LEiDA_Voxel_CodeOcean.m     Capsule entry point (steps 2-4, from pre-extracted eigenvectors)
  data/                                Put demo eigenvector + Scores files here (not tracked in git)
  results/                             Outputs written here automatically
Select_Demo_Subsample.m              One-off local script: build a balanced demo subsample for the capsule's data/
Mask_Voxels_of_Interest.m            Step 0: build a custom voxel mask
Get_EigenVectors_VoxelSpace_Server.m Step 1: extract leading eigenvectors from fMRI data
LEiDA_cluster_VoxelMNI10mm.m         Step 2: K-means clustering
Save_Occupancies_Harmonize.m         Step 2b: mode occupancy extraction + optional ComBat harmonization
LEiDA_stats_Voxel_FracOccup_ComBat.m Step 3a: statistics on mode occupancy between conditions
Scores_vs_Mode_Occupancy.m           Step 3b: pyramid-wide mode-occupancy vs. scores statistics (no figure)
Plot_FracOccup_stats.m               Step 4: statistics summary plots
Plot_ClustVoxelCentroid_Pyramid_RSNs.m  Step 5: centroid pyramid render
Choose_Relevant_Modes.m              Step 6: automatic key-mode selection
Plot_KeyModes_Slices_Stats.m         Step 7: key-mode slice + occupancy plots
Plot_Mode_TransparentBrain.m         Step 8: detailed 3D key-mode render
Plot_KeyModes_vs_Scores.m            Step 9: key-modes vs. scores figure + CSV
combat/                              ComBat site-harmonization toolbox (Johnson et al.)
utilities/                           Colormaps, MNI masks, Yeo RSN parcellation, plotting/stats helpers
```

## Usage

1. Edit the "USER INPUT" section inside `Get_EigenVectors_VoxelSpace_Server.m` (data directory, filename pattern, voxel mask, output filename) and run it to extract leading eigenvectors from your fMRI data.
2. Edit the directories and filenames at the top of `run_LEiDA_Voxel.m` (`fMRI_dir`, `leida_dir`, `Scores_Table`, `results_dir`, `file_V1`, `cluster_file`, `stats_file`) and run the script section by section to cluster, harmonize, test, and visualize the coupling modes.

See the header comments in `run_LEiDA_Voxel.m` for the full function reference, input/output descriptions, and a worked example of each step.

## Reproducing the manuscript's results

Running the pipeline on the full 2177-scan ADNI cohort with the settings below reproduces the analysis reported in the manuscript (see Methods, sections 4.5-4.6):

| Step | Parameter | Value |
|---|---|---|
| 2. Clustering | `mink`, `maxk` | 2, 20 (209 total cluster centroids, Σ K for K=2..20) |
| 2. Clustering | `replicates` | 20 (lowest total cluster-to-centroid distance retained per K) |
| 2b. Occupancies | `apply_combat` | 1 (site harmonization, preserving diagnostic group, age, education, and sex as covariates) |
| 3a. Condition statistics | test | Welch's t-test (unequal variances), two-sided |
| 3a. Condition statistics | `n_permutations` | 100,000 (group labels randomly permuted to build an empirical null per pairwise comparison) |
| 3a. Condition statistics | `n_bootstraps` | 0 (sample sizes are large enough that no within-permutation bootstrap variance estimate was needed) |
| 3a. Condition statistics | effect size | Hedges' g (sample-size-corrected, preferred over Cohen's d for the unequal group sizes here) |

Expected runtime for the full dataset with these settings: ~24 hours on a normal desktop (see [Requirements](#requirements)).

## Code Ocean capsule

[`CodeOcean_Capsule/`](CodeOcean_Capsule/) is a self-contained copy of the pipeline, structured as a standard Code Ocean capsule (`code/`, `data/`, `results/`), starting from already-extracted leading eigenvectors (i.e. skipping step 0/1, which need the raw fMRI NIfTI files). It's independent of the rest of this repository — `code/` already has its own copies of every function file it needs (including `combat/` and `utilities/`), so the folder can be uploaded to Code Ocean as-is.

1. On your own machine, with the full cohort's already-extracted eigenvectors and `Scores_ADNI` table, run [`Select_Demo_Subsample.m`](Select_Demo_Subsample.m) to pick a small demo subsample (`n_per_condition`, 60 scans per condition by default), stratified by sex and age-tertile so the groups are balanced, and save the two demo files. Example output (from an earlier run with `n_per_condition = 30`; the printed report scales with whatever `n_per_condition` you use):
   ```
   CN:  30 participants selected, age 76.2 +/- 8.7, 15 male / 15 female
   MCI: 30 participants selected, age 74.8 +/- 8.0, 15 male / 15 female
   DEM: 30 participants selected, age 73.5 +/- 7.8, 15 male / 15 female
   ```
2. Copy those two files into [`CodeOcean_Capsule/data/`](CodeOcean_Capsule/data/) (named to match `file_V1`/`Scores_Table` inside `CodeOcean_Capsule/code/run_LEiDA_Voxel_CodeOcean.m`, or edit those two lines to match your filenames).
3. Run [`CodeOcean_Capsule/code/run_LEiDA_Voxel_CodeOcean.m`](CodeOcean_Capsule/code/run_LEiDA_Voxel_CodeOcean.m) (it defaults to reading `../data/` and writing `../results/`, matching Code Ocean's convention). It runs Step 2 (cluster K=2:20), Step 2b (occupancies), Step 3a (condition statistics), Step 3b (pyramid-wide score correlations), and Step 4 (all figures, including key-mode selection and the key-modes-vs-scores figure) — Step 1 (eigenvector extraction) is the offline step you did in (1) above, since it needs the raw fMRI data that isn't part of this capsule.

See [`CodeOcean_Capsule/README.md`](CodeOcean_Capsule/README.md) for capsule-specific notes.

ComBat harmonization is off by default in the capsule (`apply_combat = 0` inside `run_LEiDA_Voxel_CodeOcean.m`) since a small demo subsample is unlikely to have enough scans per site for reliable harmonization — turn it on if your demo data spans multiple well-populated sites. If no mode survives the significance threshold on the reduced sample (likely with far fewer scans than the full study), the driver falls back to a fixed mid-K mode selection so the figure-generation steps still produce output, and logs a warning explaining why.

**Keeping the capsule in sync**: `CodeOcean_Capsule/code/` holds copies, not symlinks, of the shared pipeline functions (everything except `Get_EigenVectors_VoxelSpace_Server.m`, `Mask_Voxels_of_Interest.m`, and `run_LEiDA_Voxel.m` itself, which the capsule doesn't need). If you edit one of the shared functions at the repository root, copy the change into `CodeOcean_Capsule/code/` too.

## Notes

- **`Key_Modes_KC`'s first column (`ki`) is a POSITION into `rangeK`, not the literal number of clusters.** `Kmeans_results{ki}` and `P(:,ki,:)` are indexed by position (e.g. if `rangeK = 2:20`, `ki=1` means K=2 clusters, `ki=2` means K=3). `Choose_Relevant_Modes.m` returns `ki` correctly; when building `Key_Modes_KC` manually, convert with `ki = find(rangeK == K_wanted)` rather than using the literal K value directly — this was a real bug in earlier versions of `Plot_KeyModes_Slices_Stats.m`/`Plot_Mode_TransparentBrain.m`/`Choose_Relevant_Modes.m` (which only happened to work when `mink=1`), now fixed to use `ki` consistently everywhere.
- **Studies with only continuous scores, no discrete conditions**: skip step 3a (`LEiDA_stats_Voxel_FracOccup_ComBat`) and steps 4, 6, 7 that depend on it. Run steps 0-2b, then 3b (`Scores_vs_Mode_Occupancy`, which never needed step 3a), pick `P = P_original` or `P = P_harmonized`, choose `Key_Modes_KC` manually (e.g. `Key_Modes_KC = [2 4; 4 6]`, using the `ki` convention above), then call `Plot_Mode_TransparentBrain(results_dir, cluster_file, Key_Modes_KC)` and `Plot_KeyModes_vs_Scores(P, Scores_Table, Key_Modes_KC, results_dir, save_name)` directly — neither depends on step 3a's output.
- `Save_Occupancies_Harmonize.m` saves both the raw (`P_original`) and ComBat-harmonized (`P_harmonized`) occupancies; `run_LEiDA_Voxel.m` exposes a `use_harmonized_occupancies` toggle to pick which one `P` refers to for the rest of the pipeline (steps 3a, 3b, and 9).
- `Scores_vs_Mode_Occupancy.m` and `Plot_KeyModes_vs_Scores.m` use a fixed set of score-table column indices tailored to the ADNI `Scores_ADNI` table used in Campo et al.; adapt these indices when using a different scores table.
- **Pyramid-wide score p-values as an alternative to condition-comparison p-values**: `Scores_vs_Mode_Occupancy.m`'s `pyramid_stats_file` output can be passed directly as `stats_file` to `Plot_ClustVoxelCentroid_Pyramid_RSNs.m` — its `P_pval` has one "row" per score instead of per condition pair, so `stat_of_interest` becomes a score index (see the command-line report for which scores had significant modes). `Plot_ClustVoxelCentroid_Pyramid_RSNs.m` derives the number of conditions for its occupancy bar plots from `Index_Conditions` rather than `P_pval`'s size, so this works regardless of how many scores vs. conditions there are. The saved figure filename is automatically tagged with the condition pair or score name (e.g. `Pyramid_CN-AD_...` or `Pyramid_MoCA_...`), so different comparisons don't overwrite each other.
- `Plot_KeyModes_Slices_Stats.m` additionally takes `Scores_Table` (used for sex-based grouping via `PTGENDER`).

## Author

Joana Cabral, University of Lisbon — joanabcabral@tecnico.ulisboa.pt
