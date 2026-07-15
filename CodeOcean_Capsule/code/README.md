# code/

Entry point: **`run_LEiDA_Voxel_CodeOcean.m`** — run this script directly (it
is a plain script, not a function; it defines `data_dir`/`results_dir` at
the top, defaulting to `../data/` and `../results/`). See the top of that
file for the exact list of expected input files in `../data/`.

## What it does

`run_LEiDA_Voxel_CodeOcean.m` runs the pipeline from already-extracted
leading eigenvectors through to the final figures:

| Step | Function | Output (in `results/`) |
|---|---|---|
| 2. Cluster (K=2:20) | `LEiDA_cluster_VoxelMNI10mm.m` | `LEiDA_Clusters_VoxelMNI10mm_demo.mat` |
| 2b. Occupancies (+ optional ComBat) | `Save_Occupancies_Harmonize.m` | `LEiDA_Occupancies_demo.mat` |
| 3a. Condition statistics | `LEiDA_stats_Voxel_FracOccup_ComBat.m` | `LEiDA_Stats_FracOccup_demo.mat` |
| 3b. Score correlations (full pyramid) | `Scores_vs_Mode_Occupancy.m` | `Scores_Pyramid_Pval.mat` |
| 4. Condition-stats figures | `Plot_FracOccup_stats.m` | `Figure_FracOccup_pvalues.*`, `Figure_FracOccup_Barplot_*.*`, `Figure_FracOccup_effetcsize.*` |
| 4. Pyramid centroid figure | `Plot_ClustVoxelCentroid_Pyramid_RSNs.m` | `Centroid_Pyramid_RSNs_*.*` |
| 4. Select key modes | `Choose_Relevant_Modes.m` | (no file - returns `Key_Modes_KC`; falls back to showing the full K=4 repertoire if no mode survives significance) |
| 4. Key-mode slices figure | `Plot_KeyModes_Slices_Stats.m` | `Fig1_Key_modes_Decrease.*`, `Fig1_Key_modes_Increase.*` |
| 4. Key-mode 3D renders | `Plot_Mode_TransparentBrain.m` | `Fig3_KeyModes_TransparentBrain.*` |
| 4. Key-mode vs. scores | `Plot_KeyModes_vs_Scores.m` | `Fig_KeyModes_vs_Scores.*`, `BraVe_correlations.csv`, `Scores_Mode_Stats.mat` |

Step 1 (`Get_EigenVectors_VoxelSpace_Server.m`, which reads raw fMRI NIfTI
files) is **not** part of this capsule and is not run here — it needs the
full-size raw data. Run it offline beforehand (see the main repo
[`README.md`](../../README.md)), or use [`../../Select_Demo_Subsample.m`](../../Select_Demo_Subsample.m)
to build a demo subsample from an existing full-cohort eigenvector file.

## Dependencies (bundled)

- `combat/` — third-party ComBat site-harmonization toolbox (Johnson et al.),
  used by `Save_Occupancies_Harmonize.m` when `apply_combat=1`.
- `utilities/` — `subplot_tight.m` (tight subplot grids), the two
  permutation-test helpers used by `LEiDA_stats_Voxel_FracOccup_ComBat.m`,
  and the `.mat` assets the plotting functions load (colormap, MNI mask,
  atlas parcellations).

This folder is self-contained: every function file it needs is copied in
here (not referenced from the parent repo), so the whole `CodeOcean_Capsule/`
folder can be uploaded to Code Ocean as-is.

## Capsule-specific differences from the main repo copies

A few files here are intentionally simplified for the capsule, compared to
their counterparts one level up in the main repository:

- **`Plot_Mode_TransparentBrain.m`**: the AAL120/Desikan brain-area-list
  figures (`Fig2`/`Fig3`) are commented out, to keep the demo run's figure
  output to just the 3D renders + RSN overlap (`Fig1`). Uncomment the marked
  blocks to re-enable them.
- **`run_LEiDA_Voxel_CodeOcean.m`**: uses fewer K-means replicates and no
  permutation bootstraps by default (`replicates=20`, `n_bootstraps=0`) to
  keep runtime reasonable on a demo-sized sample; `apply_combat=0` by
  default since a small demo subsample is unlikely to have enough scans per
  site for reliable harmonization.

If you copy a fix from the main repo's function files into this folder,
keep these differences in mind rather than overwriting them wholesale.
