# data/

Put your demo input files here:

- `LEiDA_V1_all_MNI10mm_demo.mat` — leading eigenvectors (output of `Get_EigenVectors_VoxelSpace_Server.m`, subsetted with `Select_Demo_Subsample.m`): `V1_all`, `ind_voxels`, `MNI_lowres_Mask`, `data_info`, `Scan_num`, `Scan_length`.
- `Scores_ADNI_demo.mat` — the `Scores_ADNI` table (`SITE`, `AGE_AT_SCAN`, `PTGENDER`, `PTEDUCAT`, `DX_num`, `DX`, plus the score columns used by `Scores_vs_Mode_Occupancy.m`), subsetted to the same scans.

If you name the files differently, edit `file_V1`/`Scores_Table` near the top of `../code/run_LEiDA_Voxel_CodeOcean.m` to match.

This folder is intentionally empty in the repository (only this README is tracked) — the actual `.mat` data files are not committed to git.
