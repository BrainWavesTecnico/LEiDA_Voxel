# data/

demo input files:

- `LEiDA_V1_all_MNI10mm_s300demo.mat` — leading eigenvectors (output of `Get_EigenVectors_VoxelSpace_Server.m`, subsetted with `Select_Demo_Subsample.m`): `V1_all`, `ind_voxels`, `MNI_lowres_Mask`, `data_info`, `Scan_num`, `Scan_length`.
- `Scores_ADNI_s300demo.mat` — the `Scores_ADNI` table (`SITE`, `AGE_AT_SCAN`, `PTGENDER`, `PTEDUCAT`, `DX_num`, `DX`, plus the score columns used by `Scores_vs_Mode_Occupancy.m`), subsetted to the same scans.

