# Preprocessing/

The fMRI preprocessing pipeline used to prepare the ADNI resting-state data
for this repository's LEiDA analysis, before eigenvector extraction (step 1,
[`Get_EigenVectors_VoxelSpace_Server.m`](../Get_EigenVectors_VoxelSpace_Server.m)).
Preprocessing was run with **C-PAC** (Configurable Pipeline for the Analysis
of Connectomes, v1.8.7.dev1), via Docker.

## Files

- **`run_cpapc_adni_ants.sh`** — batch driver script. Runs the C-PAC Docker
  image (`fcpindi/c-pac:latest`) once per subject over a BIDS-formatted
  dataset, in configurable batches (`BATCH_SIZE`, default 100 subjects) with
  limited parallelism (`PARALLEL_JOBS`, default 3 concurrent subjects).
  Auto-detects and skips subjects already marked as completed in the log
  file, so it can be re-run repeatedly to work through a large cohort in
  batches. Edit `BIDS_DIR`, `OUTPUT_DIR`, and `CONFIG_FILE` at the top to
  match your paths before running.
- **`cpac_pipeline_with_mean_best_for_csf_0.01_0.1.yml`** — the C-PAC
  pipeline configuration used for this dataset (referenced by
  `run_cpapc_adni_ants.sh` via `CONFIG_FILE`).

## Pipeline configuration summary

**Anatomical preprocessing**
- N4 bias field correction.
- Skull-stripping with ANTs (`niworkflows-ants`, OASIS template).
- Tissue segmentation into CSF/WM/GM with FSL-FAST, using tissue priors.
- Registration to MNI152 (2mm) with ANTs (Rigid → Affine → SyN), Lanczos
  windowed-sinc interpolation.

**Functional preprocessing**
- Slice-timing correction.
- Motion correction with AFNI `3dvolreg` (two-pass), mean-volume reference.
- No fieldmap-based distortion correction (not available for this dataset).
- Functional brain masking via an anatomically-refined (dilated) mask.
- Functional-to-anatomical coregistration with FSL boundary-based
  registration (BBR, 6 DOF).
- Functional-to-MNI normalization (2mm output resolution) via the T1
  template pathway.

**Nuisance regression / denoising** (no ICA-AROMA)
- Regressors: 24-parameter motion model (raw + squared + temporal
  derivatives + squared derivatives) and mean white-matter signal (no CSF
  regressor) with 2nd-order polynomial detrending.
- Band-pass filter 0.01-0.10 Hz, applied after regression.
- CSF/ventricle signal is deliberately **not** regressed out (only white
  matter is used as the nuisance tissue signal) — the ventricle signal is
  needed intact for the brain-ventricle coupling ("BraVe") analysis this
  pipeline is built around.

**Post-processing**
- Spatial smoothing at 5mm FWHM (both smoothed and non-smoothed outputs are
  kept).

Outputs are written in MNI152 space (2mm), which is what
[`Get_EigenVectors_VoxelSpace_Server.m`](../Get_EigenVectors_VoxelSpace_Server.m)
expects as input.
