# Preprocessing/

The fMRI preprocessing pipeline used to prepare the ADNI resting-state data
for this repository's LEiDA analysis, before eigenvector extraction (step 1,
[`Get_EigenVectors_VoxelSpace_Server.m`](../Get_EigenVectors_VoxelSpace_Server.m)).
Preprocessing was run with **C-PAC** (Configurable Pipeline for the Analysis
of Connectomes, v1.8.7.dev1), via Docker.

## Files

- **`pod_to_bids.sh`** — converts the PODelirium raw data (`POD/anat`,
  `POD/func`) to BIDS, with the slice timing; re-run it to add new files.
- **`run_cpapc_adni_ants.sh`** — batch driver script. Runs the C-PAC Docker
  image (`fcpindi/c-pac:latest`) once per subject over a BIDS-formatted
  dataset, in configurable batches (`BATCH_SIZE`, default 100 subjects) with
  limited parallelism (`PARALLEL_JOBS`, default 3 concurrent subjects).
  Auto-detects and skips subjects already marked as completed in the log
  file, so it can be re-run repeatedly to work through a large cohort in
  batches. Edit `BIDS_DIR`, `OUTPUT_DIR`, and `CONFIG_FILE` at the top to
  match your paths before running. Optional settings below them choose the
  C-PAC image tag (`CPAC_IMAGE`), the Docker platform (`DOCKER_PLATFORM`,
  `linux/amd64` on Apple Silicon Macs) and the cores and memory per subject
  (`N_CPUS`, `MEM_GB`), overriding the values in the pipeline config. Left
  empty, the script behaves as before.
- **`cpac_pipeline_with_mean_best_for_csf_0.01_0.1.yml`** — the C-PAC
  pipeline configuration used for this dataset (referenced by
  `run_cpapc_adni_ants.sh` via `CONFIG_FILE`).
- **`cpac_pipeline_linear_mcflirt_0.01_0.1.yml`** — variant of the config
  above that avoids interpolation ringing (spurious spatial waves in the
  MNI-space BOLD, strongest at sharp edges such as the ventricles). It is
  identical except for:
  - functional-to-MNI resampling with `Linear` instead of
    `LanczosWindowedSinc` (`func_registration_to_template.ANTs_pipelines.interpolation`);
  - motion correction with FSL MCFLIRT instead of AFNI `3dvolreg`
    (`motion_correction.using`);
  - `pipeline_name` set to `cpac-adni-ants-linear-mcflirt`, so its outputs
    are kept separate from the original pipeline's.

- **`cpac_pipeline_linear_mcflirt_nobandpass.yml`** — the same as
  `cpac_pipeline_linear_mcflirt_0.01_0.1.yml` but **without the band-pass
  filter**: the nuisance regression (24 motion parameters, mean white
  matter, 2nd-order polynomial detrending) is unchanged, and the data is left
  unfiltered so filtering can be done later. `pipeline_name` is
  `cpac-adni-ants-linear-mcflirt-nobandpass`.

- **`cpac_pipeline_linear_mcflirt_3mm_nobandpass.yml`** — the same as
  `cpac_pipeline_linear_mcflirt_nobandpass.yml` (no band-pass filter) but
  with functional outputs in MNI space at **3 mm** instead of 2 mm, closer to
  the ~3.4 mm native EPI voxels, so less upsampling. `pipeline_name` is
  `cpac-adni-ants-linear-mcflirt-3mm-nobandpass`.

  To use any of these variants, point `CONFIG_FILE` in
  `run_cpapc_adni_ants.sh` at it.
  Use a new `OUTPUT_DIR` too: the script skips subjects already marked as
  completed in that directory's log.

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

## PODelirium: step by step

How the PODelirium data are preprocessed, on a Mac or Linux computer with
Docker. On Windows, use WSL2 with Docker Desktop and run the same commands
in the Ubuntu terminal.

1. **Install Docker Desktop** (https://www.docker.com/products/docker-desktop/).
   On an Apple Silicon Mac, turn on *Settings → General → Use Rosetta for
   x86_64/amd64 emulation*, and give Docker most of the memory in
   *Settings → Resources* (at least 26 GB per subject processed at once).
   Check it works: `docker run --rm hello-world`.
2. **Get the scripts** (this branch):

       git clone -b delirium-data-preprocessing https://github.com/BrainWavesTecnico/LEiDA_Voxel.git ~/LEiDA_Voxel
       cd ~/LEiDA_Voxel/Preprocessing

3. **Download C-PAC** (about 4 GB, once):

       docker pull --platform linux/amd64 fcpindi/c-pac:release-v1.8.7

4. **Convert the raw data to BIDS** with `pod_to_bids.sh`. The raw folder
   (`POD/anat`, `POD/func`) is only read. Do a dry run first, then add `--go`:

       bash pod_to_bids.sh /path/to/PODelirium/POD /path/to/PODelirium/POD_BIDS
       bash pod_to_bids.sh /path/to/PODelirium/POD /path/to/PODelirium/POD_BIDS --go

   It also writes the acquisition timing (TR = 2 s, interleaved slice order)
   that slice-timing correction needs. If new raw files arrive later, run it
   again: it only adds the new ones.
5. **Run C-PAC on all subjects** with `run_cpac_all.sh` (described below):

       ./run_cpac_all.sh /path/to/PODelirium/POD_BIDS /path/to/PODelirium/POD_cpac_A_3mm

   It prints a summary first (subjects, scans, cores, memory, disk space).
   Each scan takes roughly 15–30 minutes. Follow progress from another
   terminal window with
   `tail -f /path/to/PODelirium/POD_cpac_A_3mm/processing_log.txt`.
   If some subjects fail, run the same command again: finished subjects are
   skipped, and subjects with new sessions are processed again.
6. **Results** for LEiDA are in
   `POD_cpac_A_3mm/output/pipeline_cpac-adni-ants-linear-mcflirt-3mm/sub-*/ses-*/func/`:
   `*_space-MNI152NLin6ASym_desc-preproc_bold.nii.gz` (3 mm MNI, band-passed
   0.01–0.1 Hz) and its 5 mm smoothed version.

Steps of the pipeline (`cpac_pipeline_linear_mcflirt_3mm_0.01_0.1.yml`):
slice-timing correction → motion correction (FSL MCFLIRT) → coregistration
to the T1 (FSL BBR) → nuisance regression (24 motion parameters + mean white
matter, quadratic detrending; no CSF regression) → band-pass 0.01–0.1 Hz →
normalization to MNI152 with ANTs, Linear interpolation to 3 mm → 5 mm
smoothing.

Note: the data contain an alternating odd/even slice artefact (interleaved
acquisition with breathing-related head motion). It appears as stripes after
normalization and is removed by the 5 mm smoothing; unsmoothed outputs still
contain it. See `destripe_slices.py` for an optional native-space correction.

## Running a whole dataset

**`run_cpac_all.sh`** runs C-PAC with Docker on every subject (all sessions)
of a BIDS dataset, using `run_cpapc_adni_ants.sh`:

    ./run_cpac_all.sh <BIDS folder> <output folder> [config file]

The default config is `cpac_pipeline_linear_mcflirt_3mm_0.01_0.1.yml`
(slice-timing correction, MCFLIRT, nuisance regression then 0.01–0.1 Hz
band-pass, Linear resampling to 3 mm MNI, 5 mm smoothing). Before starting it
checks the dataset (and that `SliceTiming` is present when slice-timing
correction is on), the config, Docker and the free disk space. It chooses how
many subjects run at once from the cores and the memory Docker has (override
with `PARALLEL_JOBS`, `N_CPUS`, `MEM_GB`), uses `linux/amd64` on Apple
Silicon, and keeps a Mac awake while it runs. Finished subjects are skipped,
so after failures just run it again.

## Band-pass before slice-timing correction

C-PAC always runs slice-timing correction before its band-pass. To filter
first, **`bandpass_bids.sh`** band-passes every BOLD run of a BIDS dataset
in native space (AFNI `3dTproject -passband`, default 0.01–0.1 Hz, quadratic
detrending), adds each voxel's mean back so C-PAC can still mask and register
the images, copies all other files unchanged, and writes a new BIDS dataset.
Run it in the C-PAC container (see the script header), then run C-PAC on the
new dataset. **`cpac_pipeline_linear_mcflirt_3mm_0.01_0.1.yml`** is the
matching config: Linear/MCFLIRT, slice-timing correction, 3 mm functional
outputs, 5 mm smoothing, and C-PAC's own 0.01–0.1 Hz band-pass after
nuisance regression, which removes what regressing the unfiltered motion
parameters puts back outside the band.

**`cpac_pipeline_linear_mcflirt_3mm_0.01_0.1_bpbefore.yml`** (not
recommended, kept as a record of a test) sets
`bandpass_filtering_order: 'Before'`. In C-PAC 1.8.7 this band-passes the
BOLD data before nuisance regression, but the regression uses the
*unfiltered* regressors (the filtered copy is only written as an output) and
nothing is filtered afterwards, so frequencies outside 0.01–0.1 Hz are
reintroduced (Lindquist et al. 2019).

## Affine, nearest-neighbour normalization (no interpolation smoothing)

**`native_to_mni_affine_nn.sh`** maps C-PAC's native-space preprocessed BOLD
(`*_desc-preproc_bold.nii.gz` without `space-`) to MNI152 with a single
resampling step that does not mix voxel values: an **affine-only** transform
(func→T1 from C-PAC, or an ANTs rigid fit if C-PAC didn't write it; T1→MNI
from an ANTs Rigid+Affine fit, no SyN) applied with **nearest-neighbour**
interpolation, so each MNI voxel takes the value of exactly one native voxel.
C-PAC itself only offers Linear, BSpline or Lanczos for this step. Run it in
the C-PAC container after C-PAC, per subject/session; see the header of the
script for the command. Note that C-PAC's motion correction has already
resampled the native data once (MCFLIRT, trilinear) before this step.

With **`--destripe`**, it first runs **`destripe_slices.py`** on the native
data. With interleaved slice acquisition, head motion within a TR (e.g.
breathing-related nodding) gives odd and even slices different intensity
from volume to volume: uniform within each slice, alternating along the slice
axis, with a strength that can change sign along it. After normalization this
appears as stripes across MNI slices, and it survives band-pass filtering.
For each volume, the script fits each slice with a plane (an offset plus
left-right and front-back gradients: the shape a small head rotation such
as nodding produces), then fits each of these components across the 5
nearest slices with a local linear trend plus an alternating term, and
subtracts the alternating plane from that slice (`--model mean` removes
only the offset). Patterns within a slice beyond a plane, signals that vary
smoothly across slices and voxel means are kept. It saves the cleaned
native data, the subtracted components (`*_desc-destripe_artefact-*.tsv`,
slices × volumes) and prints how much variance was removed.

## Running on the ISR HPC cluster

The cluster has no Docker. [`hpc_isr/`](hpc_isr/README.md) has Slurm and
Apptainer scripts that run the same C-PAC pipeline there, one job-array task
per subject, with step-by-step instructions.
