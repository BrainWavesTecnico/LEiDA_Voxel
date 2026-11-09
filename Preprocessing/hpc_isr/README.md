# Running C-PAC on the ISR HPC cluster

Scripts to run the C-PAC preprocessing in [`../`](../README.md) on the ISR HPC
cluster ([docs](https://github.com/hpc-isr/hpc-docs)). The cluster has no
Docker: C-PAC's Docker image runs with **Apptainer**, inside **Slurm** jobs,
one job-array task per subject.

- **`submit_cpac.sh`** — run on `hpc01`. Lists the subjects of a BIDS dataset
  that aren't done yet and submits them as a job array.
- **`cpac_subject.sbatch`** — one task = one subject. Builds the C-PAC image
  (`fcpindi/c-pac:release-v1.8.7`) on the node's `/scratch` if needed, copies
  the subject there, runs C-PAC with the job's 8 cores / 32 GB, copies the
  results to `RESULTS_DIR`, and cleans up.

## Step by step

**1. Get an account and log in** (see the cluster's
[Getting an account](https://github.com/hpc-isr/hpc-docs/blob/main/docs/user-guide/01-getting-an-account.md)):

    ssh <username>@hpc01.isr.tecnico.ulisboa.pt

**2. Get this repository on the cluster** (in your home, visible from every node):

    git clone https://github.com/BrainWavesTecnico/LEiDA_Voxel.git ~/LEiDA_Voxel

**3. Copy your dataset from your computer** (run on your computer; `rsync`
resumes if interrupted):

    rsync -avP /path/on/your/computer/my_dataset/ <username>@hpc01.isr.tecnico.ulisboa.pt:data/my_dataset/raw/

Home is limited to 200 GiB (500 GiB hard). If your group has NAS storage
(`/storage/anubis-nas/<group>`), use it for large datasets and results.

**4. Make sure the data is in BIDS format.** C-PAC needs:

    BIDS/
      dataset_description.json
      sub-01/anat/sub-01_T1w.nii.gz
      sub-01/func/sub-01_task-rest_bold.nii.gz
      sub-01/func/sub-01_task-rest_bold.json   (RepetitionTime, SliceTiming)

DICOMs must be converted first (e.g. with `dcm2bids` or `heudiconv`), in a
Slurm job, not on `hpc01`.

**5. Edit the paths at the top of `submit_cpac.sh`:** `BIDS_DIR`,
`RESULTS_DIR` and `CONFIG_FILE` (by default the linear/MCFLIRT config).

**6. Test one subject:**

    cd ~/LEiDA_Voxel/Preprocessing/hpc_isr
    ./submit_cpac.sh 1
    squeue -u "$USER"

The first task on a node first builds the C-PAC image (~tens of minutes).
Follow it in `RESULTS_DIR/slurm/cpac-<jobid>_1.out`. A subject takes a few
hours. Check the outputs in `RESULTS_DIR/output/` before going on.

**7. Run all the others:**

    ./submit_cpac.sh

At most `MAX_PARALLEL` (default 4) subjects run at once. Re-run the script
any time: it only submits subjects without a done marker
(`RESULTS_DIR/.done/sub-*`), e.g. after failures.

**8. Look at failures:** `RESULTS_DIR/slurm/*.out`, `RESULTS_DIR/log/` and
`RESULTS_DIR/crash/`. `sacct -X -j <jobid>` shows `OUT_OF_MEMORY` or
`TIMEOUT`; raise `--mem` / `--time` in `cpac_subject.sbatch` if needed.

## Space on `/scratch`

Each running subject uses up to tens of GB of C-PAC working files on its
node's `/scratch`, plus the ~5 GB image once per node. Your `/scratch` quota is
100 GiB (soft) per node, so keep `MAX_PARALLEL` at about 4 unless you have
checked actual use (`du -sh /scratch/$USER` in a job on that node).
