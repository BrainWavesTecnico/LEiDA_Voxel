#!/bin/bash
# One C-PAC run per array task (one subject each), on the ISR HPC cluster.
# Submit it with submit_cpac.sh, which writes the subject list and sets the
# variables below; don't sbatch this file directly.
#
# Per task it:
#   1. builds the C-PAC Apptainer image on this node's /scratch if it isn't there
#      (once per node; Docker isn't available on the cluster);
#   2. copies that subject's BIDS folder to /scratch;
#   3. runs C-PAC with the cores and memory Slurm gave the job;
#   4. copies outputs, logs and crash files to RESULTS_DIR and cleans /scratch,
#      however the job ends.
#SBATCH --job-name=cpac
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --open-mode=append

set -euo pipefail

: "${BIDS_DIR:?set by submit_cpac.sh}"
: "${RESULTS_DIR:?set by submit_cpac.sh}"
: "${CONFIG_FILE:?set by submit_cpac.sh}"
: "${SUBJECT_LIST:?set by submit_cpac.sh}"
CPAC_IMAGE="${CPAC_IMAGE:-docker://fcpindi/c-pac:release-v1.8.7}"
SIF="/scratch/$USER/containers/$(basename "${CPAC_IMAGE//:/-}").sif"

SUB=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SUBJECT_LIST")
[ -n "$SUB" ] || { echo "No subject on line $SLURM_ARRAY_TASK_ID of $SUBJECT_LIST"; exit 1; }
DONE="$RESULTS_DIR/.done/sub-$SUB"
echo "=== sub-$SUB on $(hostname), task $SLURM_ARRAY_TASK_ID, job $SLURM_JOB_ID, $(date)"

# A requeued task that already finished must not redo the work.
if [ -f "$DONE" ]; then
    echo "sub-$SUB is already done ($DONE), skipping."
    exit 0
fi

WORK="/scratch/$USER/cpac-$SLURM_JOB_ID"
mkdir -p "$WORK"/{bids,outputs,tmp,home} "$(dirname "$SIF")" "$RESULTS_DIR/.done"

finish() {
    local status=$?
    # Outputs of all subjects share one tree: output/<pipeline_name>/sub-*/...
    for d in output log crash; do
        if [ -d "$WORK/outputs/$d" ]; then
            mkdir -p "$RESULTS_DIR/$d"
            cp -r "$WORK/outputs/$d/." "$RESULTS_DIR/$d/"
        fi
    done
    rm -rf "$WORK"
    echo "=== sub-$SUB ended with status $status at $(date)"
}
trap finish EXIT
trap 'exit 143' TERM   # time limit or preemption: still run finish

# 1. The image, once per node. The lock makes tasks that start together on one
#    node wait for a single build.
(
    flock 9
    if [ ! -f "$SIF" ]; then
        echo "Building $SIF from $CPAC_IMAGE (tens of minutes, once per node)"
        export APPTAINER_CACHEDIR="$WORK/apptainer/cache" APPTAINER_TMPDIR="$WORK/apptainer/tmp"
        mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
        apptainer build "$SIF.part" "$CPAC_IMAGE" && mv "$SIF.part" "$SIF"
        rm -rf "$WORK/apptainer" "$SIF.part"
        [ -f "$SIF" ]
    fi
) 9>"$SIF.lock"

# 2. This subject only, plus the dataset-level files C-PAC and the BIDS
#    validator read (dataset_description.json, task-*_bold.json, ...).
find "$BIDS_DIR" -maxdepth 1 -type f \( -name '*.json' -o -name '*.tsv' -o -name '.bidsignore' \) \
    -exec cp {} "$WORK/bids/" \;
cp -r "$BIDS_DIR/sub-$SUB" "$WORK/bids/"

# 3. Run C-PAC with what Slurm allocated: these flags override the cores,
#    memory and ANTs threads written in the pipeline config. Memory is kept a
#    little under the Slurm limit so the job isn't killed for going over it.
NCPU="$SLURM_CPUS_PER_TASK"
MEM_GB=$(( SLURM_MEM_PER_NODE / 1024 * 85 / 100 ))
echo "C-PAC: $NCPU cores, ${MEM_GB} GB, config $CONFIG_FILE"

srun apptainer run --cleanenv \
    --home "$WORK/home" \
    --bind "$WORK/bids:/bids:ro,$WORK/outputs:/outputs,$WORK/tmp:/tmp,$CONFIG_FILE:/pipeline_config.yml:ro" \
    "$SIF" /bids /outputs participant \
    --pipeline_file /pipeline_config.yml \
    --participant_label "$SUB" \
    --n_cpus "$NCPU" --num_ants_threads "$NCPU" --mem_gb "$MEM_GB"

touch "$DONE"
