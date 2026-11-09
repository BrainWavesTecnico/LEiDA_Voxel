#!/bin/bash
# Submit C-PAC for every subject of a BIDS dataset that isn't done yet, as a
# Slurm job array (one task per subject). Run it on hpc01:
#
#   ./submit_cpac.sh          # all remaining subjects
#   ./submit_cpac.sh 1        # only the first remaining subject (a test run)
#
# Re-running it later submits only the subjects without a done marker, e.g.
# after failures. Edit the paths below first, or set them in the environment.
set -euo pipefail

BIDS_DIR="${BIDS_DIR:-$HOME/data/my_dataset/BIDS}"
RESULTS_DIR="${RESULTS_DIR:-$HOME/data/my_dataset/cpac_linear_mcflirt}"
CONFIG_FILE="${CONFIG_FILE:-$HOME/LEiDA_Voxel/Preprocessing/cpac_pipeline_linear_mcflirt_0.01_0.1.yml}"
MAX_PARALLEL="${MAX_PARALLEL:-4}"   # subjects running at once (each 8 cores, 32 GB)
LIMIT="${1:-0}"                     # 0 = no limit

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for p in "$BIDS_DIR" "$CONFIG_FILE"; do
    [ -e "$p" ] || { echo "Not found: $p (edit the paths at the top of $0)"; exit 1; }
done
[ -f "$BIDS_DIR/dataset_description.json" ] || { echo "$BIDS_DIR has no dataset_description.json: is it a BIDS dataset?"; exit 1; }

mkdir -p "$RESULTS_DIR/slurm" "$RESULTS_DIR/.done"
LIST="$RESULTS_DIR/slurm/subjects-$(date +%Y%m%d-%H%M%S).txt"

total=0
: > "$LIST"
for d in "$BIDS_DIR"/sub-*/; do
    s=$(basename "$d"); s=${s#sub-}
    total=$((total + 1))
    [ -f "$RESULTS_DIR/.done/sub-$s" ] || echo "$s" >> "$LIST"
done
if [ "$LIMIT" -gt 0 ]; then
    head -n "$LIMIT" "$LIST" > "$LIST.tmp" && mv "$LIST.tmp" "$LIST"
fi
n=$(wc -l < "$LIST")

echo "Subjects in $BIDS_DIR: $total, done: $(ls "$RESULTS_DIR/.done" | wc -l), submitting: $n"
if [ "$n" -eq 0 ]; then
    rm -f "$LIST"
    echo "Nothing to submit."
    exit 0
fi

sbatch --array="1-$n%$MAX_PARALLEL" \
    --output="$RESULTS_DIR/slurm/%x-%A_%a.out" \
    --export="ALL,BIDS_DIR=$BIDS_DIR,RESULTS_DIR=$RESULTS_DIR,CONFIG_FILE=$CONFIG_FILE,SUBJECT_LIST=$LIST" \
    "$HERE/cpac_subject.sbatch"
echo "Subject list: $LIST"
echo "Logs: $RESULTS_DIR/slurm/   Check with: squeue -u \"\$USER\""
