#!/bin/bash
# Run C-PAC on every subject (all sessions) of a BIDS dataset with Docker,
# using run_cpapc_adni_ants.sh, with resources chosen for this computer.
#
#   ./run_cpac_all.sh <BIDS folder> <output folder> [config file]
#
# Default config: cpac_pipeline_linear_mcflirt_3mm_0.01_0.1.yml (slice timing,
# MCFLIRT, nuisance regression then 0.01-0.1 Hz band-pass, Linear resampling
# to 3 mm MNI, 5 mm smoothing). Before starting it checks the dataset, the
# config, Docker and the free disk space, and keeps the computer awake.
# Subjects already done are skipped, so after failures just run it again.
#
# Optional environment variables: PARALLEL_JOBS, N_CPUS, MEM_GB (otherwise
# set from the cores and the memory Docker has; default 8 cores, at least the
# config's num_OMP_threads, and 24 GB per subject), CPAC_IMAGE.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIDS_DIR=${1:?usage: run_cpac_all.sh <BIDS folder> <output folder> [config file]}
OUTPUT_DIR=${2:?usage: run_cpac_all.sh <BIDS folder> <output folder> [config file]}
CONFIG_FILE=${3:-$HERE/cpac_pipeline_linear_mcflirt_3mm_0.01_0.1.yml}
CPAC_IMAGE=${CPAC_IMAGE:-fcpindi/c-pac:release-v1.8.7}
BIDS_DIR=$(cd "$BIDS_DIR" && pwd); CONFIG_FILE=$(cd "$(dirname "$CONFIG_FILE")" && pwd)/$(basename "$CONFIG_FILE")
mkdir -p "$OUTPUT_DIR"; OUTPUT_DIR=$(cd "$OUTPUT_DIR" && pwd)

fail() { echo "ERROR: $*"; exit 1; }

# --- Checks -------------------------------------------------------------
[ -f "$BIDS_DIR/dataset_description.json" ] || fail "$BIDS_DIR has no dataset_description.json"
[ -f "$CONFIG_FILE" ] || fail "config not found or not a file: $CONFIG_FILE"
grep -q "pipeline_name:" "$CONFIG_FILE" || fail "$CONFIG_FILE does not look like a C-PAC config"
[ -x "$HERE/run_cpapc_adni_ants.sh" ] || chmod +x "$HERE/run_cpapc_adni_ants.sh"

n_sub=$(find "$BIDS_DIR" -maxdepth 1 -type d -name 'sub-*' | wc -l | tr -d ' ')
n_bold=$(find "$BIDS_DIR" -name '*_bold.nii*' | wc -l | tr -d ' ')
[ "$n_sub" -gt 0 ] || fail "no sub-* folders in $BIDS_DIR"

if grep -q '"SliceTiming"' "$BIDS_DIR"/*_bold.json 2>/dev/null || \
   find "$BIDS_DIR" -name '*_bold.json' -exec grep -l '"SliceTiming"' {} + 2>/dev/null | grep -q .; then
    stc_info="found"
else
    stc_info="NOT found"
fi
if grep -A1 'slice_timing_correction:' "$CONFIG_FILE" | grep -q 'run: \[On\]' && [ "$stc_info" != "found" ]; then
    fail "slice-timing correction is on in the config, but no SliceTiming in the BIDS JSON files"
fi

docker info >/dev/null 2>&1 || fail "Docker is not running: open Docker Desktop and try again"
docker image inspect "$CPAC_IMAGE" >/dev/null 2>&1 || { echo "Pulling $CPAC_IMAGE ..."; docker pull --platform linux/amd64 "$CPAC_IMAGE"; }

# --- Resources ------------------------------------------------------------
if [ "$(uname)" = Darwin ]; then NCORES=$(sysctl -n hw.ncpu); else NCORES=$(nproc); fi
DOCKER_MEM_GB=$(( $(docker info --format '{{.MemTotal}}') / 1024 / 1024 / 1024 ))
N_CPUS=${N_CPUS:-8}
# C-PAC refuses to start if the config's num_OMP_threads exceeds the cores
# per subject, so never give a subject fewer cores than that.
OMP=$(grep -m1 'num_OMP_threads:' "$CONFIG_FILE" | awk '{print $2}')
case "$OMP" in ''|*[!0-9]*) OMP=1 ;; esac
[ "$N_CPUS" -ge "$OMP" ] || N_CPUS=$OMP
[ "$N_CPUS" -le "$NCORES" ] || N_CPUS=$NCORES
[ "$N_CPUS" -ge "$OMP" ] || fail "the config asks for num_OMP_threads: $OMP, but this computer has only $NCORES cores"
MEM_GB=${MEM_GB:-24}   # C-PAC refuses to start a run if a step is estimated to need more
if [ -z "${PARALLEL_JOBS:-}" ]; then
    by_mem=$(( (DOCKER_MEM_GB - 4) / MEM_GB ))     # leave ~4 GB for Docker itself
    by_cpu=$(( NCORES / N_CPUS ))
    PARALLEL_JOBS=$(( by_mem < by_cpu ? by_mem : by_cpu ))
    [ "$PARALLEL_JOBS" -ge 1 ] || PARALLEL_JOBS=1
fi

free_gb=$(df -Pk "$OUTPUT_DIR" | awk 'NR==2 {print int($4/1024/1024)}')
need_gb=$(( n_sub * 3 ))
platform=""; [ "$(uname -m)" = arm64 ] && platform="linux/amd64"

echo "================================================"
echo "BIDS dataset : $BIDS_DIR"
echo "               $n_sub subjects, $n_bold BOLD runs, SliceTiming $stc_info"
echo "Output       : $OUTPUT_DIR"
echo "Config       : $CONFIG_FILE"
echo "              ($(grep -m1 'pipeline_name:' "$CONFIG_FILE" | sed 's/^ *//'))"
echo "Image        : $CPAC_IMAGE ${platform:+(platform $platform)}"
echo "Computer     : $NCORES cores, Docker has ${DOCKER_MEM_GB} GB"
echo "Per subject  : $N_CPUS cores, $MEM_GB GB;  $PARALLEL_JOBS subjects at a time"
echo "Disk         : ${free_gb} GB free (roughly ${need_gb} GB needed)"
echo "================================================"
[ "$free_gb" -ge "$need_gb" ] || echo "WARNING: free disk space may not be enough."
[ "$PARALLEL_JOBS" -le 1 ] && [ "$DOCKER_MEM_GB" -lt 20 ] && \
    echo "WARNING: Docker has little memory. Raise it in Docker Desktop > Settings > Resources."

# --- Run --------------------------------------------------------------------
RUN="$OUTPUT_DIR/run_cpac_batch.sh"
sed -e "s#^BIDS_DIR=.*#BIDS_DIR=\"$BIDS_DIR\"#" \
    -e "s#^OUTPUT_DIR=.*#OUTPUT_DIR=\"$OUTPUT_DIR\"#" \
    -e "s#^CONFIG_FILE=.*#CONFIG_FILE=\"$CONFIG_FILE\"#" \
    -e "s#^BATCH_SIZE=.*#BATCH_SIZE=100000#" \
    -e "s#^PARALLEL_JOBS=.*#PARALLEL_JOBS=$PARALLEL_JOBS#" \
    -e "s#^CPAC_IMAGE=.*#CPAC_IMAGE=\"$CPAC_IMAGE\"#" \
    -e "s#^DOCKER_PLATFORM=.*#DOCKER_PLATFORM=\"$platform\"#" \
    -e "s#^N_CPUS=.*#N_CPUS=\"$N_CPUS\"#" \
    -e "s#^MEM_GB=.*#MEM_GB=\"$MEM_GB\"#" \
    "$HERE/run_cpapc_adni_ants.sh" > "$RUN"
chmod +x "$RUN"

echo "Starting at $(date). Logs: $OUTPUT_DIR/logs/sub-*.log"
echo "Follow progress in another window with:  tail -f \"$OUTPUT_DIR/processing_log.txt\""
if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -i "$RUN"          # keep the Mac awake while it runs
else
    "$RUN"
fi

ok=$(grep "Successfully completed:" "$OUTPUT_DIR/processing_log.txt" 2>/dev/null | awk '{print $4}' | sort -u | wc -l | tr -d ' ')
echo
echo "Finished at $(date): $ok of $n_sub subjects done."
if [ "$ok" -lt "$n_sub" ]; then
    echo "Some subjects failed (see $OUTPUT_DIR/logs/). Run this script again to retry them."
fi
