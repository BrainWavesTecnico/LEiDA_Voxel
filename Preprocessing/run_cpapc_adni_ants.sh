#!/bin/bash
# Configuration
BIDS_DIR="/home/vmiguel/Desktop/ADNI_ZIPS/fmri/BIDS_NII"
OUTPUT_DIR="/home/vmiguel/Desktop/cpac_outputs_ants_mean_v4"
CONFIG_FILE="/home/vmiguel/Desktop/pipelines/cpac_pipeline_with_mean_best_for_csf_0.01_0.1.yml"
BATCH_SIZE=100
PARALLEL_JOBS=3
LOG_FILE="$OUTPUT_DIR/processing_log.txt"

# Create output directory
mkdir -p $OUTPUT_DIR
docker run --rm \
  -v $OUTPUT_DIR:/outputs \
  alpine:latest \
  sh -c "mkdir -p /outputs/output /outputs/working /outputs/logs /outputs/crash && chown -R $(id -u):$(id -g) /outputs"

# Auto-detect already processed subjects
if [ -f "$LOG_FILE" ]; then
    DONE=($(grep "✓ Successfully completed:" "$LOG_FILE" | awk '{print $4}' | sort -u))
else
    DONE=()
fi

# Get all subjects from BIDS
all_subjects=($(ls -d $BIDS_DIR/sub-* | xargs -n 1 basename | sed 's/sub-//'))

# Filter out already done
SUBJECTS=()
for sub in "${all_subjects[@]}"; do
    skip=false
    for done_sub in "${DONE[@]}"; do
        if [ "$sub" == "$done_sub" ]; then
            skip=true
            break
        fi
    done
    if [ "$skip" == false ]; then
        SUBJECTS+=("$sub")
    fi
done

# Limit to batch size
if [ ${#SUBJECTS[@]} -gt $BATCH_SIZE ]; then
    SUBJECTS=("${SUBJECTS[@]:0:$BATCH_SIZE}")
fi

total_subjects=${#SUBJECTS[@]}
total_done=${#DONE[@]}
total_all=${#all_subjects[@]}

echo "================================================" | tee -a "$LOG_FILE"
echo "Already processed: $total_done / $total_all" | tee -a "$LOG_FILE"
echo "This batch: $total_subjects subjects" | tee -a "$LOG_FILE"
echo "Remaining after this batch: $((total_all - total_done - total_subjects))" | tee -a "$LOG_FILE"
echo "Parallel jobs: $PARALLEL_JOBS" | tee -a "$LOG_FILE"
echo "Started at: $(date)" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"

if [ $total_subjects -eq 0 ]; then
    echo "No new subjects to process. All done!" | tee -a "$LOG_FILE"
    exit 0
fi

# Function to process a single subject
process_subject() {
    local subject=$1
    local subject_num=$2
    local subject_log="$OUTPUT_DIR/logs/sub-${subject}.log"

    echo "" | tee -a "$subject_log"
    echo "================================================" | tee -a "$subject_log"
    echo "Processing subject $subject_num/$total_subjects: $subject" | tee -a "$subject_log"
    echo "Date: $(date)" | tee -a "$subject_log"
    echo "================================================" | tee -a "$subject_log"

    docker run --rm \
      -v $BIDS_DIR:/bids:ro \
      -v $OUTPUT_DIR:/outputs \
      -v $CONFIG_FILE:/pipeline_config.yml:ro \
      fcpindi/c-pac:latest \
      /bids /outputs participant \
      --pipeline_file /pipeline_config.yml \
      --participant_label $subject >> "$subject_log" 2>&1

    if [ $? -eq 0 ]; then
        echo "✓ Successfully completed: $subject at $(date)" | tee -a "$subject_log" >> "$LOG_FILE"
        return 0
    else
        echo "✗ Failed: $subject at $(date)" | tee -a "$subject_log" >> "$LOG_FILE"
        return 1
    fi
}

# Process subjects with parallelization
for ((i=0; i<total_subjects; i++)); do
    subject="${SUBJECTS[$i]}"

    process_subject "$subject" "$((i+1))" &

    job_count=$(jobs -r | wc -l)

    while [ $job_count -ge $PARALLEL_JOBS ]; do
        sleep 5
        job_count=$(jobs -r | wc -l)
    done

    echo "Launched: sub-$subject ($((i+1))/$total_subjects) - Active jobs: $job_count" | tee -a "$LOG_FILE"
done

echo "Waiting for all jobs to complete..." | tee -a "$LOG_FILE"
wait
echo "All parallel jobs completed!" | tee -a "$LOG_FILE"

# Fix permissions
docker run --rm \
  -v $OUTPUT_DIR:/outputs \
  alpine:latest \
  chown -R $(id -u):$(id -g) /outputs >> "$LOG_FILE" 2>&1

# Summary
new_done=($(grep "✓ Successfully completed:" "$LOG_FILE" | awk '{print $4}'))
remaining=$((total_all - ${#new_done[@]}))
echo "Total processed: ${#new_done[@]} / $total_all"

echo "" | tee -a "$LOG_FILE"
echo "================================================" | tee -a "$LOG_FILE"
echo "Batch complete at: $(date)" | tee -a "$LOG_FILE"
echo "Total processed: ${#new_done[@]} / $total_all" | tee -a "$LOG_FILE"
echo "Remaining: $remaining" | tee -a "$LOG_FILE"
if [ $remaining -gt 0 ]; then
    echo "Run this script again for the next batch." | tee -a "$LOG_FILE"
else
    echo "All subjects completed!" | tee -a "$LOG_FILE"
fi
echo "================================================" | tee -a "$LOG_FILE"