#!/bin/bash
# Convert the PODelirium raw data to BIDS for C-PAC. The raw folder is only
# read, never changed.
#
#   bash pod_to_bids.sh <raw POD folder> <BIDS folder>            # dry run
#   bash pod_to_bids.sh <raw POD folder> <BIDS folder> --go       # do it
#
# Raw layout:  POD/anat/<ID>_T1.nii             (raw T1; _bc/_bf versions ignored:
#                                                C-PAC does its own bias correction)
#              POD/func/<ID>_refMRI_d5.nii       -> session "pre"
#              POD/func/<ID>_post_refMRI_d5.nii  -> session "post"
# BIDS layout: sub-<ID>/ses-<pre|post>/anat/sub-<ID>_ses-<ses>_T1w.nii
#              sub-<ID>/ses-<pre|post>/func/sub-<ID>_ses-<ses>_task-rest_bold.nii
# The T1 is copied into every session that has fMRI (C-PAC needs one per session).
#
# It can be re-run after new raw files arrive: new files are added, files
# that changed in the raw folder (e.g. re-transferred because the first copy
# was incomplete) replace their BIDS copy, and identical files are left as
# they are. Added or replaced files make run_cpac_all.sh process that
# subject again.
#
# Acquisition (from the collaborators): TR = 2 s, 50 axial slices, single-band,
# interleaved ascending with odd slices first; slice 1 is the bottom slice
# (RPI). The slice times below assume the slices are evenly spread over the
# TR (2 s / 50 = 0.04 s per slice). They are written to task-rest_bold.json.
set -u

SRC=${1:?usage: pod_to_bids.sh <raw POD folder> <BIDS folder> [--go]}
DST=${2:?usage: pod_to_bids.sh <raw POD folder> <BIDS folder> [--go]}
SRC=${SRC%/}; DST=${DST%/}
GO=0; [ "${3:-}" = "--go" ] && GO=1

[ -d "$SRC/anat" ] && [ -d "$SRC/func" ] || { echo "No anat/ and func/ in $SRC"; exit 1; }

# Copy src to dst. Returns 0 = copied new, 2 = replaced (dst differed from
# src, e.g. src was re-transferred), 1 = dst already identical. Copied files
# get the current time, so run_cpac_all.sh sees they are newer than the
# subject's outputs and processes the subject again.
# On a Mac (APFS) `cp -c` makes a clone that uses no extra disk space.
copy() {
    local status=0
    if [ -e "$2" ]; then
        cmp -s "$1" "$2" && return 1
        status=2
    fi
    if [ $GO = 1 ]; then
        mkdir -p "$(dirname "$2")"
        rm -f "$2"
        cp -c "$1" "$2" 2>/dev/null || cp "$1" "$2"
        touch "$2"
    fi
    return $status
}

n_new=0; n_had=0; n_upd=0; problems=0
for f in "$SRC"/func/*.nii "$SRC"/func/*.nii.gz; do
    [ -e "$f" ] || continue
    name=$(basename "$f")
    ext=".nii"; case "$name" in *.nii.gz) ext=".nii.gz";; esac
    case "$name" in
        *_post_refMRI_d5$ext) id=${name%_post_refMRI_d5$ext}; ses=post ;;
        *_refMRI_d5$ext)      id=${name%_refMRI_d5$ext};      ses=pre ;;
        *) echo "SKIPPED (name not recognised): func/$name"; problems=$((problems+1)); continue ;;
    esac
    sub=$(echo "$id" | tr -cd '[:alnum:]')   # BIDS labels: letters and digits only
    t1=""
    for c in "$SRC/anat/${id}_T1.nii" "$SRC/anat/${id}_T1.nii.gz"; do [ -e "$c" ] && t1=$c; done
    if [ -z "$t1" ]; then
        echo "SKIPPED (no anat/${id}_T1.nii for it): func/$name"; problems=$((problems+1)); continue
    fi
    t1ext=".nii"; case "$t1" in *.nii.gz) t1ext=".nii.gz";; esac
    base="$DST/sub-$sub/ses-$ses"
    copy "$t1" "$base/anat/sub-${sub}_ses-${ses}_T1w$t1ext"; rt=$?
    [ $rt = 2 ] && echo "UPDATED sub-$sub ses-$ses T1  <-  anat/$(basename "$t1") (differed from the BIDS copy)"
    copy "$f" "$base/func/sub-${sub}_ses-${ses}_task-rest_bold$ext"; rf=$?
    case $rf in
        0) echo "NEW     sub-$sub ses-$ses  <-  func/$name + anat/$(basename "$t1")"; n_new=$((n_new+1)) ;;
        2) echo "UPDATED sub-$sub ses-$ses  <-  func/$name (differed from the BIDS copy)"; n_upd=$((n_upd+1)) ;;
        *) n_had=$((n_had+1)) ;;
    esac
    [ $rt = 2 ] && [ $rf = 1 ] && n_upd=$((n_upd+1))
done

# T1s without any fMRI (not included)
for t1 in "$SRC"/anat/*_T1.nii "$SRC"/anat/*_T1.nii.gz; do
    [ -e "$t1" ] || continue
    id=$(basename "$t1"); id=${id%_T1.nii.gz}; id=${id%_T1.nii}
    ls "$SRC/func/${id}_"*refMRI_d5.nii* >/dev/null 2>&1 || echo "NOTE: anat/$(basename "$t1") has no fMRI, not included"
done

if [ $GO = 1 ]; then
    mkdir -p "$DST"
    cat > "$DST/dataset_description.json" <<'EOF'
{
  "Name": "PODelirium",
  "BIDSVersion": "1.8.0",
  "DatasetType": "raw"
}
EOF
    # Slice times: odd slices (1,3,...,49) at 0, 0.04, ..., 0.96 s, then even
    # slices (2,...,50) at 1.00, ..., 1.96 s. Listed in slice order 1..50.
    cat > "$DST/task-rest_bold.json" <<'EOF'
{
  "TaskName": "rest",
  "RepetitionTime": 2.0,
  "SliceEncodingDirection": "k",
  "SliceTiming": [0.0, 1.0, 0.04, 1.04, 0.08, 1.08, 0.12, 1.12, 0.16, 1.16, 0.2, 1.2, 0.24, 1.24, 0.28, 1.28, 0.32, 1.32, 0.36, 1.36, 0.4, 1.4, 0.44, 1.44, 0.48, 1.48, 0.52, 1.52, 0.56, 1.56, 0.6, 1.6, 0.64, 1.64, 0.68, 1.68, 0.72, 1.72, 0.76, 1.76, 0.8, 1.8, 0.84, 1.84, 0.88, 1.88, 0.92, 1.92, 0.96, 1.96]
}
EOF
fi

echo
echo "fMRI runs: $n_new new, $n_upd updated (changed since the last conversion), $n_had unchanged, $problems skipped"
if [ $GO = 1 ]; then
    echo "BIDS folder: $DST ($(ls -d "$DST"/sub-* 2>/dev/null | wc -l | tr -d ' ') subjects)"
else
    echo "Dry run only: nothing was written. If this looks right, add --go"
fi
