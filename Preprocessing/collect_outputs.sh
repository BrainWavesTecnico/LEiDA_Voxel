#!/bin/bash
# Copy the main C-PAC outputs of each scan into a small BIDS-style folder,
# ready to compress and share:
#
#   sub-XX/ses-YY/anat/  T1 in MNI space (skull-stripped) + brain, GM, WM
#                        and CSF masks in MNI space
#   sub-XX/ses-YY/func/  final preprocessed fMRI in MNI space (unsmoothed);
#                        with --smoothed also the smoothed version
#   (plus the matching .json sidecars)
#
#   bash collect_outputs.sh <C-PAC output folder> <destination folder> [--smoothed] [--zip]
#
# <C-PAC output folder> is the one given to run_cpac_all.sh (it contains
# output/pipeline_*/). Only reads it. --zip also makes <destination>.zip.
# Lists any scan with a missing file, and scans in the BIDS input that have
# no output at all.
set -u

SRC=${1:?usage: collect_outputs.sh <C-PAC output folder> <destination folder> [--smoothed] [--zip]}
DST=${2:?usage: collect_outputs.sh <C-PAC output folder> <destination folder> [--smoothed] [--zip]}
SRC=${SRC%/}; DST=${DST%/}
SMOOTHED=0; ZIP=0
for a in "$@"; do
    [ "$a" = --smoothed ] && SMOOTHED=1
    [ "$a" = --zip ] && ZIP=1
done

PIPE=$(ls -d "$SRC"/output/pipeline_* 2>/dev/null | head -1)
[ -n "$PIPE" ] || { echo "No output/pipeline_* folder in $SRC"; exit 1; }
SPACE="space-MNI152NLin6ASym"
mkdir -p "$DST"

# copy file (and its .json sidecar) into dest dir; on a Mac `cp -c` makes a
# clone that uses no extra disk space until the copy is changed
put() {
    local f=$1 d=$2
    mkdir -p "$d"
    cp -c "$f" "$d/" 2>/dev/null || cp "$f" "$d/"
    local j="${f%.nii.gz}.json"
    [ -f "$j" ] && { cp -c "$j" "$d/" 2>/dev/null || cp "$j" "$d/"; }
}
pick() { ls $1 2>/dev/null | head -1; }

n_scan=0; n_ok=0; problems=""
for sesdir in "$PIPE"/sub-*/ses-*/; do
    [ -d "$sesdir" ] || continue
    sesdir=${sesdir%/}
    ses=$(basename "$sesdir"); sub=$(basename "$(dirname "$sesdir")")
    n_scan=$((n_scan+1))
    out="$DST/$sub/$ses"
    missing=""

    for pat in "${SPACE}_desc-preproc_T1w" "${SPACE}_desc-brain_mask" \
               "${SPACE}_label-GM_mask" "${SPACE}_label-WM_mask" "${SPACE}_label-CSF_mask"; do
        f=$(pick "$sesdir/anat/*_${pat}.nii.gz")
        if [ -n "$f" ]; then put "$f" "$out/anat"; else missing="$missing anat:$pat"; fi
    done

    f=$(pick "$sesdir/func/*_${SPACE}_desc-preproc_bold.nii.gz")
    if [ -n "$f" ]; then put "$f" "$out/func"; else missing="$missing func:${SPACE}_desc-preproc_bold"; fi
    if [ $SMOOTHED = 1 ]; then
        f=$(ls "$sesdir"/func/*_${SPACE}_*bold.nii.gz 2>/dev/null | grep -E 'desc-sm|fwhm|smooth' | head -1)
        if [ -n "$f" ]; then put "$f" "$out/func"; else missing="$missing func:smoothed"; fi
    fi

    if [ -z "$missing" ]; then
        n_ok=$((n_ok+1))
    else
        problems="$problems
  $sub $ses: missing$missing"
    fi
done

# Scans in the BIDS input without any output (the BIDS folder used for the run)
nooutput=""
bids=$(grep -m1 '^BIDS_DIR=' "$SRC/run_cpac_batch.sh" 2>/dev/null | cut -d'"' -f2)
if [ -n "$bids" ] && [ -d "$bids" ]; then
    for s in "$bids"/sub-*/ses-*/; do
        s=${s%/}; ses=$(basename "$s"); sub=$(basename "$(dirname "$s")")
        [ -d "$PIPE/$sub/$ses" ] || nooutput="$nooutput $sub/$ses"
    done
fi

cat > "$DST/README.txt" <<EOF
Preprocessed PODelirium data (C-PAC 1.8.7, $(basename "$PIPE" | sed 's/^pipeline_//')).
Collected on $(date).

Per scan (sub-XX/ses-YY):
  anat/*_${SPACE}_desc-preproc_T1w   T1, skull-stripped, in MNI152 space (2 mm)
  anat/*_${SPACE}_desc-brain_mask    brain mask, MNI (2 mm)
  anat/*_${SPACE}_label-GM/WM/CSF_mask  tissue masks, MNI (2 mm)
  func/*_${SPACE}_desc-preproc_bold  final fMRI in MNI space, 3 mm:
       slice-timing correction, motion correction (MCFLIRT), regression of
       24 motion parameters + mean white matter (no CSF), band-pass
       0.01-0.1 Hz, linear interpolation to MNI, no spatial smoothing
EOF
[ $SMOOTHED = 1 ] && echo "  func/*smoothed*: the same, with 5 mm FWHM spatial smoothing" >> "$DST/README.txt"
echo "Scripts and details: https://github.com/BrainWavesTecnico/LEiDA_Voxel/tree/delirium-data-preprocessing/Preprocessing" >> "$DST/README.txt"

echo "Scans found:    $n_scan  (complete: $n_ok)"
echo "Subjects:       $(ls -d "$DST"/sub-* 2>/dev/null | wc -l | tr -d ' ')"
[ -n "$problems" ] && echo "Scans with missing files:$problems"
[ -n "$nooutput" ] && echo "Scans in the BIDS input with no output:$nooutput"
echo "Size:           $(du -sh "$DST" | cut -f1)   ($DST)"

if [ $ZIP = 1 ]; then
    ( cd "$(dirname "$DST")" && zip -qr "$(basename "$DST").zip" "$(basename "$DST")" ) && \
        echo "Zip file:       $DST.zip ($(du -sh "$DST.zip" | cut -f1))"
fi
