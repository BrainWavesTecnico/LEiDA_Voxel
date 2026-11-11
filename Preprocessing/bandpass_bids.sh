#!/bin/bash
# Band-pass filter every BOLD run of a BIDS dataset in native space, BEFORE
# any other processing (so C-PAC's slice-timing correction then runs on the
# filtered data), and write the result as a new BIDS dataset.
#
#   - filter: AFNI 3dTproject -passband <low> <high> (default 0.01-0.1 Hz),
#     with quadratic detrending (-polort 2)
#   - each voxel's temporal mean is added back, so C-PAC can still use the
#     image intensities for brain masking and registration
#   - all other files (T1w, JSON sidecars incl. SliceTiming, TSV) are copied
#     unchanged
#
# Runs inside the C-PAC container (AFNI). From a Mac:
#
#   docker run --rm --platform linux/amd64 \
#     -v /Users/user/Documents/Research/PODelirium:/data \
#     -v ~/LEiDA_Voxel/Preprocessing:/scripts:ro \
#     --entrypoint bash fcpindi/c-pac:release-v1.8.7 \
#     /scripts/bandpass_bids.sh /data/POD_BIDS /data/POD_BIDS_bp
#
# Usage: bandpass_bids.sh <input BIDS folder> <output BIDS folder> [low_hz] [high_hz]
# Already filtered runs are skipped, so it can be re-run safely.
set -euo pipefail

IN=${1:?usage: bandpass_bids.sh <in_bids> <out_bids> [low_hz] [high_hz]}
OUT=${2:?usage: bandpass_bids.sh <in_bids> <out_bids> [low_hz] [high_hz]}
LO=${3:-0.01}
HI=${4:-0.1}
IN=${IN%/}; OUT=${OUT%/}
[ -f "$IN/dataset_description.json" ] || { echo "ERROR: $IN has no dataset_description.json"; exit 1; }
[ "$IN" != "$OUT" ] || { echo "ERROR: output folder must differ from input"; exit 1; }
mkdir -p "$OUT"

# 1. Everything except the BOLD images, copied as is.
(cd "$IN" && find . -type f ! -name '*_bold.nii' ! -name '*_bold.nii.gz') | while read -r f; do
    mkdir -p "$OUT/$(dirname "$f")"
    cp -p "$IN/$f" "$OUT/$f"
done

# 2. BOLD runs: band-pass, add the mean back.
n=0; done_before=0
while read -r f; do
    rel=${f#"$IN"/}
    base=${rel%.gz}; base=${base%.nii}
    dst="$OUT/$base.nii.gz"
    if [ -f "$dst" ]; then done_before=$((done_before + 1)); continue; fi
    mkdir -p "$(dirname "$dst")"
    tmp=$(mktemp -d)
    TR=$(3dinfo -tr "$f")
    echo "== $rel  (TR $TR s, passband $LO-$HI Hz)"
    3dTstat -mean -prefix "$tmp/mean.nii.gz" "$f" >/dev/null 2>&1
    3dTproject -input "$f" -TR "$TR" -passband "$LO" "$HI" -polort 2 -prefix "$tmp/bp.nii.gz" >/dev/null 2>&1
    3dcalc -a "$tmp/bp.nii.gz" -b "$tmp/mean.nii.gz" -expr 'a+b' -float -prefix "$tmp/out.nii.gz" >/dev/null 2>&1
    3drefit -TR "$TR" "$tmp/out.nii.gz" >/dev/null 2>&1 || true
    mv "$tmp/out.nii.gz" "$dst"
    rm -rf "$tmp"
    n=$((n + 1))
done < <(find "$IN" -type f \( -name '*_bold.nii' -o -name '*_bold.nii.gz' \) | sort)

echo
echo "Filtered $n BOLD run(s), $done_before already done. Output: $OUT"
