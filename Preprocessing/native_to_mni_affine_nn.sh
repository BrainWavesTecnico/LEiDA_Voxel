#!/bin/bash
# Map C-PAC's native-space preprocessed BOLD to MNI space with ONE resampling
# step that does not mix voxel values:
#   - affine only (rotation, translation, scaling, shear; no nonlinear warp)
#   - nearest-neighbour: every MNI voxel takes the value of exactly one
#     native voxel, so no smoothing or averaging is introduced.
#
# Transform chain, applied in a single antsApplyTransforms call:
#   native BOLD --(C-PAC's func->T1 matrix, or an ANTs rigid fit)--> T1
#               --(ANTs Rigid+Affine fit of the T1 brain to MNI152)--> MNI
#
# Runs inside the C-PAC container (it uses ANTs, FSL, AFNI and c3d from it).
# From a Mac, for one subject/session:
#
#   docker run --rm --platform linux/amd64 \
#     -v /path/to/cpac_output_dir:/data \
#     -v ~/LEiDA_Voxel/Preprocessing:/scripts:ro \
#     --entrypoint bash fcpindi/c-pac:release-v1.8.7 \
#     /scripts/native_to_mni_affine_nn.sh /data/output/pipeline_<name>/sub-C01_ses-pre 2
#
# Arguments: <C-PAC subject/session output folder (has anat/ and func/)> [MNI voxel size in mm, default 2]
# Writes to <that folder>/affine_nn/.
set -euo pipefail

SES_DIR=${1:?usage: native_to_mni_affine_nn.sh <cpac sub/ses output folder> [res_mm]}
RES=${2:-2}
ANAT="$SES_DIR/anat"; FUNC="$SES_DIR/func"
OUT="$SES_DIR/affine_nn"; W="$OUT/work"
mkdir -p "$W"

pick() {  # first file matching a pattern, excluding template-space files
    local f
    f=$(ls $1 2>/dev/null | grep -v 'space-' | head -1 || true)
    echo "$f"
}
need() { [ -n "$2" ] && [ -f "$2" ] || { echo "ERROR: no $1 found. Files there:"; ls "$3"; exit 1; }; }

BOLD=$(pick "$FUNC/*_desc-preproc_bold.nii.gz");          need "native preprocessed BOLD (*_desc-preproc_bold.nii.gz)" "$BOLD" "$FUNC"
# T1 brain: C-PAC's *_desc-brain_T1w if present, otherwise the head image
# (*_desc-head_T1w) times the brain mask (*_desc-brain_mask), both native.
T1=$(pick "$ANAT/*_desc-brain_T1w.nii.gz")
if [ -z "$T1" ]; then
    HEAD=$(pick "$ANAT/*_desc-head_T1w.nii.gz");          need "native T1 head (*_desc-head_T1w.nii.gz)" "$HEAD" "$ANAT"
    MASK=$(pick "$ANAT/*_desc-brain_mask.nii.gz");        need "native brain mask (*_desc-brain_mask.nii.gz)" "$MASK" "$ANAT"
    T1="$W/T1w_brain.nii.gz"
    fslmaths "$HEAD" -mas "$MASK" "$T1"
fi
MEAN=$(pick "$FUNC/*_desc-mean_bold.nii.gz"); [ -n "$MEAN" ] || MEAN=$(pick "$FUNC/*_sbref.nii.gz")
MAT=$(pick "$FUNC/*from-bold_to-T1w*xfm.mat")

echo "Native BOLD : $BOLD"
echo "T1 brain    : $T1"
echo "Mean BOLD   : ${MEAN:-none}"
echo "func->T1    : ${MAT:-none (will fit an ANTs rigid transform)}"
TR=$(3dinfo -tr "$BOLD")
echo "TR          : $TR s, output voxel size: ${RES} mm"

TPL="$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz"

# 1. T1 brain -> MNI152 brain: rigid then affine (12 parameters), no SyN.
if [ ! -f "$W/t1_to_mni_0GenericAffine.mat" ]; then
    echo "== Affine registration of the T1 to MNI152"
    antsRegistration --dimensionality 3 --float 0 \
        --output "[$W/t1_to_mni_,$W/t1_in_mni_check.nii.gz]" \
        --interpolation Linear --use-histogram-matching 0 \
        --winsorize-image-intensities "[0.005,0.995]" \
        --initial-moving-transform "[$TPL,$T1,1]" \
        --transform "Rigid[0.1]" --metric "MI[$TPL,$T1,1,32,Regular,0.25]" \
        --convergence "[1000x500x250x100,1e-6,10]" --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
        --transform "Affine[0.1]" --metric "MI[$TPL,$T1,1,32,Regular,0.25]" \
        --convergence "[1000x500x250x100,1e-6,10]" --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox
fi

# 2. func -> T1 as an ITK transform.
if [ -n "$MAT" ]; then
    echo "== Converting C-PAC's func->T1 matrix to ITK"
    fslroi "$BOLD" "$W/bold_vol0" 0 1          # only the grid of this image is used
    c3d_affine_tool -ref "$T1" -src "$W/bold_vol0.nii.gz" "$MAT" -fsl2ras -oitk "$W/func_to_t1.txt"
    FUNC2T1="$W/func_to_t1.txt"
else
    [ -n "$MEAN" ] || { echo "ERROR: no func->T1 matrix and no mean BOLD to fit one. Files:"; ls "$FUNC"; exit 1; }
    echo "== Rigid registration of the mean BOLD to the T1 (ANTs, mutual information)"
    antsRegistration --dimensionality 3 --float 0 \
        --output "[$W/func_to_t1_]" --interpolation Linear \
        --initial-moving-transform "[$T1,$MEAN,1]" \
        --transform "Rigid[0.1]" --metric "MI[$T1,$MEAN,1,32,Regular,0.25]" \
        --convergence "[1000x500x250,1e-6,10]" --shrink-factors 4x2x1 --smoothing-sigmas 2x1x0vox
    FUNC2T1="$W/func_to_t1_0GenericAffine.mat"
fi

# 3. Output grid: MNI152 at the requested voxel size (only the grid is used).
if [ "$RES" = "2" ]; then
    REF="$TPL"
else
    ResampleImageBySpacing 3 "$TPL" "$W/ref_${RES}mm.nii.gz" "$RES" "$RES" "$RES" 0 0 1
    REF="$W/ref_${RES}mm.nii.gz"
fi

# 4. One nearest-neighbour resampling of the whole time series.
NAME=$(basename "$BOLD" _desc-preproc_bold.nii.gz)
DEST="$OUT/${NAME}_space-MNI152affine_res-${RES}mm_desc-preprocNN_bold.nii.gz"
echo "== Applying func->T1->MNI (affine) with nearest neighbour"
antsApplyTransforms -d 3 -e 3 -n NearestNeighbor --float 1 \
    -i "$BOLD" -r "$REF" -o "$DEST" \
    -t "$W/t1_to_mni_0GenericAffine.mat" -t "$FUNC2T1"
3drefit -TR "$TR" "$DEST" >/dev/null 2>&1 || true

# Same transform for the mean BOLD and the T1, to check the alignment.
[ -n "$MEAN" ] && antsApplyTransforms -d 3 -n NearestNeighbor -i "$MEAN" -r "$REF" \
    -o "$OUT/${NAME}_space-MNI152affine_res-${RES}mm_desc-meanNN_bold.nii.gz" \
    -t "$W/t1_to_mni_0GenericAffine.mat" -t "$FUNC2T1"
antsApplyTransforms -d 3 -n NearestNeighbor -i "$T1" -r "$REF" \
    -o "$OUT/${NAME%%_task-*}_space-MNI152affine_res-${RES}mm_desc-brainNN_T1w.nii.gz" \
    -t "$W/t1_to_mni_0GenericAffine.mat"

echo
echo "Done: $DEST"
3dinfo -tr -nt -ni -nj -nk "$DEST" | awk '{print "TR " $1 " s, " $2 " volumes, grid " $3 "x" $4 "x" $5}'
