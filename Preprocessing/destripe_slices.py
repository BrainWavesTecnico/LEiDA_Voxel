#!/usr/bin/env python
"""Remove the alternating odd/even slice artefact from a native-space 4D BOLD.

With interleaved slice acquisition, head motion within a TR (e.g. breathing-
related nodding) gives odd and even slices different intensity, changing
from volume to volume. It is uniform within each slice and alternates along
the slice axis, with a strength that can vary (and change sign) along it.
After normalization it shows up as stripes across MNI slices.

Method (per volume, on the voxel-demeaned data, so voxel means are kept):
  1. m[z, t]: mean signal of each slice z inside the brain mask.
  2. For each slice, fit the 5 nearest slices with
         m = c0 + c1 * (z' - z) + a * (-1)^z'
     (a local linear trend across slices plus an alternating term).
     a * (-1)^z is the artefact of slice z at that volume.
  3. Subtract it uniformly from every brain voxel of slice z.
Only the slice-uniform alternating component is removed: patterns within a
slice and signals that vary smoothly across slices are kept.

Usage:
  destripe_slices.py <in_bold.nii.gz> <out_bold.nii.gz> [--report <prefix>]
Writes <prefix>_artefact.tsv (slices x volumes, the subtracted signal) and
prints how much variance was removed.
"""
import sys
import numpy as np
import nibabel as nib

WINDOW = 5          # slices used in each local fit
MIN_VOXELS = 50     # slices with fewer brain voxels are left unchanged


def destripe(data, mask):
    """data: X x Y x Z x T float array; mask: X x Y x Z bool.
    Returns (cleaned data, artefact Z x T, list of corrected slices)."""
    X, Y, Z, T = data.shape
    demeaned = data - data.mean(axis=3, keepdims=True)
    nvox = mask.reshape(-1, Z).sum(axis=0)
    valid = np.where(nvox >= MIN_VOXELS)[0]

    m = np.zeros((Z, T))
    for z in valid:
        m[z] = demeaned[:, :, z, :][mask[:, :, z]].mean(axis=0)

    art = np.zeros((Z, T))
    if len(valid) < WINDOW:
        return data, art, []
    for i, z in enumerate(valid):
        lo = min(max(i - WINDOW // 2, 0), len(valid) - WINDOW)
        win = valid[lo:lo + WINDOW]
        design = np.column_stack([np.ones(len(win)), win - z, (-1.0) ** win])
        coef = np.linalg.lstsq(design, m[win], rcond=None)[0]   # 3 x T
        art[z] = coef[2] * (-1.0) ** z

    cleaned = data.copy()
    for z in valid:
        sl = cleaned[:, :, z, :]
        sl[mask[:, :, z]] -= art[z]
    return cleaned, art, list(valid)


def main(argv):
    if len(argv) < 3:
        sys.exit(__doc__)
    src, dst = argv[1], argv[2]
    prefix = argv[argv.index('--report') + 1] if '--report' in argv else None

    img = nib.load(src)
    data = np.asarray(img.dataobj, dtype=np.float64)
    if data.ndim != 4:
        sys.exit('Expected a 4D image: %s' % src)
    mask = data.std(axis=3) > 0

    cleaned, art, valid = destripe(data, mask)

    # How much was removed
    Z = data.shape[2]
    sm_before = np.zeros((Z, data.shape[3])); sm_after = np.zeros_like(sm_before)
    d0 = data - data.mean(axis=3, keepdims=True); d1 = cleaned - cleaned.mean(axis=3, keepdims=True)
    for z in valid:
        sm_before[z] = d0[:, :, z, :][mask[:, :, z]].mean(axis=0)
        sm_after[z] = d1[:, :, z, :][mask[:, :, z]].mean(axis=0)
    v = np.array(valid)

    def neighbour_r(s):
        r = [np.corrcoef(s[a], s[b])[0, 1] for a, b in zip(v[:-1], v[1:]) if b == a + 1]
        r2 = [np.corrcoef(s[a], s[b])[0, 1] for a, b in zip(v[:-2], v[2:]) if b == a + 2]
        return np.nanmean(r), np.nanmean(r2)

    vox_var = d0[mask].var(axis=1).sum()
    removed_vox = 1 - d1[mask].var(axis=1).sum() / vox_var
    removed_slice = 1 - sm_after[v].var() / sm_before[v].var()
    rb, rb2 = neighbour_r(sm_before); ra, ra2 = neighbour_r(sm_after)
    print('Slices corrected: %d of %d (others have < %d brain voxels)' % (len(valid), Z, MIN_VOXELS))
    print('Variance removed: %.1f%% of all brain-voxel variance, %.1f%% of slice-mean variance'
          % (100 * removed_vox, 100 * removed_slice))
    print('Slice-mean correlation, adjacent / two apart: before %.2f / %.2f, after %.2f / %.2f'
          % (rb, rb2, ra, ra2))

    out = nib.Nifti1Image(cleaned.astype(np.float32), img.affine, img.header)
    out.set_data_dtype(np.float32)
    out.to_filename(dst)
    if prefix:
        np.savetxt(prefix + '_artefact.tsv', art, delimiter='\t', fmt='%.6g',
                   header='odd/even slice artefact subtracted: rows = slices (z, from 0), columns = volumes')


if __name__ == '__main__':
    main(sys.argv)
