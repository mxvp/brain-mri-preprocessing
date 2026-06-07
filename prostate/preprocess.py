"""Per-patient preprocessing pipeline.

Steps (in order):
    1. N4 bias correction on T2 (and ADC per config)
    2. Prostate gland segmentation on T2 → mask in T2 voxel space
    3. Rigid register ADC and DWI → T2
    4. Crop all volumes + mask around the prostate-mask bbox + margin
    5. Resample to target voxel spacing (PI-CAI default 0.5×0.5×3 mm)
    6. Per-modality normalization
    7. Save as a stacked NIfTI (4 channels: T2, DWI, ADC, mask)

Per-patient runtime: ~5–20 sec without nnU-Net seg, dominated by N4 +
registration. With nnU-Net seg add ~5 sec on GPU, ~1 min on CPU.
"""
from __future__ import annotations

import logging
from pathlib import Path

import nibabel as nib
import numpy as np

from .patients import PatientRecord
from . import segment, registration, normalize as norm

log = logging.getLogger(__name__)


def _n4(path: Path, cfg: dict, mask: "object" = None) -> "object":
    """Run SimpleITK N4 bias field correction. Returns a SimpleITK.Image."""
    import SimpleITK as sitk
    img = sitk.ReadImage(str(path), sitk.sitkFloat32)
    n4 = sitk.N4BiasFieldCorrectionImageFilter()
    n4.SetMaximumNumberOfIterations(list(cfg["n4"]["n_iterations"]))
    shrink = cfg["n4"]["shrink_factor"]
    img_small = sitk.Shrink(img, [shrink] * img.GetDimension())
    mask_small = (sitk.Shrink(mask, [shrink] * img.GetDimension())
                  if mask is not None else
                  sitk.OtsuThreshold(img_small, 0, 1, 200))
    n4.Execute(img_small, mask_small)
    bias = n4.GetLogBiasFieldAsImage(img)
    return img / sitk.Exp(bias)


def _bbox_from_mask(mask: np.ndarray) -> tuple[slice, slice, slice]:
    """Tight bounding box of a binary mask, as slices."""
    coords = np.argwhere(mask > 0)
    if not len(coords):
        # Fallback: full volume
        s = mask.shape
        return slice(0, s[0]), slice(0, s[1]), slice(0, s[2])
    lo = coords.min(axis=0)
    hi = coords.max(axis=0) + 1
    return tuple(slice(int(l), int(h)) for l, h in zip(lo, hi))


def _expand_bbox(bbox, margin_vox, shape):
    """Pad a bbox by per-axis voxel margins, clipped to volume shape."""
    return tuple(
        slice(max(0, b.start - m), min(s, b.stop + m))
        for b, m, s in zip(bbox, margin_vox, shape)
    )


def _resample_to_spacing(sitk_img, target_spacing, interp_label=False):
    """Resample a SimpleITK.Image to target voxel spacing."""
    import SimpleITK as sitk
    in_spacing = np.array(sitk_img.GetSpacing())
    in_size = np.array(sitk_img.GetSize())
    target_spacing = np.array(target_spacing, dtype=float)
    out_size = (in_size * in_spacing / target_spacing).round().astype(int).tolist()
    interp = sitk.sitkNearestNeighbor if interp_label else sitk.sitkLinear
    return sitk.Resample(
        sitk_img, out_size, sitk.Transform(), interp,
        sitk_img.GetOrigin(), target_spacing.tolist(),
        sitk_img.GetDirection(), 0, sitk_img.GetPixelID(),
    )


def preprocess_patient(record: PatientRecord, cfg: dict, output_dir: Path) -> Path:
    """Run the full pipeline for one patient. Returns the output NIfTI path."""
    import SimpleITK as sitk

    out_path = output_dir / f"{record.patient_id}.nii.gz"
    if out_path.exists():
        log.info(f"{record.patient_id}: output exists, skipping")
        return out_path

    if not record.is_complete():
        log.warning(f"{record.patient_id}: missing {record.missing()}, skipping")
        return None

    log.info(f"{record.patient_id}: starting")
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. N4 on T2 (and optionally ADC)
    t2_n4 = _n4(record.t2, cfg) if "t2" in cfg["n4"]["apply_to"] else \
            sitk.ReadImage(str(record.t2), sitk.sitkFloat32)
    adc_in = (_n4(record.adc, cfg) if "adc" in cfg["n4"]["apply_to"] else
              sitk.ReadImage(str(record.adc), sitk.sitkFloat32))
    dwi_in = sitk.ReadImage(str(record.dwi), sitk.sitkFloat32)

    # 2. Prostate gland mask, computed on T2 voxel grid.
    # We write the (possibly N4-corrected) T2 to a temp file so the seg backend
    # can re-read it through the same code path it'd use in production.
    tmp_t2 = output_dir / f"_tmp_{record.patient_id}_t2.nii.gz"
    sitk.WriteImage(t2_n4, str(tmp_t2))
    mask_arr = segment.segment_prostate(tmp_t2, cfg)
    tmp_t2.unlink()

    # SimpleITK image for the mask, sharing T2's geometry
    mask_img = sitk.GetImageFromArray(mask_arr.astype(np.uint8))
    mask_img.CopyInformation(t2_n4)

    # 3. Register DWI and ADC to T2
    adc_reg = registration.register_to_t2(record.adc, record.t2,
                                          cfg["registration"]["interpolator"])
    dwi_reg = registration.register_to_t2(record.dwi, record.t2,
                                          cfg["registration"]["interpolator"])

    # 4. Crop around the prostate mask + margin
    margin_mm = np.array(cfg["crop"]["margin_mm"], dtype=float)
    spacing = np.array(t2_n4.GetSpacing())            # (x, y, z)
    margin_vox = (margin_mm / spacing).round().astype(int).tolist()
    bbox = _expand_bbox(_bbox_from_mask(mask_arr.transpose(2, 1, 0)),
                        margin_vox[::-1], mask_arr.shape)
    # SimpleITK uses (z, y, x) ordering for arrays; convert bbox to ITK index/size
    z_slc, y_slc, x_slc = bbox
    crop_index = [int(x_slc.start), int(y_slc.start), int(z_slc.start)]
    crop_size  = [int(x_slc.stop - x_slc.start),
                  int(y_slc.stop - y_slc.start),
                  int(z_slc.stop - z_slc.start)]

    def _crop(img):
        return sitk.RegionOfInterest(img, crop_size, crop_index)

    t2_c, adc_c, dwi_c, mask_c = (_crop(t2_n4), _crop(adc_reg),
                                  _crop(dwi_reg), _crop(mask_img))

    # 5. Resample all to target spacing
    sp = cfg["resample"]["voxel_spacing"]
    t2_r, adc_r, dwi_r = (_resample_to_spacing(im, sp) for im in (t2_c, adc_c, dwi_c))
    mask_r = _resample_to_spacing(mask_c, sp, interp_label=True)

    # 6. Per-modality normalization
    to_arr = sitk.GetArrayFromImage
    t2_arr  = norm.apply(to_arr(t2_r),  cfg["normalize"]["t2"])
    adc_arr = norm.apply(to_arr(adc_r), cfg["normalize"]["adc"])
    dwi_arr = norm.apply(to_arr(dwi_r), cfg["normalize"]["dwi"])
    mask_arr_r = to_arr(mask_r).astype(np.float32)

    # 7. Save as 4-channel NIfTI ordered [t2, dwi, adc, mask] along the 4th axis
    channels = cfg["output"]["channels"]
    arrs = {"t2": t2_arr, "dwi": dwi_arr, "adc": adc_arr, "mask": mask_arr_r}
    stack = np.stack([arrs[c] for c in channels], axis=-1).astype(np.float32)
    # nibabel expects (x, y, z, c); SimpleITK arrays are (z, y, x) so reorder
    stack = np.moveaxis(stack, [0, 1, 2], [2, 1, 0])

    affine = np.eye(4)
    affine[:3, :3] = np.diag(sp + [1])[:3, :3]
    out_img = nib.Nifti1Image(stack, affine)
    out_img.header.set_xyzt_units("mm", "sec")
    nib.save(out_img, out_path)
    log.info(f"{record.patient_id}: wrote {out_path}  shape={stack.shape}")
    return out_path
