"""Canonical PI-CAI preprocessing: resample + center-crop.

After this step every volume is at the same voxel grid (0.5 × 0.5 × 3.0 mm,
160 × 160 × 20 voxels = 80 × 80 × 60 mm physical FOV, prostate-centered).
This is the same preprocessing the official picai_prep tool applies before
training a downstream model — it normalises away the per-center FOV / matrix
variance that's intrinsic to the raw acquisitions.

  - PI-CAI tight-FOV T2 (~192 mm FOV, 384²)     → center-cropped to 80×80
  - PI-CAI wide-FOV T2  (~350 mm FOV, 1024²)    → resampled then cropped to 80×80
  - ADC / DWI (~250+ mm FOV, native)            → resampled then cropped to 80×80
  - PROMISE12 / Prostate158 / TCIA-Biopsy       → same treatment
  - Smaller-than-target inputs (rare)            → zero-padded out to 80×80×60

Input: the `<output_root>` produced by `prostate organize` (with manifest).
Output: a sibling dir `<preprocessed_root>` with the same flat layout, plus
its own `manifest.csv` (the same rows as input, paths updated to point at
the preprocessed files).
"""
from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# picai_prep canonical config (matrix is z, y, x).
DEFAULT_SPACING_XYZ      = (0.5, 0.5, 3.0)
DEFAULT_MATRIX_SIZE_XYZ  = (160, 160, 20)


def _resample_to_spacing(img, target_spacing_xyz):
    """SimpleITK resample to a uniform voxel spacing. Linear interpolation."""
    import SimpleITK as sitk
    in_spacing = np.array(img.GetSpacing())
    in_size    = np.array(img.GetSize())
    out_spacing = np.array(target_spacing_xyz, dtype=float)
    out_size = np.maximum(
        (in_size * in_spacing / out_spacing).round().astype(int),
        1,
    ).tolist()
    return sitk.Resample(
        img, out_size, sitk.Transform(), sitk.sitkLinear,
        img.GetOrigin(), out_spacing.tolist(),
        img.GetDirection(), 0, img.GetPixelID(),
    )


def _center_crop_or_pad(img, target_size_xyz):
    """Center-crop the input to target_size; zero-pad axes that fall short."""
    import SimpleITK as sitk
    in_size = list(img.GetSize())          # (x, y, z)
    target  = list(target_size_xyz)        # (x, y, z)

    # Pad first if any dim is smaller than target.
    pad_lo = [0, 0, 0]; pad_hi = [0, 0, 0]
    for i, (s, t) in enumerate(zip(in_size, target)):
        if s < t:
            total = t - s
            pad_lo[i] = total // 2
            pad_hi[i] = total - pad_lo[i]
    if any(pad_lo) or any(pad_hi):
        img = sitk.ConstantPad(img, pad_lo, pad_hi, 0)
        in_size = list(img.GetSize())

    # Now every dim ≥ target — center-crop.
    crop_lo = [(s - t) // 2 for s, t in zip(in_size, target)]
    return sitk.RegionOfInterest(img, target, crop_lo)


def preprocess_one(
    src_path: Path,
    dst_path: Path,
    spacing_xyz: tuple,
    matrix_size_xyz: tuple,
) -> None:
    """Resample → center crop/pad → write. Idempotent (skip if dst exists)."""
    if dst_path.exists():
        return
    import SimpleITK as sitk
    img = sitk.ReadImage(str(src_path), sitk.sitkFloat32)
    img = _resample_to_spacing(img, spacing_xyz)
    img = _center_crop_or_pad(img, matrix_size_xyz)
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(img, str(dst_path))


def run(
    organized_root: Path,
    preprocessed_root: Path,
    spacing_xyz: tuple = DEFAULT_SPACING_XYZ,
    matrix_size_xyz: tuple = DEFAULT_MATRIX_SIZE_XYZ,
    max_workers: int = 8,
    filters: dict | None = None,
) -> Path:
    """Run the full preprocessing pass and write a parallel manifest."""
    organized_root    = Path(organized_root)
    preprocessed_root = Path(preprocessed_root)
    in_manifest = pd.read_csv(organized_root / "manifest.csv")
    filters = filters or {}
    log.info(f"loaded {len(in_manifest)} rows from {organized_root}/manifest.csv")
    log.info(f"target spacing: {spacing_xyz}   matrix: {matrix_size_xyz}")
    log.info(f"output root:    {preprocessed_root}")

    # Drop wide-FOV rows: naive geometric center crop assumes the prostate
    # is at the image center, which holds for prostate-focused acquisitions
    # but NOT for body-coil whole-pelvis scans (the prostate sits low-
    # posterior in the body, off-center in the image). These rows are flagged
    # `wide_fov=True` by qc.py — only set for T2 modality where the bias is
    # severe (~350 mm native FOV).
    if "wide_fov" in in_manifest.columns:
        n_wide = int(in_manifest["wide_fov"].fillna(False).astype(bool).sum())
        if n_wide > 0:
            log.info(f"dropping {n_wide} rows flagged wide_fov "
                     f"(prostate off-center → center crop misses it)")
            in_manifest = in_manifest[~in_manifest["wide_fov"].fillna(False).astype(bool)].reset_index(drop=True)
            log.info(f"  remaining: {len(in_manifest)} rows to preprocess")

    # Optional subset filters — keep only rows matching the requested
    # source / modality / view values. Each filter is a list; empty/None = no filter.
    for col, values in [
        ("source",   filters.get("sources")),
        ("modality", filters.get("modalities")),
        ("view",     filters.get("views")),
    ]:
        if values:
            before = len(in_manifest)
            in_manifest = in_manifest[in_manifest[col].isin(values)].reset_index(drop=True)
            log.info(f"  filter {col} in {values}: {before} → {len(in_manifest)}")

    preprocessed_root.mkdir(parents=True, exist_ok=True)
    ok = fail = 0
    failures: list[dict] = []

    def _one(row):
        src = organized_root / row["path"]
        dst = preprocessed_root / row["path"]
        try:
            preprocess_one(src, dst, spacing_xyz, matrix_size_xyz)
            return row.name, True, None
        except Exception as e:
            return row.name, False, str(e)

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_one, row) for _, row in in_manifest.iterrows()]
        for i, fut in enumerate(as_completed(futures), 1):
            idx, success, err = fut.result()
            if success:
                ok += 1
            else:
                fail += 1
                failures.append({"row_idx": idx, "error": err,
                                 "path": in_manifest.at[idx, "path"]})
                log.error(f"row {idx} ({in_manifest.at[idx, 'path']}): {err[:80]}")
            if i % 500 == 0:
                log.info(f"  progress: {i}/{len(in_manifest)}  ok={ok} fail={fail}")
    log.info(f"  done: ok={ok}  fail={fail}")

    # Manifest for the preprocessed dir — same rows, paths now relative to
    # preprocessed_root (which is identical to organized_root). Drop failed.
    failed_idx = {f["row_idx"] for f in failures}
    out_manifest = in_manifest[~in_manifest.index.isin(failed_idx)].copy().reset_index(drop=True)
    out_manifest.to_csv(preprocessed_root / "manifest.csv", index=False)
    log.info(f"wrote manifest: {preprocessed_root / 'manifest.csv'}  ({len(out_manifest)} rows)")

    if failures:
        pd.DataFrame(failures).to_csv(preprocessed_root / "dropped.csv", index=False)
        log.warning(f"wrote dropped.csv ({len(failures)} failures)")

    return preprocessed_root / "manifest.csv"
