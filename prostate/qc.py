"""Post-organize QC + canonicalization.

Two passes over the existing organized output:

  1. For every file, compute the in-plane physical FOV (mm) and flag files
     whose FOV is unusually wide as `wide_fov=True`. PI-CAI is multi-center —
     some sites use a body-coil whole-pelvis protocol (~350 mm FOV, prostate
     occupies <15% of the image) while the standard prostate-focused protocol
     is ~150-200 mm FOV. The two are different image distributions and would
     hurt a generative model trained on the mix. Flag is non-destructive so
     downstream consumers can filter at load time.

  2. Optionally rewrite each NIfTI to a single canonical orientation (LPS by
     default, since 3 of 4 sources are already LPS). Updates the affine so
     the physical anatomy stays in place — only the storage convention
     changes.

Run via `python -m prostate qc`. Idempotent: re-runs skip files whose
manifest entry already matches the on-disk header.
"""
from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# Wide-FOV threshold. PI-CAI prostate-focused T2 protocols are ~150-200 mm
# in-plane; whole-pelvis T2 surveys are ~350 mm. 250 mm is a clean midpoint.
# This threshold applies ONLY to T2 — diffusion sequences (ADC/DWI) are
# acquired with body coils and have a natively larger FOV that isn't anomalous.
WIDE_FOV_THRESHOLD_MM = 250.0


def _file_properties(nifti_path: Path) -> dict:
    """Cheap properties: shape, voxel spacing, FOV, orientation, nonzero %."""
    img = nib.load(nifti_path)
    sx, sy, sz = (float(z) for z in img.header.get_zooms()[:3])
    nx, ny, nz_dim = img.shape
    fov_xy_max = float(max(sx * nx, sy * ny))
    d = np.asarray(img.dataobj)
    nonzero_pct = float(100.0 * (d != 0).mean())
    orient = "".join(nib.aff2axcodes(img.affine))
    return {
        "physical_fov_xy_mm": round(fov_xy_max, 1),
        "orient":             orient,
        "nonzero_pct":        round(nonzero_pct, 2),
    }


QC_COLUMNS = ["physical_fov_xy_mm", "orient", "nonzero_pct", "wide_fov"]


def add_qc_columns(manifest: pd.DataFrame, output_root: Path,
                   max_workers: int = 8) -> pd.DataFrame:
    """For every row, compute QC fields and merge them into the manifest.

    Existing QC columns are overwritten so re-running is idempotent.
    """
    manifest = manifest.drop(columns=[c for c in QC_COLUMNS if c in manifest.columns])

    def _row_qc(idx_path):
        idx, rel_path = idx_path
        try:
            props = _file_properties(output_root / rel_path)
            return idx, props, None
        except Exception as e:
            return idx, {}, str(e)

    results: dict = {}
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_row_qc, (i, p))
                   for i, p in manifest["path"].items()]
        for i, fut in enumerate(as_completed(futures), 1):
            idx, props, err = fut.result()
            if err:
                log.error(f"row {idx} ({manifest.at[idx, 'path']}): {err}")
                continue
            results[idx] = props
            if i % 1000 == 0:
                log.info(f"  qc progress: {i}/{len(manifest)}")
    log.info(f"  qc complete: {len(results)}/{len(manifest)} rows processed")

    qc_df = pd.DataFrame.from_dict(results, orient="index")
    return manifest.join(qc_df)


def canonicalize_orientation_inplace(
    manifest: pd.DataFrame,
    output_root: Path,
    target: str = "LPS",
    max_workers: int = 8,
) -> pd.DataFrame:
    """Rewrite NIfTI files whose orientation differs from `target`.

    Uses SimpleITK's DICOMOrient — preserves the physical anatomy by
    reindexing voxels and updating the affine consistently.
    """
    import SimpleITK as sitk

    targets = manifest[manifest["orient"] != target]
    log.info(f"  canonicalizing {len(targets)} / {len(manifest)} files → {target}")

    def _rewrite(rel_path: str) -> tuple[str, str | None]:
        try:
            full = output_root / rel_path
            img = sitk.ReadImage(str(full))
            img = sitk.DICOMOrient(img, target)
            sitk.WriteImage(img, str(full))
            return rel_path, None
        except Exception as e:
            return rel_path, str(e)

    ok = fail = 0
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(_rewrite, p) for p in targets["path"]]
        for i, fut in enumerate(as_completed(futures), 1):
            _, err = fut.result()
            if err is None:
                ok += 1
            else:
                fail += 1
            if i % 500 == 0:
                log.info(f"  rewrite progress: {i}/{len(targets)}  ok={ok} fail={fail}")
    log.info(f"  rewrite complete: ok={ok} fail={fail}")

    # Re-tag orient column for the rows we rewrote
    rewritten_idx = targets.index
    manifest.loc[rewritten_idx, "orient"] = target
    return manifest


def run(output_root: Path,
        target_orientation: str = "LPS",
        canonicalize: bool = True,
        max_workers: int = 8) -> Path:
    """Top-level entry: read manifest, add QC, canonicalize, write back."""
    manifest_path = output_root / "manifest.csv"
    manifest = pd.read_csv(manifest_path)
    log.info(f"loaded manifest: {len(manifest)} rows")

    log.info(f"\nStep 1: QC properties (FOV, orient, nonzero)")
    manifest = add_qc_columns(manifest, output_root, max_workers=max_workers)

    # Drop rows whose QC failed (corrupted files unable to be read). Save the
    # dropped list to dropped.csv for audit; leave the files on disk so the
    # user can inspect / replace.
    corrupted_mask = manifest["physical_fov_xy_mm"].isna()
    if corrupted_mask.any():
        dropped = manifest[corrupted_mask].copy()
        dropped["dropped_reason"] = "qc_unreadable"
        dropped_path = output_root / "dropped.csv"
        if dropped_path.exists():
            prior = pd.read_csv(dropped_path)
            dropped = pd.concat([prior, dropped], ignore_index=True).drop_duplicates("path")
        dropped.to_csv(dropped_path, index=False)
        log.warning(f"  dropping {int(corrupted_mask.sum())} unreadable rows → {dropped_path}")
        manifest = manifest[~corrupted_mask].reset_index(drop=True)

    # Wide-FOV is only meaningful for T2 (where it flags the body-coil
    # whole-pelvis acquisitions). ADC/DWI are natively acquired at wider FOV,
    # so flagging them by the same threshold is misleading.
    manifest["wide_fov"] = (
        (manifest["modality"] == "t2") &
        (manifest["physical_fov_xy_mm"] > WIDE_FOV_THRESHOLD_MM)
    )

    if canonicalize:
        log.info(f"\nStep 2: canonicalize orientation → {target_orientation}")
        manifest = canonicalize_orientation_inplace(
            manifest, output_root, target=target_orientation, max_workers=max_workers)

    manifest.to_csv(manifest_path, index=False)
    log.info(f"\nupdated manifest: {manifest_path}")

    # Summary
    log.info(f"\nwide_fov counts (T2 only — for ADC/DWI use physical_fov_xy_mm):")
    t2 = manifest[manifest["modality"] == "t2"]
    for source, sub in t2.groupby("source"):
        n_wide = int(sub["wide_fov"].sum())
        log.info(f"  {source:<14s} {n_wide:>5d} flagged wide / {len(sub):>5d} T2")
    log.info(f"\norientation totals (after canonicalize):")
    for orient, n in manifest["orient"].value_counts().items():
        log.info(f"  {orient:<5s} {n}")

    return manifest_path
