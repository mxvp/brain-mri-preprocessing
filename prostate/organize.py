"""Walk every configured source, materialize each volume as a NIfTI in a
harmonized output layout, emit a manifest with train/val splits.

Output layout (under cfg.paths.output_root):
    t2_axial/<source>_<patient>.nii.gz       (combined T2 axial pool)
    t2_sagittal/<source>_<patient>.nii.gz    (PI-CAI only currently)
    t2_coronal/<source>_<patient>.nii.gz     (PI-CAI only currently)
    adc/<source>_<patient>.nii.gz
    dwi/<source>_<patient>.nii.gz
    manifest.csv

The manifest has one row per volume with: source, patient_id, modality, view,
path (relative to output_root), src_format, n_slices, split ('train'|'val').
"""
from __future__ import annotations

import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from . import convert, sources, splits

log = logging.getLogger(__name__)


def _subdir_for(modality: str, view: str) -> str:
    if modality == "t2":
        return f"t2_{view}"          # t2_axial / t2_sagittal / t2_coronal
    return modality                  # adc / dwi


def organize_all(cfg: dict, limit_per_source: int | None = None) -> Path:
    output_root = Path(cfg["paths"]["output_root"]).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    log.info(f"output root: {output_root}")

    rows: list[dict] = []
    for source_name, source_cfg in cfg["sources"].items():
        if not source_cfg.get("enabled", True):
            log.info(f"{source_name}: disabled")
            continue
        root = Path(source_cfg["root"]).expanduser()
        if not root.exists():
            log.warning(f"{source_name}: root {root} not found, skipping")
            continue

        walker = sources.WALKERS[source_name]
        records = list(walker(root))
        if limit_per_source:
            records = records[:limit_per_source]
        log.info(f"{source_name}: discovered {len(records)} volumes")

        min_slices = int(source_cfg.get("min_slices", 5))
        max_workers = int(cfg.get("workers", 8))
        kept = dropped = failed = 0

        def _materialize_one(rec):
            out_name = f"{rec.source}_{rec.patient_id}.nii.gz"
            subdir = _subdir_for(rec.modality, rec.view)
            dst_path = output_root / subdir / out_name
            try:
                convert.materialize_nifti(rec.src_path, rec.src_format, dst_path)
                slices = convert.n_slices(dst_path)
            except Exception as e:
                return ("failed", rec, dst_path, str(e), 0)
            if slices < min_slices:
                dst_path.unlink(missing_ok=True)
                return ("dropped", rec, dst_path, "", slices)
            return ("kept", rec, dst_path, "", slices)

        # Per-task timeout — SimpleITK occasionally hangs on malformed MetaImages.
        # 60 sec is generous (typical conversion is <1 sec) but bounded.
        task_timeout = float(cfg.get("task_timeout_sec", 60.0))
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            fut_to_rec = {pool.submit(_materialize_one, r): r for r in records}
            for i, fut in enumerate(as_completed(fut_to_rec), 1):
                rec = fut_to_rec[fut]
                try:
                    status, rec, dst_path, err, slices = fut.result(timeout=task_timeout)
                except Exception as e:
                    log.error(f"{rec.src_path.name}: {type(e).__name__} ({str(e)[:60]})")
                    failed += 1
                    continue
                if status == "kept":
                    rows.append({
                        "source":      rec.source,
                        "patient_id":  rec.patient_id,
                        "modality":    rec.modality,
                        "view":        rec.view,
                        "path":        str(dst_path.relative_to(output_root)),
                        "src_format":  rec.src_format,
                        "n_slices":    slices,
                    })
                    kept += 1
                elif status == "dropped":
                    dropped += 1
                else:
                    log.error(f"{rec.src_path.name}: failed ({err})")
                    failed += 1
                if i % 500 == 0:
                    log.info(f"  progress: {i}/{len(records)}  kept={kept} dropped={dropped} failed={failed}")
        log.info(f"  → kept {kept}, dropped {dropped} (< {min_slices} slices), failed {failed}")

    manifest = pd.DataFrame(rows).sort_values(
        ["source", "patient_id", "modality", "view"]).reset_index(drop=True)

    if manifest.empty:
        log.warning("manifest is empty — no volumes written.")
        return output_root / "manifest.csv"

    # Train/val split — patient-grouped, source-stratified
    val_fraction = float(cfg["split"]["val_fraction"])
    seed = int(cfg["split"]["seed"])
    log.info(f"\nSplit: val_fraction={val_fraction}, seed={seed}")
    manifest["split"] = splits.patient_grouped_split(
        manifest, val_fraction=val_fraction, seed=seed)
    splits.verify_no_leakage(manifest, manifest["split"])

    manifest_path = output_root / "manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    log.info(f"\nwrote manifest: {manifest_path}")

    # Summary
    log.info("\nFinal counts per output dir:")
    by_dir = manifest.groupby(["modality", "view"]).size().sort_values(ascending=False)
    for (mod, view), n in by_dir.items():
        log.info(f"  {_subdir_for(mod, view):<14s} {n}")
    log.info(f"\nPer-source totals:")
    for source, n in manifest["source"].value_counts().items():
        log.info(f"  {source:<14s} {n}")
    log.info(f"\nSplit totals:")
    for split_name, n in manifest["split"].value_counts().items():
        log.info(f"  {split_name:<14s} {n}")

    return manifest_path
