"""Warp shipped tumor segmentation masks into our preprocessed space.

UPENN-GBM and UCSF-PDGM both ship expert/automated tumor segmentations drawn on
their as-downloaded volumes. Those volumes are BraTS-framed, and our pipeline
re-registers them to SRI24_SKULLSTRIPPED for canonical placement, so the shipped
masks do NOT overlay on our preprocessed output.

The pipeline discarded its registration matrices (brainles_preprocessing writes
them into a temp dir that is deleted on exit, because preprocess.py never passed
`temp_folder`). Rather than re-run preprocessing, this script recovers the
mapping: both endpoints still exist on disk, so registering the raw center
modality onto its preprocessed counterpart reproduces the composite
(CoM shift + affine) transform, which is then applied to the mask.

Masks are warped with `genericLabel` interpolation — linear would produce
fractional labels at every boundary.

Usage:
    python warp_masks.py upenn <raw_root> <preprocessed_dir> <output_dir>
    python warp_masks.py ucsf  <raw_root> <preprocessed_dir> <output_dir>
    python warp_masks.py upenn <raw_root> <preprocessed_dir> <output_dir> --workers 12
"""
from __future__ import annotations

import argparse
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

log = logging.getLogger(__name__)


# Per-cohort layout. `subject_dir` -> subject_id, and how to find the three
# inputs we need: raw reference image, raw mask, preprocessed reference image.
#
# The reference modality should be the CENTER modality used during
# preprocessing (t1 for both cohorts) — that is the volume that was actually
# registered to the atlas, so recovering its transform is exact. Moving
# modalities were co-registered onto it and share the same output grid, so the
# recovered transform applies to the whole subject.
COHORTS = {
    "upenn": {
        "subject_glob":  "UPENN-GBM-*",
        "subject_id":    lambda d: d.name,
        "raw_ref":       lambda d: d / f"{d.name}_T1.nii.gz",
        "raw_ref_alt":   lambda d: d / f"{d.name}_T1GD.nii.gz",
        # Priority order. UPENN-GBM ships automated segmentations for every
        # subject and expert-corrected ones for a subset; prefer the corrected
        # version wherever it exists. `mask_kind` in the report records which
        # was actually used, so the split is auditable after the fact.
        "raw_masks": lambda d: [
            ("manual",    d / f"{d.name}_segm.nii.gz"),
            ("automated", d / f"{d.name}_automated_approx_segm.nii.gz"),
        ],
        "prep_ref":      lambda sid: f"UPenn_{sid}_t1_preprocessed.nii.gz",
        "prep_ref_alt":  lambda sid: f"UPenn_{sid}_t1gd_preprocessed.nii.gz",
        "out_name":      lambda sid: f"UPenn_{sid}_segm_preprocessed.nii.gz",
    },
    "ucsf": {
        # Dirs are "UCSF-PDGM-0004" on disk; some TCIA copies use a
        # "_nifti" suffix, so match both and strip it in subject_id.
        "subject_glob":  "UCSF-PDGM-*",
        "subject_id":    lambda d: d.name.replace("_nifti", ""),
        "raw_ref":       lambda d: d / f"{d.name.replace('_nifti','')}_T1_bias.nii.gz",
        "raw_ref_alt":   lambda d: d / f"{d.name.replace('_nifti','')}_T1c_bias.nii.gz",
        "raw_masks": lambda d: [
            ("tumor_seg", d / f"{d.name.replace('_nifti','')}_tumor_segmentation.nii.gz"),
        ],
        "prep_ref":      lambda sid: f"{sid}_t1_preprocessed.nii.gz",
        "prep_ref_alt":  lambda sid: f"{sid}_t1c_preprocessed.nii.gz",
        "out_name":      lambda sid: f"{sid}_segm_preprocessed.nii.gz",
    },
}


def _resolve(primary, alt):
    """Prefer the center modality; fall back to the contrast-enhanced one."""
    if primary.exists():
        return primary
    if alt.exists():
        return alt
    return None


def _resolve_mask(spec, subject_dir):
    """First existing mask in the cohort's priority order. Returns (kind, path)."""
    for kind, path in spec["raw_masks"](subject_dir):
        if path.exists():
            return kind, path
    return None, None


def warp_one(subject_dir: Path, prep_dir: Path, out_dir: Path, cohort: str) -> dict:
    """Recover the raw->preprocessed transform and apply it to the mask."""
    import ants
    import numpy as np

    spec = COHORTS[cohort]
    sid = spec["subject_id"](subject_dir)
    out_path = out_dir / spec["out_name"](sid)
    if out_path.exists():
        return {"subject": sid, "status": "skipped (exists)"}

    raw_ref = _resolve(spec["raw_ref"](subject_dir), spec["raw_ref_alt"](subject_dir))
    if raw_ref is None:
        return {"subject": sid, "status": "no raw reference image"}
    mask_kind, raw_mask = _resolve_mask(spec, subject_dir)
    if raw_mask is None:
        return {"subject": sid, "status": "no mask"}

    # Match the preprocessed reference to whichever raw modality we resolved,
    # so we register like-for-like rather than across contrasts.
    is_alt = raw_ref == spec["raw_ref_alt"](subject_dir)
    prep_ref = prep_dir / (spec["prep_ref_alt"](sid) if is_alt else spec["prep_ref"](sid))
    if not prep_ref.exists():
        prep_ref = prep_dir / spec["prep_ref"](sid)
        if not prep_ref.exists():
            return {"subject": sid, "status": "no preprocessed reference"}

    fixed = ants.image_read(str(prep_ref))
    moving = ants.image_read(str(raw_ref))
    reg = ants.registration(fixed=fixed, moving=moving, type_of_transform="Affine")

    mask = ants.image_read(str(raw_mask))
    warped = ants.apply_transforms(
        fixed=fixed, moving=mask,
        transformlist=reg["fwdtransforms"],
        interpolator="genericLabel",
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    ants.image_write(warped, str(out_path))

    # Sanity metrics: labels must survive, and the tumor must land in the brain.
    m_in, m_out = mask.numpy(), warped.numpy()
    labels_in = sorted(int(v) for v in np.unique(m_in) if v != 0)
    labels_out = sorted(int(v) for v in np.unique(m_out) if v != 0)
    tum = m_out > 0
    inside = float(100 * (tum & (fixed.numpy() != 0)).sum() / max(tum.sum(), 1))
    return {
        "subject": sid,
        "status": "ok",
        "mask_kind": mask_kind,
        "ref_modality": raw_ref.name,
        "labels_in": labels_in,
        "labels_out": labels_out,
        "labels_match": labels_in == labels_out,
        "pct_inside_brain": round(inside, 2),
        "vox_in": int((m_in > 0).sum()),
        "vox_out": int(tum.sum()),
    }


def dry_run(cohort: str, raw_root: Path, prep_dir: Path) -> int:
    """Resolve every input path without registering anything.

    Catches the failure modes that would otherwise only surface after a long
    run (or not at all): a subject_glob that matches nothing, mask filenames
    that differ from what this cohort spec expects, post-op subjects with no
    preprocessed counterpart.
    """
    spec = COHORTS[cohort]
    subject_dirs = sorted(
        d for d in raw_root.glob(spec["subject_glob"])
        if d.is_dir() and not d.name.startswith(".")
    )
    if not subject_dirs:
        log.error(f"glob {spec['subject_glob']!r} matched NO dirs under {raw_root}")
        log.error("  the cohort spec's subject_glob likely doesn't match this layout.")
        log.error(f"  what's actually there: {[d.name for d in sorted(raw_root.iterdir())[:5] if d.is_dir()]}")
        return 1

    counts = {"ok": 0, "no raw ref": 0, "no mask": 0, "no preprocessed ref": 0}
    examples = {k: [] for k in counts}
    kind_counts: dict[str, int] = {}
    for d in subject_dirs:
        sid = spec["subject_id"](d)
        raw_ref = _resolve(spec["raw_ref"](d), spec["raw_ref_alt"](d))
        if raw_ref is None:
            counts["no raw ref"] += 1; examples["no raw ref"].append(sid); continue
        mask_kind, mask_path = _resolve_mask(spec, d)
        if mask_path is None:
            counts["no mask"] += 1; examples["no mask"].append(sid); continue
        kind_counts[mask_kind] = kind_counts.get(mask_kind, 0) + 1
        is_alt = raw_ref == spec["raw_ref_alt"](d)
        prep = prep_dir / (spec["prep_ref_alt"](sid) if is_alt else spec["prep_ref"](sid))
        if not prep.exists() and not (prep_dir / spec["prep_ref"](sid)).exists():
            counts["no preprocessed ref"] += 1; examples["no preprocessed ref"].append(sid); continue
        counts["ok"] += 1
        if len(examples["ok"]) < 3:
            examples["ok"].append(f"{sid}  [{raw_ref.name} -> {prep.name}]")

    log.info(f"{cohort} DRY RUN: {len(subject_dirs)} subject dirs under {raw_root}")
    for k, n in counts.items():
        if not n:
            continue
        log.info(f"  {k:<22s} {n:>5d}")
        for e in examples[k][:3]:
            log.info(f"      e.g. {e}")
        if len(examples[k]) > 3:
            log.info(f"      ... and {len(examples[k])-3} more")
    if kind_counts:
        log.info(f"  mask flavours found: {kind_counts}")
    if counts["ok"] == 0:
        log.error("nothing would be processed — fix paths before the real run.")
        return 1
    return 0


def run(cohort: str, raw_root: Path, prep_dir: Path, out_dir: Path, workers: int = 4):
    spec = COHORTS[cohort]
    subject_dirs = sorted(
        d for d in raw_root.glob(spec["subject_glob"])
        if d.is_dir() and not d.name.startswith(".")
    )
    if not subject_dirs:
        raise SystemExit(
            f"glob {spec['subject_glob']!r} matched no dirs under {raw_root} — "
            f"nothing to do. Run with --dry-run to inspect the layout."
        )
    log.info(f"{cohort}: {len(subject_dirs)} subject dirs under {raw_root}")
    log.info(f"  preprocessed: {prep_dir}")
    log.info(f"  output:       {out_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    results, ok, failed = [], 0, 0

    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(warp_one, d, prep_dir, out_dir, cohort): d
            for d in subject_dirs
        }
        for i, fut in enumerate(as_completed(futures), 1):
            try:
                r = fut.result()
            except Exception as e:
                r = {"subject": futures[fut].name, "status": f"error: {e}"}
            results.append(r)
            if r["status"] == "ok":
                ok += 1
                if not r["labels_match"]:
                    log.warning(f"  {r['subject']}: labels {r['labels_in']} -> {r['labels_out']}")
                if r["pct_inside_brain"] < 95:
                    log.warning(f"  {r['subject']}: only {r['pct_inside_brain']}% of mask inside brain")
            elif r["status"].startswith(("no ", "error")):
                failed += 1
                log.warning(f"  {r['subject']}: {r['status']}")
            if i % 50 == 0:
                log.info(f"  progress: {i}/{len(subject_dirs)}  ok={ok} failed={failed}")

    log.info(f"done: ok={ok}  failed={failed}  skipped={len(results)-ok-failed}")

    import csv
    report = out_dir / "warp_report.csv"
    keys = ["subject", "status", "mask_kind", "ref_modality", "labels_match",
            "pct_inside_brain", "vox_in", "vox_out"]
    with open(report, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(results, key=lambda r: r["subject"]))
    log.info(f"wrote {report}")
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Warp shipped tumor masks into preprocessed space")
    p.add_argument("cohort", choices=sorted(COHORTS))
    p.add_argument("raw_root", type=Path, help="Raw dataset root (dir of per-subject dirs)")
    p.add_argument("preprocessed_dir", type=Path, help="Our preprocessed output dir for this cohort")
    p.add_argument("output_dir", type=Path, help="Where to write warped masks")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--dry-run", action="store_true",
                   help="Resolve all input paths and report what would happen, without registering")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
    if args.dry_run:
        raise SystemExit(dry_run(args.cohort, args.raw_root, args.preprocessed_dir))
    run(args.cohort, args.raw_root, args.preprocessed_dir, args.output_dir, args.workers)
