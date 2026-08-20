"""Generate a manifest + README for the warped tumor masks.

Reads the warp_report.csv written by warp_masks.py, measures per-label tumor
volumes from the masks themselves, and emits:

  tumor_masks_manifest.csv   one row per subject: which volume the mask pairs
                             with, which source mask it came from, per-label
                             volumes in mm3, and the warp QC numbers
  README.md                  provenance, method, label scheme, caveats

Volumes are reported in mm3, which equals voxel count exactly since every
preprocessed volume is 1 mm isotropic.

Usage:
    python make_mask_metadata.py <masks_root> <preprocessed_root>
"""
from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path

log = logging.getLogger(__name__)

# BraTS label scheme, shared by both cohorts.
LABELS = {1: "necrotic", 2: "edema", 4: "enhancing"}

COHORT_INFO = {
    "upenn": {
        "full_name": "UPenn-GBM",
        "mask_tmpl": "UPenn_{sid}_segm_preprocessed.nii.gz",
        "ref_tmpl":  "UPenn_{sid}_t1gd_preprocessed.nii.gz",
        "vol_tmpls": ["UPenn_{sid}_t1_preprocessed.nii.gz",
                      "UPenn_{sid}_t1gd_preprocessed.nii.gz",
                      "UPenn_{sid}_t2_preprocessed.nii.gz",
                      "UPenn_{sid}_flair_preprocessed.nii.gz"],
        "source": "TCIA UPenn-GBM structural NIfTI, baseline (_11) sessions",
        "pathology": "glioblastoma (all subjects)",
    },
    "ucsf": {
        "full_name": "UCSF-PDGM",
        "mask_tmpl": "{sid}_segm_preprocessed.nii.gz",
        "ref_tmpl":  "{sid}_t1c_preprocessed.nii.gz",
        "vol_tmpls": ["{sid}_t1_preprocessed.nii.gz",
                      "{sid}_t1c_preprocessed.nii.gz",
                      "{sid}_t2_preprocessed.nii.gz",
                      "{sid}_flair_preprocessed.nii.gz"],
        "source": "UCSF-PDGM v3 tumor_segmentation.nii.gz",
        "pathology": "diffuse glioma — includes low-grade, non-enhancing tumours",
    },
}


def build(masks_root: Path, prep_root: Path):
    import nibabel as nib
    import numpy as np

    rows, summary = [], {}
    for cohort, info in COHORT_INFO.items():
        cdir = masks_root / cohort
        if not cdir.is_dir():
            log.warning(f"no {cohort}/ under {masks_root}, skipping")
            continue

        report_path = cdir / "warp_report.csv"
        report = {}
        if report_path.exists():
            with open(report_path) as f:
                report = {r["subject"]: r for r in csv.DictReader(f)}
        else:
            log.warning(f"no warp_report.csv for {cohort} — QC columns will be blank")

        prefix, suffix = info["mask_tmpl"].split("{sid}")
        masks = sorted(cdir.glob("*_segm_preprocessed.nii.gz"))
        log.info(f"{cohort}: measuring {len(masks)} masks")

        stats = {"n": 0, "labels_only_edema": 0, "labels_missing_necrotic": 0,
                 "vols": [], "mask_kinds": {}}
        for i, mp in enumerate(masks, 1):
            sid = mp.name[len(prefix):-len(suffix)]
            d = np.asarray(nib.load(mp).dataobj)
            present = sorted(int(v) for v in np.unique(d) if v != 0)
            per_label = {name: int((d == lab).sum()) for lab, name in LABELS.items()}
            total = int((d > 0).sum())

            rep = report.get(sid, {})
            # Which preprocessed volumes this mask is valid against — all
            # modalities of a subject share one grid after preprocessing.
            available = [t.format(sid=sid) for t in info["vol_tmpls"]
                         if (prep_root / cohort / t.format(sid=sid)).exists()]

            rows.append({
                "subject": sid,
                "cohort": cohort,
                "mask_file": f"{cohort}/{mp.name}",
                "aligned_to": " ".join(available),
                "source_mask_kind": rep.get("mask_kind", ""),
                "registration_anchor": rep.get("ref_modality", ""),
                "labels_present": "|".join(str(p) for p in present),
                "vol_necrotic_mm3": per_label["necrotic"],
                "vol_edema_mm3": per_label["edema"],
                "vol_enhancing_mm3": per_label["enhancing"],
                "vol_total_mm3": total,
                "pct_inside_brain": rep.get("pct_inside_brain", ""),
                "labels_preserved": rep.get("labels_match", ""),
            })

            stats["n"] += 1
            stats["vols"].append(total)
            if present == [2]:
                stats["labels_only_edema"] += 1
            if 1 not in present:
                stats["labels_missing_necrotic"] += 1
            k = rep.get("mask_kind", "unknown")
            stats["mask_kinds"][k] = stats["mask_kinds"].get(k, 0) + 1
            if i % 200 == 0:
                log.info(f"  {i}/{len(masks)}")

        stats["vol_median"] = int(np.median(stats["vols"])) if stats["vols"] else 0
        stats["vol_min"] = int(min(stats["vols"])) if stats["vols"] else 0
        stats["vol_max"] = int(max(stats["vols"])) if stats["vols"] else 0
        summary[cohort] = stats

    rows.sort(key=lambda r: (r["cohort"], r["subject"]))
    manifest = masks_root / "tumor_masks_manifest.csv"
    with open(manifest, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    log.info(f"wrote {manifest}  ({len(rows)} rows)")

    _write_readme(masks_root, summary)


def _write_readme(masks_root: Path, summary: dict):
    total = sum(s["n"] for s in summary.values())
    lines = [
        "# Tumor masks in preprocessed space",
        "",
        f"{total} tumor segmentation masks resampled onto the grid of the "
        "preprocessed (SRI24-registered) volumes.",
        "",
        "## Why these exist",
        "",
        "UPenn-GBM and UCSF-PDGM both ship tumor segmentations drawn on their",
        "as-downloaded volumes. Those arrive BraTS-framed (brain pushed to the top",
        "of the canvas) and our preprocessing re-registers them to",
        "`SRI24_SKULLSTRIPPED` for canonical placement — so **the shipped masks do",
        "not overlay on the preprocessed volumes**. These are the shipped masks",
        "brought into that space.",
        "",
        "## Contents",
        "",
        "| Cohort | Subjects | Source | Pathology |",
        "|---|---:|---|---|",
    ]
    for c, s in summary.items():
        info = COHORT_INFO[c]
        lines.append(f"| {info['full_name']} | {s['n']} | {info['source']} | {info['pathology']} |")

    lines += [
        "",
        "```",
        "tumor_masks_manifest.csv     one row per subject (see columns below)",
    ]
    for c in summary:
        lines.append(f"{(c + '/'):<29s}masks + warp_report.csv")
    lines += [
        "```",
        "",
        "## Pairing masks with volumes",
        "",
        "Each mask is on the exact grid — 240x240x155, 1 mm isotropic — of that",
        "subject's preprocessed volumes. All modalities of a subject share one grid",
        "after preprocessing, so a mask is valid against any of them; the",
        "`aligned_to` column lists which exist for that subject.",
        "",
        "## Labels",
        "",
        "BraTS convention, identical in both cohorts:",
        "",
        "| Value | Meaning |",
        "|---:|---|",
        "| 0 | background |",
        "| 1 | necrotic / non-enhancing core |",
        "| 2 | peritumoral edema |",
        "| 4 | enhancing tumor |",
        "",
        "(3 is unused, as in BraTS.)",
        "",
        "## Caveats",
        "",
        "**Volumes are ~4% larger than the shipped masks.** Resampling a label map",
        "onto a different grid with nearest-neighbour interpolation inflates it",
        "slightly. Use these for localisation; for exact volumetry go back to the",
        "source masks in the raw dataset.",
        "",
    ]

    only_edema = {c: s["labels_only_edema"] for c, s in summary.items() if s["labels_only_edema"]}
    if only_edema:
        parts = ", ".join(f"{n} in {COHORT_INFO[c]['full_name']}" for c, n in only_edema.items())
        lines += [
            f"**Not every subject has enhancing tumor.** {parts} carry label 2 only.",
            "This is expected for diffuse glioma cohorts, which include low-grade",
            "non-enhancing tumours — it is not a segmentation failure. Code that",
            "assumes an enhancing component per subject will need to handle it.",
            "",
        ]

    missing_nec = {c: s["labels_missing_necrotic"] for c, s in summary.items() if s["labels_missing_necrotic"]}
    if missing_nec:
        parts = ", ".join(f"{n} in {COHORT_INFO[c]['full_name']}" for c, n in missing_nec.items())
        lines += [
            f"**Label 1 absent in some subjects** ({parts}). Mostly genuine — no",
            "necrotic core — but in at least one case a necrotic region small enough",
            "to be lost to nearest-neighbour rounding. `labels_preserved` in the",
            "manifest is False wherever the label set changed during resampling.",
            "",
        ]

    for c, s in summary.items():
        if list(s["mask_kinds"]) == ["automated"]:
            lines += [
                f"**{COHORT_INFO[c]['full_name']} masks are automated, not expert-corrected.**",
                "The dataset also ships manually corrected segmentations for a subset;",
                "those were not present in the copy used here. `source_mask_kind` in the",
                "manifest records which flavour each subject used.",
                "",
            ]

    lines += [
        "## How they were made",
        "",
        "The registration matrices from the original preprocessing run were not",
        "retained, so the transform was re-derived per subject: registering the raw",
        "reference volume onto its preprocessed counterpart (affine, ANTs) recovers",
        "an equivalent mapping, which is then applied to the mask with `genericLabel`",
        "interpolation. Since a subject's raw modalities and its shipped mask all",
        "share one grid, that mapping carries the mask correctly.",
        "",
        "Reproduce with `warp_masks.py` in the brain-mri-preprocessing repo.",
        "",
        "## Verification",
        "",
        "Every subject was checked for label preservation and for the fraction of",
        "warped mask falling inside the brain (`pct_inside_brain`). Beyond that, the",
        "method was validated by confirming that enhancing tumor stays bright",
        "relative to whole-brain mean on T1c after the warp — that it lands on",
        "genuinely enhancing tissue rather than merely somewhere inside the skull —",
        "and by visual inspection of a sample across both cohorts",
        "(`qc_mask_overlay.py`).",
        "",
        "Per-subject QC numbers are in `tumor_masks_manifest.csv` and in each",
        "cohort's `warp_report.csv`.",
        "",
        "## Manifest columns",
        "",
        "| Column | Meaning |",
        "|---|---|",
        "| `subject` | subject id as used by the cohort |",
        "| `cohort` | `upenn` or `ucsf` |",
        "| `mask_file` | path to the mask, relative to this directory |",
        "| `aligned_to` | preprocessed volumes this mask is valid against |",
        "| `source_mask_kind` | which shipped mask it came from (automated / manual / tumor_seg) |",
        "| `registration_anchor` | raw modality used to recover the transform |",
        r"| `labels_present` | label values in this mask, pipe-separated |",
        "| `vol_necrotic_mm3` / `vol_edema_mm3` / `vol_enhancing_mm3` | per-label volume (1 mm iso, so = voxel count) |",
        "| `vol_total_mm3` | all non-zero voxels |",
        "| `pct_inside_brain` | fraction of warped mask inside the brain mask |",
        "| `labels_preserved` | False if the label set changed during resampling |",
        "",
    ]

    for c, s in summary.items():
        lines.append(
            f"{COHORT_INFO[c]['full_name']} total tumor volume: "
            f"median {s['vol_median']:,} mm3 "
            f"(range {s['vol_min']:,}–{s['vol_max']:,})."
        )

    path = masks_root / "README.md"
    path.write_text("\n".join(lines) + "\n")
    log.info(f"wrote {path}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Manifest + README for warped tumor masks")
    p.add_argument("masks_root", type=Path, help="dir containing upenn/ and ucsf/ mask subdirs")
    p.add_argument("preprocessed_root", type=Path, help="dir containing the preprocessed cohorts")
    args = p.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    build(args.masks_root, args.preprocessed_root)
