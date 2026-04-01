"""Diagnose brain MRI datasets: what modalities, preprocessing state, and atlas space each contains.

Groups files by dataset and preprocessing state, showing a compact overview of what you have.

Usage:
    python diagnosis.py data/ -r
    python diagnosis.py data/ -r --json > diagnosis.json
"""

import argparse
import json
import logging
import re
from collections import defaultdict
from pathlib import Path

import nibabel as nib
import numpy as np

log = logging.getLogger(__name__)

# Known atlas shapes at ~1mm iso
ATLAS_SHAPES = {
    (240, 240, 155): "SRI24",
    (182, 218, 182): "MNI152_1mm",
    (182, 218, 160): "MNI152_1mm_crop",
    (176, 208, 176): "Talairach_T88",
    (91, 109, 91): "MNI152_2mm",
}

MODALITY_PATTERNS = [
    (re.compile(r"(?i)flair"), "FLAIR"),
    (re.compile(r"(?i)t1[_\-]?gd|T1GD|t1Gd|ce[_\-]GD"), "T1Gd"),
    (re.compile(r"(?i)(?<![a-z])t1c(?:_bias)?\.nii|[_\-]t1c[_\.\-]|_T1c_|_T1c\."), "T1c"),
    (re.compile(r"(?i)[_\-]t2[_\-w\.]|^t2[_\-w\.]|(?<![a-z])t2(?:w)?(?=[\._\-]|$)"), "T2"),
    (re.compile(r"(?i)[_\-]t1[_\-w\.]|^t1[_\-w\.]|(?<![a-z])t1(?:w)?(?=[\._\-]|$)|mprage|fspgr|spgr|bravo|ffe|tfe|mpr[\-_]"), "T1"),
]

SEG_PATTERN = re.compile(r"(?i)seg|GlistrBoost|_mask|_label|parenchyma|tumor|fseg")


def _infer_modality(name: str) -> str:
    if SEG_PATTERN.search(name):
        return "SEG"
    for pat, mod in MODALITY_PATTERNS:
        if pat.search(name):
            return mod
    return "?"


def _detect_atlas(shape, zooms) -> str:
    s = tuple(shape[:3])
    if all(0.9 <= float(z) <= 1.1 for z in zooms[:3]):
        if s in ATLAS_SHAPES:
            return ATLAS_SHAPES[s]
    return "native"


def _analyze_file(path: Path) -> dict:
    img = nib.load(path)
    data = img.get_fdata()
    shape = tuple(data.shape)
    zooms = img.header.get_zooms()
    data_sq = data.squeeze()
    nz_frac = np.count_nonzero(data_sq) / data_sq.size

    zooms_3 = [float(z) for z in zooms[:3]]
    is_iso = max(zooms_3) / min(zooms_3) < 1.1 if min(zooms_3) > 0 else False
    is_3d = all(s >= 50 for s in shape[:3]) and all(z <= 3.0 for z in zooms_3)

    # Normalization
    nz = data_sq[data_sq != 0]
    if len(nz) == 0:
        norm = "empty"
    elif abs(float(nz.mean())) < 1 and 0.3 < float(nz.std()) < 2:
        norm = "z-score"
    elif float(nz.min()) >= 0 and float(nz.max()) <= 1.01:
        norm = "[0,1]"
    elif float(nz.min()) < -1 and float(nz.max()) < 10:
        norm = "z-score_shifted"
    elif float(nz.mean()) > 50:
        norm = "raw"
    else:
        norm = "other"

    # Detect format
    if path.name.endswith(".hdr"):
        fmt = "analyze"
    elif path.name.endswith(".mgz"):
        fmt = "mgz"
    else:
        fmt = "nifti"

    return {
        "file": str(path),
        "format": fmt,
        "modality": _infer_modality(path.name),
        "shape": list(shape[:3]),
        "voxel_mm": [round(float(z), 2) for z in zooms[:3]],
        "orientation": "".join(nib.aff2axcodes(img.affine)),
        "atlas": _detect_atlas(shape, zooms),
        "stripped": nz_frac < 0.35,
        "nz_pct": round(nz_frac * 100, 1),
        "is_3d": is_3d,
        "is_iso": is_iso,
        "norm": norm,
        "has_neg": bool(data_sq.min() < 0),
    }


def _classify_variant(d: dict) -> str:
    """Classify a file as raw, preprocessed, or segmentation."""
    if d["modality"] == "SEG":
        return "segmentation"
    if d["atlas"] != "native":
        return "registered"
    if d["stripped"]:
        return "skull-stripped"
    return "raw"


def print_dataset_report(diagnoses: list[dict]):
    """Group by dataset, then by variant (raw/preprocessed), show compact summary."""

    # Group by dataset
    datasets = defaultdict(list)
    for d in diagnoses:
        parts = Path(d["file"]).parts
        try:
            idx = parts.index("data")
            dataset = parts[idx + 1]
        except (ValueError, IndexError):
            dataset = "unknown"
        datasets[dataset].append(d)

    for ds_name in sorted(datasets):
        ds_files = datasets[ds_name]
        print(f"\n{'='*80}")
        print(f"  {ds_name}  ({len(ds_files)} files)")
        print(f"{'='*80}")

        # Sub-group by variant
        variants = defaultdict(list)
        for d in ds_files:
            variants[_classify_variant(d)].append(d)

        for variant in ["raw", "skull-stripped", "registered", "segmentation"]:
            if variant not in variants:
                continue
            vfiles = variants[variant]
            print(f"\n  [{variant.upper()}] ({len(vfiles)} files)")

            # Group by (modality, shape, atlas, norm) to collapse similar files
            groups = defaultdict(list)
            for d in vfiles:
                key = (d["modality"], str(d["shape"]), d["atlas"], d["norm"],
                       d["orientation"], str(d["voxel_mm"]), d["stripped"], d.get("format", "nifti"))
                groups[key].append(d)

            for (mod, shape, atlas, norm, orient, vox, stripped, fmt), gfiles in sorted(groups.items()):
                strip_str = f", {gfiles[0]['nz_pct']}% nz" if stripped else ""
                iso_str = "iso" if gfiles[0]["is_iso"] else "aniso"
                neg_str = ", NEG" if any(f["has_neg"] for f in gfiles) else ""
                td_str = "" if all(f["is_3d"] for f in gfiles) else ", NOT_3D"
                fmt_str = f", {fmt}" if fmt != "nifti" else ""

                print(f"    {mod:<8} {len(gfiles):>3}x  {shape}  @ {vox} mm ({iso_str})  "
                      f"orient={orient}  atlas={atlas}  norm={norm}{strip_str}{neg_str}{td_str}{fmt_str}")

                # Show example paths (max 3)
                for f in gfiles[:3]:
                    print(f"             {f['file']}")
                if len(gfiles) > 3:
                    print(f"             ... +{len(gfiles)-3} more")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Diagnose brain MRI datasets")
    parser.add_argument("path", type=Path, help="Directory to scan")
    parser.add_argument("--recursive", "-r", action="store_true", help="Scan recursively")
    parser.add_argument("--json", action="store_true", help="Output raw JSON")
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING, format="%(message)s")

    if args.path.is_file():
        files = [args.path]
    elif args.recursive:
        files = sorted(
            list(args.path.rglob("*.nii*")) +
            list(args.path.rglob("*.hdr")) +
            list(args.path.rglob("*.mgz"))
        )
    else:
        files = sorted(
            list(args.path.glob("*.nii*")) +
            list(args.path.glob("*.hdr")) +
            list(args.path.glob("*.mgz"))
        )

    # NIfTI + Analyze (.hdr) + MGZ
    files = [f for f in files
             if f.name.endswith(".nii.gz") or f.name.endswith(".nii")
             or f.name.endswith(".hdr") or f.name.endswith(".mgz")]

    if not files:
        print(f"No NIfTI files found in {args.path}")
        raise SystemExit(1)

    print(f"Scanning {len(files)} NIfTI files...")

    diagnoses = []
    for i, f in enumerate(files):
        try:
            diagnoses.append(_analyze_file(f))
        except Exception as e:
            print(f"  FAILED: {f}: {e}")

    if args.json:
        def _fix(obj):
            if isinstance(obj, dict): return {k: _fix(v) for k, v in obj.items()}
            if isinstance(obj, list): return [_fix(v) for v in obj]
            if isinstance(obj, (np.bool_,)): return bool(obj)
            if isinstance(obj, (np.integer,)): return int(obj)
            if isinstance(obj, (np.floating,)): return float(obj)
            return obj
        print(json.dumps(_fix(diagnoses), indent=2))
    else:
        print_dataset_report(diagnoses)
