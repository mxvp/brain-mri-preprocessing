"""Helper utilities for manifest management and dataset inspection.

Usage:
    python helpers.py split data/staging/ppmi/manifest.json 5
    python helpers.py merge data/staging/ppmi/ manifest_merged.json
    python helpers.py count data/staging/ppmi/manifest.json data/preprocessed/ppmi/
"""

import argparse
import json
import math
import logging
from pathlib import Path

log = logging.getLogger(__name__)


def split_manifest(manifest_path: Path, n_parts: int):
    """Split a manifest.json into N parts for parallel SLURM jobs."""
    with open(manifest_path) as f:
        subjects = json.load(f)

    chunk = math.ceil(len(subjects) / n_parts)
    parent = manifest_path.parent
    stem = manifest_path.stem

    for i in range(n_parts):
        part = subjects[i * chunk : (i + 1) * chunk]
        if not part:
            continue
        out = parent / f"{stem}_part{i + 1}.json"
        with open(out, "w") as f:
            json.dump(part, f, indent=2)
        vols = sum(1 + len(s.get("moving", [])) for s in part)
        print(f"  {out.name}: {len(part)} subjects, {vols} volumes")

    print(f"\nSplit {len(subjects)} subjects into {n_parts} parts")


def merge_manifests(manifest_dir: Path, output_path: Path):
    """Merge all manifest_part*.json files back into one."""
    parts = sorted(manifest_dir.glob("manifest_part*.json"))
    if not parts:
        print(f"No manifest_part*.json found in {manifest_dir}")
        return

    subjects = []
    for p in parts:
        with open(p) as f:
            subjects.extend(json.load(f))

    with open(output_path, "w") as f:
        json.dump(subjects, f, indent=2)

    print(f"Merged {len(parts)} parts -> {output_path.name} ({len(subjects)} subjects)")


def count_progress(manifest_path: Path, output_dir: Path):
    """Count how many subjects/volumes are done vs expected."""
    with open(manifest_path) as f:
        subjects = json.load(f)

    done = 0
    total = 0
    failed = []

    for s in subjects:
        sid = s["subject_id"]
        mods = [s["center"]["modality"]] + [m["modality"] for m in s.get("moving", [])]
        for mod in mods:
            total += 1
            out = output_dir / f"{sid}_{mod}_preprocessed.nii.gz"
            if out.exists():
                done += 1
            else:
                failed.append(f"{sid}_{mod}")

    print(f"{done} / {total} volumes done ({total - done} remaining)")
    if failed and len(failed) <= 20:
        print(f"Missing: {', '.join(failed)}")
    elif failed:
        print(f"Missing: {len(failed)} volumes (too many to list)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocessing helpers")
    sub = parser.add_subparsers(dest="command")

    sp = sub.add_parser("split", help="Split manifest into N parts")
    sp.add_argument("manifest", type=Path)
    sp.add_argument("n_parts", type=int)

    mg = sub.add_parser("merge", help="Merge manifest parts back together")
    mg.add_argument("manifest_dir", type=Path)
    mg.add_argument("output", type=Path)

    ct = sub.add_parser("count", help="Count preprocessing progress")
    ct.add_argument("manifest", type=Path)
    ct.add_argument("output_dir", type=Path)

    args = parser.parse_args()

    if args.command == "split":
        split_manifest(args.manifest, args.n_parts)
    elif args.command == "merge":
        merge_manifests(args.manifest_dir, args.output)
    elif args.command == "count":
        count_progress(args.manifest, args.output_dir)
    else:
        parser.print_help()
