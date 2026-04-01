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


def _is_subject_done(subject: dict, output_dir: Path) -> bool:
    """Check if all outputs for a subject already exist."""
    sid = subject["subject_id"]
    center_mod = subject["center"]["modality"]
    if not (output_dir / f"{sid}_{center_mod}_preprocessed.nii.gz").exists():
        return False
    for m in subject.get("moving", []):
        if not (output_dir / f"{sid}_{m['modality']}_preprocessed.nii.gz").exists():
            return False
    return True


def split_manifest(manifest_path: Path, n_parts: int, output_dir: Path = None):
    """Split a manifest.json into N parts for parallel SLURM jobs.
    If output_dir is given, only includes subjects that aren't done yet."""
    with open(manifest_path) as f:
        subjects = json.load(f)

    if output_dir and output_dir.exists():
        remaining = [s for s in subjects if not _is_subject_done(s, output_dir)]
        print(f"  {len(subjects)} total, {len(subjects) - len(remaining)} done, {len(remaining)} remaining")
        subjects = remaining

    if not subjects:
        print("  Nothing to process!")
        return

    chunk = math.ceil(len(subjects) / n_parts)
    parent = manifest_path.parent

    for i in range(n_parts):
        part = subjects[i * chunk : (i + 1) * chunk]
        if not part:
            continue
        out = parent / f"manifest_part{i + 1}.json"
        with open(out, "w") as f:
            json.dump(part, f, indent=2)
        vols = sum(1 + len(s.get("moving", [])) for s in part)
        print(f"  {out.name}: {len(part)} subjects, {vols} volumes")

    print(f"\nSplit {len(subjects)} remaining subjects into {n_parts} parts")


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

    sp = sub.add_parser("split", help="Split manifest into N parts (only unfinished subjects)")
    sp.add_argument("manifest", type=Path)
    sp.add_argument("n_parts", type=int)
    sp.add_argument("--output-dir", type=Path, help="Preprocessed output dir to check for done subjects")

    mg = sub.add_parser("merge", help="Merge manifest parts back together")
    mg.add_argument("manifest_dir", type=Path)
    mg.add_argument("output", type=Path)

    ct = sub.add_parser("count", help="Count preprocessing progress")
    ct.add_argument("manifest", type=Path)
    ct.add_argument("output_dir", type=Path)

    args = parser.parse_args()

    if args.command == "split":
        split_manifest(args.manifest, args.n_parts, output_dir=args.output_dir)
    elif args.command == "merge":
        merge_manifests(args.manifest_dir, args.output)
    elif args.command == "count":
        count_progress(args.manifest, args.output_dir)
    else:
        parser.print_help()
