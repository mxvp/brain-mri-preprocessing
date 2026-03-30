"""Find and prepare brain MRI files from a dataset for preprocessing.

Handles per-dataset file discovery, format conversion, and subject grouping.
Outputs a manifest.json with per-subject center + moving modalities.

Usage:
    python prepare.py ixi data/IXI staging/IXI
    python prepare.py ppmi data/PPMI staging/PPMI
    python prepare.py --list
"""

import argparse
import json
import logging
from pathlib import Path

from datasets import REGISTRY

log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare brain MRI files from a dataset")
    parser.add_argument("dataset", nargs="?", help=f"Dataset name: {', '.join(REGISTRY.keys())}")
    parser.add_argument("input", nargs="?", type=Path, help="Input directory (dataset root)")
    parser.add_argument("output", nargs="?", type=Path, help="Output directory for prepared files")
    parser.add_argument("--list", action="store_true", help="List available datasets")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.list:
        print("Available datasets:")
        for name, cls in REGISTRY.items():
            doc = (cls.__doc__ or "").split("\n")[0].strip()
            print(f"  {name:12s} {doc}")
        raise SystemExit(0)

    if not args.dataset or not args.input or not args.output:
        parser.error("Required: dataset, input, output (or --list)")

    if args.dataset not in REGISTRY:
        parser.error(f"Unknown dataset '{args.dataset}'. Use --list to see options.")

    dataset = REGISTRY[args.dataset]()
    subjects = dataset.prepare(args.input, args.output)

    args.output.mkdir(parents=True, exist_ok=True)
    manifest = args.output / "manifest.json"
    with open(manifest, "w") as f:
        json.dump(subjects, f, indent=2)

    total_volumes = sum(1 + len(s["moving"]) for s in subjects)
    log.info(f"Wrote {len(subjects)} subjects ({total_volumes} volumes) to {manifest}")
    log.info(f"Next: python preprocess.py --manifest {manifest} --output <output_dir>")
