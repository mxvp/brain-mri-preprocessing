"""Find and prepare T1w files from a dataset for preprocessing.

Handles per-dataset file discovery and format conversion (Analyze, DICOM -> NIfTI).
Outputs a file list ready for preprocess.py or preprocess_slurm.sh.

Usage:
    python prepare.py ixi data/IXI staging/IXI
    python prepare.py ppmi data/PPMI/PPMI staging/PPMI
    python prepare.py --list
"""

import argparse
import logging
from pathlib import Path

from datasets import REGISTRY

log = logging.getLogger(__name__)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare T1w files from a dataset")
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
    files = dataset.prepare(args.input, args.output)

    filelist = args.output / "files.txt"
    args.output.mkdir(parents=True, exist_ok=True)
    with open(filelist, "w") as f:
        for p in files:
            f.write(f"{p.resolve()}\n")

    log.info(f"Wrote {len(files)} paths to {filelist}")
    log.info(f"Next: python preprocess.py {args.output} <output_dir> --batch")
    log.info(f"  or: sbatch preprocess_slurm.sh {filelist} <output_dir>")
