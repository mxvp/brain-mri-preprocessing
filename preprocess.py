"""Preprocess brain MRI volumes: N4 + skull-strip + affine registration to SRI24.

Wraps brainles-preprocessing.

Usage:
    python preprocess.py input.nii.gz output.nii.gz
    python preprocess.py input.nii.gz output.nii.gz --device cpu
    python preprocess.py input_dir/ output_dir/ --batch
    python preprocess.py --filelist files.txt --output output_dir/
"""

import argparse
import logging
import time
from pathlib import Path

from brainles_preprocessing.brain_extraction import HDBetExtractor
from brainles_preprocessing.constants import Atlas
from brainles_preprocessing.modality import CenterModality
from brainles_preprocessing.n4_bias_correction import SitkN4BiasCorrector
from brainles_preprocessing.preprocessor import AtlasCentricPreprocessor
from brainles_preprocessing.registration import ANTsRegistrator


log = logging.getLogger(__name__)


def preprocess_volume(input_path: Path, output_path: Path, device: str = "0"):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    log.info(f"Processing {input_path.name}")
    t0 = time.time()

    center = CenterModality(
        modality_name="t1",
        input_path=input_path,
        n4_bias_correction=True,
        raw_bet_output_path=output_path,
    )

    preprocessor = AtlasCentricPreprocessor(
        center_modality=center,
        moving_modalities=[],
        registrator=ANTsRegistrator(
            registration_params={"type_of_transform": "Affine"}
        ),
        brain_extractor=HDBetExtractor(),
        n4_bias_corrector=SitkN4BiasCorrector(),
        atlas_image_path=Atlas.SRI24,
        use_gpu=(device != "cpu"),
    )

    log_dir = output_path.parent / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{input_path.stem}.log"

    preprocessor.run(
        log_file=log_file,
    )
    log.info(f"Done {input_path.name} -> {output_path.name} ({time.time() - t0:.1f}s)")


def _output_name(input_path: Path) -> str:
    name = input_path.name
    for ext in (".nii.gz", ".nii"):
        if name.endswith(ext):
            return name[: -len(ext)] + "_preprocessed.nii.gz"
    return name + "_preprocessed.nii.gz"


def preprocess_batch(files: list[Path], output_dir: Path, device: str = "0"):
    output_dir.mkdir(parents=True, exist_ok=True)
    if not files:
        log.warning("No files to process")
        return

    log.info(f"Processing {len(files)} volumes")
    failed = []
    for i, f in enumerate(files):
        out = output_dir / _output_name(f)
        if out.exists():
            log.info(f"[{i+1}/{len(files)}] Skipping {f.name} (already exists)")
            continue
        log.info(f"[{i+1}/{len(files)}] {f.name}")
        try:
            preprocess_volume(f, out, device=device)
        except Exception:
            log.exception(f"Failed: {f.name}")
            failed.append(f.name)

    if failed:
        log.warning(f"{len(failed)}/{len(files)} failed: {failed}")
    else:
        log.info(f"All {len(files)} volumes processed successfully")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess brain MRI for GBM-MAE")
    parser.add_argument("input", nargs="?", type=Path, help="Input .nii.gz file or directory (with --batch)")
    parser.add_argument("output", type=Path, help="Output .nii.gz file or directory")
    parser.add_argument("--batch", action="store_true", help="Process all .nii* files in input dir")
    parser.add_argument("--filelist", type=Path, help="Text file with one .nii.gz path per line")
    parser.add_argument("--device", default="0", help="GPU device ID or 'cpu' (default: 0)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.filelist:
        files = [Path(line.strip()) for line in open(args.filelist) if line.strip()]
        preprocess_batch(files, args.output, device=args.device)
    elif args.batch:
        if not args.input:
            parser.error("input directory required with --batch")
        files = sorted(args.input.glob("*.nii*"))
        preprocess_batch(files, args.output, device=args.device)
    else:
        if not args.input:
            parser.error("input file required")
        preprocess_volume(args.input, args.output, device=args.device)
