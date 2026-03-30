"""Preprocess brain MRI volumes: N4 + skull-strip + affine registration to SRI24.

Wraps brainles-preprocessing. Supports per-subject multi-modality preprocessing
where T1 is registered to atlas and other modalities (T2, FLAIR) are co-registered
through T1.

Usage:
    python preprocess.py --manifest manifest.json output_dir/
    python preprocess.py input.nii.gz output.nii.gz
    python preprocess.py input.nii.gz output.nii.gz --device cpu
"""

import argparse
import json
import logging
import time
from pathlib import Path

import nibabel as nib
import numpy as np
from brainles_preprocessing.brain_extraction import HDBetExtractor
from brainles_preprocessing.constants import Atlas
from brainles_preprocessing.modality import CenterModality, Modality
from brainles_preprocessing.n4_bias_correction import SitkN4BiasCorrector
from brainles_preprocessing.preprocessor import AtlasCentricPreprocessor
from brainles_preprocessing.registration import ANTsRegistrator


log = logging.getLogger(__name__)


def _clip_negatives(path: Path) -> Path:
    """Clip negative values to 0. Returns path to cleaned file (in-place if needed)."""
    img = nib.load(path)
    data = img.get_fdata()
    if data.min() < 0:
        data = np.clip(data, 0, None)
        cleaned = path.parent / f"_cleaned_{path.name}"
        nib.save(nib.Nifti1Image(data.astype(np.float32), img.affine, img.header), cleaned)
        log.info(f"Clipped negatives in {path.name}")
        return cleaned
    return path


def _output_name(input_path: Path) -> str:
    name = input_path.name
    for ext in (".nii.gz", ".nii"):
        if name.endswith(ext):
            return name[: -len(ext)] + "_preprocessed.nii.gz"
    return name + "_preprocessed.nii.gz"


def preprocess_subject(subject: dict, output_dir: Path, device: str = "0"):
    """Preprocess a single subject with center + moving modalities."""
    subject_id = subject["subject_id"]
    center_info = subject["center"]
    moving_info = subject.get("moving", [])

    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # Check if all outputs already exist
    center_out = output_dir / _output_name(Path(center_info["path"]))
    moving_outs = [output_dir / _output_name(Path(m["path"])) for m in moving_info]
    if center_out.exists() and all(o.exists() for o in moving_outs):
        log.info(f"Skipping {subject_id} (all outputs exist)")
        return

    log.info(f"Processing {subject_id} ({1 + len(moving_info)} modalities)")
    t0 = time.time()

    # Clip negatives
    center_path = _clip_negatives(Path(center_info["path"]))

    center = CenterModality(
        modality_name=center_info["modality"],
        input_path=center_path,
        n4_bias_correction=True,
        raw_bet_output_path=center_out,
    )

    moving = []
    for m_info, m_out in zip(moving_info, moving_outs):
        m_path = _clip_negatives(Path(m_info["path"]))
        moving.append(Modality(
            modality_name=m_info["modality"],
            input_path=m_path,
            n4_bias_correction=True,
            raw_bet_output_path=m_out,
        ))

    preprocessor = AtlasCentricPreprocessor(
        center_modality=center,
        moving_modalities=moving,
        registrator=ANTsRegistrator(
            registration_params={"type_of_transform": "Affine"}
        ),
        brain_extractor=HDBetExtractor(),
        n4_bias_corrector=SitkN4BiasCorrector(),
        atlas_image_path=Atlas.SRI24,
        use_gpu=(device != "cpu"),
    )

    preprocessor.run(log_file=log_dir / f"{subject_id}.log")

    # Clean up temp clipped files
    for p in [Path(center_info["path"])] + [Path(m["path"]) for m in moving_info]:
        cleaned = p.parent / f"_cleaned_{p.name}"
        if cleaned.exists():
            cleaned.unlink()

    elapsed = time.time() - t0
    log.info(f"Done {subject_id} ({elapsed:.1f}s)")


def preprocess_manifest(manifest_path: Path, output_dir: Path, device: str = "0"):
    """Preprocess all subjects from a manifest.json."""
    with open(manifest_path) as f:
        subjects = json.load(f)

    output_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Processing {len(subjects)} subjects from {manifest_path.name}")

    failed = []
    for i, subject in enumerate(subjects):
        log.info(f"[{i+1}/{len(subjects)}] {subject['subject_id']}")
        try:
            preprocess_subject(subject, output_dir, device=device)
        except Exception:
            log.exception(f"Failed: {subject['subject_id']}")
            failed.append(subject["subject_id"])

    if failed:
        log.warning(f"{len(failed)}/{len(subjects)} failed: {failed}")
    else:
        log.info(f"All {len(subjects)} subjects processed successfully")


# Keep single-file mode for quick tests
def preprocess_volume(input_path: Path, output_path: Path, device: str = "0"):
    """Preprocess a single volume (no co-registration)."""
    subject = {
        "subject_id": input_path.stem,
        "center": {"modality": "t1", "path": str(input_path)},
        "moving": [],
    }
    preprocess_subject(subject, output_path.parent, device=device)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Preprocess brain MRI")
    parser.add_argument("input", nargs="?", type=Path, help="Input .nii.gz file")
    parser.add_argument("output", nargs="?", type=Path, help="Output .nii.gz file or directory")
    parser.add_argument("--manifest", type=Path, help="manifest.json from prepare.py")
    parser.add_argument("--device", default="0", help="GPU device ID or 'cpu' (default: 0)")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.manifest:
        if not args.output:
            parser.error("--output required with --manifest")
        preprocess_manifest(args.manifest, args.output, device=args.device)
    elif args.input:
        if not args.output:
            parser.error("output path required")
        preprocess_volume(args.input, args.output, device=args.device)
    else:
        parser.error("Provide --manifest or input file")
