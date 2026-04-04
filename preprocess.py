"""Preprocess brain MRI volumes: N4 + skull-strip + affine registration to SRI24.

Wraps brainles-preprocessing. Supports per-subject multi-modality preprocessing
where T1 is registered to atlas and other modalities (T2, FLAIR) are co-registered
through T1.

Usage:
    python preprocess.py --manifest manifest.json --output output_dir/
    python preprocess.py input.nii.gz --output output.nii.gz
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
from brainles_preprocessing.preprocessor import AtlasCentricPreprocessor, NativeSpacePreprocessor
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


def _com_align(input_path: Path, atlas_enum) -> Path:
    """Create a copy with origin shifted so brain CoM aligns with atlas CoM.
    Returns path to the aligned copy (caller should clean up)."""
    import ants
    import brainles_preprocessing.registration as _reg

    atlas_name = {Atlas.SRI24: "sri24.nii", Atlas.SRI24_SKULLSTRIPPED: "sri24_skullstripped.nii"}
    atlas_paths = list(Path(_reg.__file__).parent.rglob(atlas_name.get(atlas_enum, "sri24.nii")))
    if not atlas_paths:
        log.warning("Atlas not found for CoM alignment, skipping")
        return input_path
    atlas = ants.image_read(str(atlas_paths[0]))
    moving = ants.image_read(str(input_path))

    def _phys_com(img):
        arr = img.numpy()
        thresh = np.percentile(arr[arr > 0], 10) if arr.max() > 0 else 0
        com_vox = np.argwhere(arr > thresh).mean(axis=0)
        direction = np.array(img.direction).reshape(3, 3)
        return np.array(img.origin) + direction @ (com_vox * np.array(img.spacing))

    shift = _phys_com(atlas) - _phys_com(moving)
    new_origin = [moving.origin[i] + shift[i] for i in range(3)]
    moving.set_origin(new_origin)

    aligned = input_path.parent / f"_com_aligned_{input_path.name}"
    ants.image_write(moving, str(aligned))
    log.info(f"CoM aligned {input_path.name} (shift: [{shift[0]:.1f}, {shift[1]:.1f}, {shift[2]:.1f}])")
    return aligned


def _output_name(subject_id: str, modality: str) -> str:
    return f"{subject_id}_{modality}_preprocessed.nii.gz"


def preprocess_subject(subject: dict, output_dir: Path, device: str = "0"):
    """Preprocess a single subject with center + moving modalities."""
    subject_id = subject["subject_id"]
    center_info = subject["center"]
    moving_info = subject.get("moving", [])

    output_dir.mkdir(parents=True, exist_ok=True)
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    # Check if all outputs already exist
    center_out = output_dir / _output_name(subject_id, center_info["modality"])
    moving_outs = [output_dir / _output_name(subject_id, m["modality"]) for m in moving_info]
    if center_out.exists() and all(o.exists() for o in moving_outs):
        log.info(f"Skipping {subject_id} (all outputs exist)")
        return

    log.info(f"Processing {subject_id} ({1 + len(moving_info)} modalities)")
    t0 = time.time()

    pre_registered = subject.get("pre_registered", False)
    pre_skull_stripped = subject.get("pre_skull_stripped", False)
    com_align = subject.get("com_align", False) or pre_skull_stripped

    # Clip negatives
    center_path = _clip_negatives(Path(center_info["path"]))

    # CoM alignment if requested (helps with defaced/cropped inputs)
    if com_align and not pre_skull_stripped:
        atlas_enum = Atlas.SRI24
        center_path = _com_align(center_path, atlas_enum)

    center = CenterModality(
        modality_name=center_info["modality"],
        input_path=center_path,
        n4_bias_correction=True,
        raw_bet_output_path=center_out,
    )

    moving = []
    for m_info, m_out in zip(moving_info, moving_outs):
        m_path = _clip_negatives(Path(m_info["path"]))
        if com_align and not pre_skull_stripped:
            m_path = _com_align(m_path, atlas_enum)
        moving.append(Modality(
            modality_name=m_info["modality"],
            input_path=m_path,
            n4_bias_correction=True,
            raw_bet_output_path=m_out,
        ))

    if pre_registered:
        # Already in SRI24 space (e.g. OASIS) — skip atlas registration
        preprocessor = NativeSpacePreprocessor(
            center_modality=center,
            moving_modalities=moving,
            brain_extractor=None if pre_skull_stripped else HDBetExtractor(),
            n4_bias_corrector=SitkN4BiasCorrector(),
            use_gpu=(device != "cpu"),
        )
    elif pre_skull_stripped:
        # Already skull-stripped but needs atlas registration (e.g. SCHIZO)
        # CoM-align for better registration convergence
        # Use raw_skull_output_path (not raw_bet) to skip brain extraction
        aligned_center = _com_align(center_path, Atlas.SRI24_SKULLSTRIPPED)
        center = CenterModality(
            modality_name=center_info["modality"],
            input_path=aligned_center,
            n4_bias_correction=True,
            raw_skull_output_path=center_out,
        )
        moving = []
        for m_info, m_out in zip(moving_info, moving_outs):
            m_path = _clip_negatives(Path(m_info["path"]))
            aligned_m = _com_align(m_path, Atlas.SRI24_SKULLSTRIPPED)
            moving.append(Modality(
                modality_name=m_info["modality"],
                input_path=aligned_m,
                n4_bias_correction=True,
                raw_skull_output_path=m_out,
            ))
        preprocessor = AtlasCentricPreprocessor(
            center_modality=center,
            moving_modalities=moving,
            registrator=ANTsRegistrator(
                registration_params={"type_of_transform": "Affine"}
            ),
            brain_extractor=None,
            n4_bias_corrector=SitkN4BiasCorrector(),
            atlas_image_path=Atlas.SRI24_SKULLSTRIPPED,
            use_gpu=(device != "cpu"),
        )
    elif subject.get("use_ss_atlas", False):
        # Defaced input — register to SS atlas for better convergence, then skull-strip (e.g. ABIDE)
        preprocessor = AtlasCentricPreprocessor(
            center_modality=center,
            moving_modalities=moving,
            registrator=ANTsRegistrator(
                registration_params={"type_of_transform": "Affine"}
            ),
            brain_extractor=HDBetExtractor(),
            n4_bias_corrector=SitkN4BiasCorrector(),
            atlas_image_path=Atlas.SRI24_SKULLSTRIPPED,
            use_gpu=(device != "cpu"),
        )
    else:
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

    # Clean up temp files (clipped + CoM-aligned)
    for p in [Path(center_info["path"])] + [Path(m["path"]) for m in moving_info]:
        for prefix in ("_cleaned_", "_com_aligned_"):
            tmp = p.parent / f"{prefix}{p.name}"
            if tmp.exists():
                tmp.unlink()

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
    parser.add_argument("--output", "-o", type=Path, help="Output directory or file")
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
            parser.error("--output required")
        preprocess_volume(args.input, args.output, device=args.device)
    else:
        parser.error("Provide --manifest or input file")
