"""Prostate gland mask. Single entry point: `segment_prostate(t2_path, cfg)`.

The PI-CAI / community standard is the Radboud nnU-Net checkpoint at
DIAGNijmegen/AbdomenMRUS-prostate-segmentation (Dice ≈ 0.96 internal /
0.82 external). Don't use TotalSegmentator — it under-segments badly
(Dice ≈ 0.15).

This module supports two backends:

  - 'stub'   — returns a fixed-fraction center bbox mask. Lets the rest of
               the pipeline run end-to-end without the 1.5 GB nnU-Net
               checkpoint. Use during scaffolding; replace before training.
  - 'nnunet' — runs the Radboud nnU-Net model. Requires the checkpoint
               path in cfg['segmentation']['nnunet_checkpoint'] and the
               nnunetv2 package installed.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)


def _stub_mask(shape: tuple[int, int, int], fraction: float) -> np.ndarray:
    """Half-FOV (or `fraction`-FOV) box centered on the image. Just a placeholder."""
    mask = np.zeros(shape, dtype=np.uint8)
    lo = [int((1 - fraction) / 2 * s) for s in shape]
    hi = [int((1 + fraction) / 2 * s) for s in shape]
    mask[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]] = 1
    return mask


def segment_prostate(t2_path: Path, cfg: dict) -> np.ndarray:
    """Return a 3D uint8 mask matching T2 in **SimpleITK array order (z, y, x)**.

    Returning in SimpleITK convention keeps downstream consumers (registration,
    crop, resample) uniform — all the heavy lifting in preprocess.py uses
    SimpleITK images, where `GetArrayFromImage` gives (z, y, x).
    """
    backend = cfg["segmentation"]["backend"]
    import SimpleITK as sitk
    t2 = sitk.ReadImage(str(t2_path))
    size_xyz = t2.GetSize()                 # SimpleITK image size is (x, y, z)
    arr_shape_zyx = (size_xyz[2], size_xyz[1], size_xyz[0])

    if backend == "stub":
        frac = float(cfg["segmentation"].get("stub_box_fraction", 0.5))
        log.warning(f"Using STUB prostate mask ({frac*100:.0f}% center box). "
                    f"Replace with the Radboud nnU-Net checkpoint before training.")
        return _stub_mask(arr_shape_zyx, frac)
    if backend == "nnunet":
        return _nnunet_mask(t2_path, cfg)
    raise ValueError(f"Unknown segmentation backend: {backend!r}")


_NNUNET_PREDICTOR = None  # cached predictor; nnU-Net init is slow


def _get_nnunet_predictor(model_dir: Path, cfg: dict):
    """Lazy-init nnUNetv2 predictor; cached across patients in a run."""
    global _NNUNET_PREDICTOR
    if _NNUNET_PREDICTOR is not None:
        return _NNUNET_PREDICTOR

    import torch
    from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log.info(f"nnU-Net device: {device}")

    folds = cfg["segmentation"].get("nnunet_folds", "all")
    ckpt_name = cfg["segmentation"].get("nnunet_checkpoint_name", "checkpoint_final.pth")

    p = nnUNetPredictor(
        tile_step_size=0.5,
        use_gaussian=True,
        use_mirroring=True,
        perform_everything_on_device=(device.type == "cuda"),
        device=device,
        verbose=False,
        allow_tqdm=False,
    )
    p.initialize_from_trained_model_folder(
        str(model_dir),
        use_folds=tuple(folds) if isinstance(folds, list) else (folds,),
        checkpoint_name=ckpt_name,
    )
    _NNUNET_PREDICTOR = p
    return p


def _nnunet_mask(t2_path: Path, cfg: dict) -> np.ndarray:
    """Run the Radboud nnU-Net prostate-gland model on a T2w NIfTI.

    Returns the binary gland mask in SimpleITK array order (z, y, x).

    Requires:
      - `nnunetv2` installed (`uv pip install nnunetv2`)
      - The Radboud model directory at cfg['segmentation']['nnunet_checkpoint'].
        See README for how to get the weights from their Docker container.
    """
    ckpt = cfg["segmentation"].get("nnunet_checkpoint")
    if not ckpt:
        raise RuntimeError(
            "segmentation.backend = 'nnunet' but no nnunet_checkpoint path "
            "configured. Set segmentation.nnunet_checkpoint in preprocess.yaml."
        )
    model_dir = Path(ckpt)
    if not model_dir.exists():
        raise FileNotFoundError(f"nnU-Net model dir not found: {model_dir}")

    predictor = _get_nnunet_predictor(model_dir, cfg)

    import SimpleITK as sitk
    # nnU-Net v2's image-reader path is the easiest to drive: it reads a
    # NIfTI, handles spacing/orientation, and gives us a mask back as a
    # numpy array in (z, y, x) order — same convention as the rest of the
    # pipeline expects.
    img = sitk.ReadImage(str(t2_path), sitk.sitkFloat32)
    arr = sitk.GetArrayFromImage(img)[None, ...]   # (C=1, z, y, x)
    spacing_zyx = list(img.GetSpacing())[::-1]     # SITK gives (x,y,z)
    properties = {
        "sitk_stuff": {
            "spacing": img.GetSpacing(),
            "origin":  img.GetOrigin(),
            "direction": img.GetDirection(),
        },
        "spacing": spacing_zyx,
        "shape_after_cropping_and_before_resampling": arr.shape[1:],
    }

    pred = predictor.predict_single_npy_array(
        input_image=arr.astype("float32"),
        image_properties=properties,
        segmentation_previous_stage=None,
        output_file_truncated=None,
        save_or_return_probabilities=False,
    )
    # Radboud model is binary {0 = background, 1 = prostate}
    return (pred > 0).astype(np.uint8)
