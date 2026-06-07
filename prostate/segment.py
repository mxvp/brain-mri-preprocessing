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


def _nnunet_mask(t2_path: Path, cfg: dict) -> np.ndarray:
    """Run the Radboud nnU-Net prostate-gland model on a T2w volume.

    Stub implementation — wire up to your local nnU-Net install. The
    expected interface is:

        from nnunetv2.inference.predict_from_raw_data import nnUNetPredictor
        predictor = nnUNetPredictor(...)
        predictor.initialize_from_trained_model_folder(checkpoint, ...)
        pred = predictor.predict_single_npy_array(t2_array[None, None, ...], properties)
        return pred[0].astype(np.uint8)

    Until the user provides a checkpoint path, this raises so the calling
    code reverts to the stub.
    """
    ckpt = cfg["segmentation"].get("nnunet_checkpoint")
    if not ckpt:
        raise RuntimeError(
            "segmentation.backend = 'nnunet' but no checkpoint path configured. "
            "Set segmentation.nnunet_checkpoint in preprocess.yaml."
        )
    raise NotImplementedError(
        "nnU-Net backend not yet wired up. See _nnunet_mask docstring for "
        "the expected interface. Download the Radboud checkpoint from "
        "https://github.com/DIAGNijmegen/AbdomenMRUS-prostate-segmentation"
    )
