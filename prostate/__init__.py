"""Prostate MRI preprocessing pipeline.

Per-patient: T2w + ADC + high-b DWI from dcm2niix output → cropped, registered,
normalized, prostate-mask-aware NIfTI ready for foundation-model encoding.

Follows the PI-CAI / Radboud nnU-Net community standard:
  N4 bias correction → prostate gland segmentation → rigid register DWI+ADC
  to T2w → crop around gland → resample to 0.5×0.5×3mm → per-modality
  normalize (T2/DWI z-score; ADC clip-then-z) → save 4-channel volume.

Config-driven via configs/preprocess.yaml. Idempotent steps; re-runs skip
existing outputs.
"""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIGS_DIR = PACKAGE_ROOT / "configs"
DATA_ROOT = REPO_ROOT / "data" / "prostate"

__all__ = ["REPO_ROOT", "PACKAGE_ROOT", "CONFIGS_DIR", "DATA_ROOT"]
