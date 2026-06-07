"""Per-modality intensity normalization.

PI-CAI convention:
  T2w + DWI: z-score over nonzero voxels (i.e. over the masked prostate region
             once we've cropped, since outside the mask is zero).
  ADC:       clip to 0.5-99.5 percentile, then z-score. ADC has hard-to-control
             outliers from CSF/bone/artifacts that wreck plain z-scoring.
"""
from __future__ import annotations

import numpy as np


def zscore_nonzero(arr: np.ndarray) -> np.ndarray:
    """Z-score using mean/std computed over nonzero voxels only. Zeros stay zero."""
    arr = arr.astype(np.float32, copy=True)
    mask = arr != 0
    if mask.sum() == 0:
        return arr
    m, s = arr[mask].mean(), arr[mask].std()
    if s == 0:
        return arr - m
    out = np.zeros_like(arr, dtype=np.float32)
    out[mask] = (arr[mask] - m) / s
    return out


def clip_then_zscore(
    arr: np.ndarray,
    clip_low: float = 0.5,
    clip_high: float = 99.5,
) -> np.ndarray:
    """Clip to percentile range, then z-score over nonzero voxels."""
    arr = arr.astype(np.float32, copy=True)
    mask = arr != 0
    if mask.sum() == 0:
        return arr
    lo, hi = np.percentile(arr[mask], [clip_low, clip_high])
    arr[mask] = np.clip(arr[mask], lo, hi)
    return zscore_nonzero(arr)


def apply(arr: np.ndarray, spec: dict) -> np.ndarray:
    """Dispatch on `spec['method']`."""
    method = spec["method"]
    if method == "zscore_nonzero":
        return zscore_nonzero(arr)
    if method == "clip_then_zscore":
        return clip_then_zscore(arr, spec.get("clip_low", 0.5),
                                spec.get("clip_high", 99.5))
    raise ValueError(f"unknown normalize method: {method!r}")
