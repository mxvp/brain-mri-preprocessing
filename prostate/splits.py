"""Patient-grouped, source-stratified train/val split.

A patient's volumes (T2 + ADC + DWI + multiple views if any) always land in
the same split — never some in train and some in val. Avoids the most common
form of leakage in medical-imaging datasets.

The split is per-source so each source contributes the same val fraction.
A fixed seed makes runs reproducible.
"""
from __future__ import annotations

import logging

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


def patient_grouped_split(
    manifest: pd.DataFrame,
    val_fraction: float = 0.05,
    seed: int = 0,
) -> pd.Series:
    """Return a Series of 'train'/'val' aligned with manifest.index.

    Splits are made per (source, patient_id) tuple, stratified by source so
    val composition mirrors the cohort.
    """
    rng = np.random.default_rng(seed)
    labels = pd.Series(index=manifest.index, dtype="<U5")

    for source, group in manifest.groupby("source"):
        patients = sorted(group["patient_id"].unique().tolist())
        rng.shuffle(patients)
        n_val = max(1, int(round(len(patients) * val_fraction)))
        val_patients = set(patients[:n_val])
        in_val = group["patient_id"].isin(val_patients)
        labels.loc[group.index[in_val]] = "val"
        labels.loc[group.index[~in_val]] = "train"
        log.info(f"  {source:<12s} {len(patients):>5d} patients → "
                 f"train: {len(patients) - n_val} patients / {int((~in_val).sum())} files  "
                 f"val: {n_val} patients / {int(in_val.sum())} files")

    return labels


def verify_no_leakage(manifest: pd.DataFrame, split: pd.Series) -> None:
    """Sanity check: no (source, patient_id) appears in both train and val."""
    df = manifest.copy()
    df["__split"] = split
    per_patient = df.groupby(["source", "patient_id"])["__split"].nunique()
    leaked = per_patient[per_patient > 1]
    if len(leaked):
        raise AssertionError(
            f"split leakage: {len(leaked)} (source, patient_id) tuples are in "
            f"both train and val. First few:\n{leaked.head()}"
        )
