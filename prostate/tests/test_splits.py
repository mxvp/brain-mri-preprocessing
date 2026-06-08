"""Tests for patient-grouped train/val split."""
import numpy as np
import pandas as pd
import pytest

from prostate.splits import patient_grouped_split, verify_no_leakage


def _fake_manifest():
    rows = []
    for source, n_patients in [("picai", 100), ("tcia_biopsy", 50), ("promise12", 20)]:
        for p in range(n_patients):
            for mod in ("t2", "adc", "dwi"):
                rows.append({"source": source,
                             "patient_id": f"p{p:04d}",
                             "modality": mod,
                             "view": "axial"})
    return pd.DataFrame(rows)


def test_no_patient_in_both_splits():
    df = _fake_manifest()
    split = patient_grouped_split(df, val_fraction=0.1, seed=42)
    verify_no_leakage(df, split)


def test_val_fraction_approximately_respected():
    df = _fake_manifest()
    split = patient_grouped_split(df, val_fraction=0.1, seed=42)
    for source in df["source"].unique():
        sub = df[df["source"] == source]
        sub_split = split.loc[sub.index]
        val_patients = sub.loc[sub_split == "val", "patient_id"].nunique()
        all_patients = sub["patient_id"].nunique()
        # within ±2 patients of the requested 10%
        target = int(round(all_patients * 0.1))
        assert abs(val_patients - target) <= 1


def test_split_is_deterministic_with_seed():
    df = _fake_manifest()
    s1 = patient_grouped_split(df, val_fraction=0.1, seed=42)
    s2 = patient_grouped_split(df, val_fraction=0.1, seed=42)
    assert (s1.values == s2.values).all()


def test_different_seeds_yield_different_splits():
    df = _fake_manifest()
    s1 = patient_grouped_split(df, val_fraction=0.1, seed=0)
    s2 = patient_grouped_split(df, val_fraction=0.1, seed=1)
    assert (s1.values != s2.values).any()


def test_all_modalities_of_a_patient_land_together():
    df = _fake_manifest()
    split = patient_grouped_split(df, val_fraction=0.1, seed=42)
    df2 = df.copy()
    df2["__split"] = split
    per_patient = df2.groupby(["source", "patient_id"])["__split"].nunique()
    assert (per_patient == 1).all()


def test_leakage_detection_catches_bad_split():
    df = _fake_manifest()
    bad_split = pd.Series(["train"] * len(df), index=df.index)
    # Force two rows of the same (source, patient_id) into different splits.
    bad_split.iloc[0] = "val"
    with pytest.raises(AssertionError, match="leakage"):
        verify_no_leakage(df, bad_split)
