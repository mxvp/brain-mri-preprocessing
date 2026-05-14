"""Imaging ↔ expression pairing.

Given an imaging subject inventory (from cohorts.load_imaging_subjects) and a
sample_meta table (from matrix.build_expression_matrix), emits a long-format
pairing table:

    imaging_id | submitter_id | sample_id | cohort | sample_type | file_id

This is the canonical join key for downstream modeling. One imaging subject
can map to multiple sample rows (primary + recurrent), which is fine — the
modeling code decides how to collapse if needed.
"""
from __future__ import annotations

import pandas as pd


def make_pairs(imaging: pd.DataFrame, sample_meta: pd.DataFrame) -> pd.DataFrame:
    """Inner-join imaging subjects with expression samples on submitter_id."""
    s = sample_meta.reset_index()
    df = imaging.merge(s, on="submitter_id", how="inner",
                       suffixes=("_imaging", "_expression"))
    cols = [
        "cohort", "imaging_id", "submitter_id",
        "sample_id", "sample_submitter", "sample_type",
        "file_id", "project_id",
    ]
    return df[[c for c in cols if c in df.columns]].drop_duplicates()


def coverage_summary(imaging: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    """One-row-per-cohort summary of pairing coverage."""
    rows = []
    for cohort, sub in imaging.groupby("cohort"):
        paired_subjects = pairs[pairs["cohort"] == cohort]["imaging_id"].nunique()
        rows.append({
            "cohort":              cohort,
            "imaging_subjects":    sub["imaging_id"].nunique(),
            "paired_subjects":     paired_subjects,
            "paired_samples":      int((pairs["cohort"] == cohort).sum()),
            "has_public":          bool(sub["has_public_expression"].iloc[0]),
        })
    return pd.DataFrame(rows).sort_values("paired_subjects", ascending=False)
