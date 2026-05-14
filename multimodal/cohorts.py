"""Cohort-aware extraction of canonical subject IDs from imaging metadata.

Each cohort has its own filename / metadata convention. This module exposes a
single function `load_imaging_subjects(cfg)` that takes the parsed cohorts.yaml
and returns a long-format DataFrame:

    cohort | imaging_id | submitter_id | project_id

`submitter_id` is the join key for the molecular repository (GDC submitter_id
for TCGA / CPTAC; the cohort's own ID for UPENN / UCSF).
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from . import REPO_ROOT, CONFIGS_DIR


def load_cohorts_config(path: Path | None = None) -> dict[str, Any]:
    path = path or (CONFIGS_DIR / "cohorts.yaml")
    with open(path) as f:
        return yaml.safe_load(f)


# ---------- per-cohort parsers ---------------------------------------------

def _parse_tcga(spec: dict, root: Path) -> pd.DataFrame:
    """TCGA: read tcga_master.csv, optionally filter by dx."""
    csv = root / spec["metadata_csv"]
    df = pd.read_csv(csv, low_memory=False)
    for k, v in (spec.get("filter") or {}).items():
        df = df[df[k] == v]
    return pd.DataFrame({
        "imaging_id":    df["subject_id"].astype(str),
        "submitter_id":  df["subject_id"].astype(str),
    })


def _parse_cptac(spec: dict, root: Path) -> pd.DataFrame:
    """CPTAC: scan a directory of `sub-C3L00016_0000.nii.gz` filenames.

    The GDC submitter is the dashed form `C3L-00016`.
    """
    d = Path(spec["dir"])
    if not d.exists():
        return pd.DataFrame(columns=["imaging_id", "submitter_id"])
    rows = []
    for f in sorted(d.glob(spec.get("glob", "*.nii.gz"))):
        m = re.match(r"sub-(C3[LN])(\d+)_\d+\.nii\.gz", f.name)
        if not m:
            continue
        rows.append({
            "imaging_id":   f.name.replace(".nii.gz", ""),
            "submitter_id": f"{m.group(1)}-{m.group(2)}",
        })
    return pd.DataFrame(rows)


def _parse_filelist(spec: dict, root: Path, pattern: re.Pattern, cohort_key: str) -> pd.DataFrame:
    """Generic filelist parser. Expects lines like `data/preprocessed/<cohort>/<filename>`."""
    f = root / spec["filelist"]
    if not f.exists():
        return pd.DataFrame(columns=["imaging_id", "submitter_id"])
    subjects: set[str] = set()
    with open(f) as fh:
        for line in fh:
            parts = line.strip().split("/")
            if len(parts) < 3 or parts[2].lower() != cohort_key:
                continue
            m = pattern.search(parts[-1])
            if m:
                subjects.add(m.group(1))
    df = pd.DataFrame({"submitter_id": sorted(subjects)})
    df["imaging_id"] = df["submitter_id"]
    return df


def _parse_upenn(spec: dict, root: Path) -> pd.DataFrame:
    return _parse_filelist(spec, root, re.compile(r"(UPENN-GBM-\d+)"), "upenn")


def _parse_ucsf(spec: dict, root: Path) -> pd.DataFrame:
    return _parse_filelist(spec, root, re.compile(r"(UCSF-PDGM-\d+)"), "ucsf")


PARSERS = {
    "tcga":  _parse_tcga,
    "cptac": _parse_cptac,
    "upenn": _parse_upenn,
    "ucsf":  _parse_ucsf,
}


# ---------- public API ------------------------------------------------------

def load_imaging_subjects(
    cfg: dict | None = None,
    cohorts: list[str] | None = None,
    repo_root: Path | None = None,
) -> pd.DataFrame:
    """Build a long DataFrame of imaging subjects across cohorts.

    Args:
        cfg: parsed cohorts.yaml (loaded fresh if None)
        cohorts: subset of cohort keys to include (default = all)
        repo_root: brain-mri-preprocessing root (default = inferred)

    Returns:
        DataFrame with columns [cohort, imaging_id, submitter_id, project_id].
    """
    cfg = cfg or load_cohorts_config()
    repo_root = repo_root or REPO_ROOT
    rows: list[pd.DataFrame] = []
    for key, spec in cfg["cohorts"].items():
        if cohorts and key not in cohorts:
            continue
        parser_key = spec["imaging"]["parser"]
        df = PARSERS[parser_key](spec["imaging"], repo_root)
        if df.empty:
            continue
        df["cohort"] = key
        df["project_id"] = spec["expression"].get("project_id")
        df["has_public_expression"] = spec["expression"].get("has_public", False)
        rows.append(df)
    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    cols = ["cohort", "imaging_id", "submitter_id", "project_id", "has_public_expression"]
    return out[[c for c in cols if c in out.columns]]
