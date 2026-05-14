"""Load paired (imaging latents, expression) for downstream modeling.

Single entry point: `load_paired_dataset(latents_path)`. Joins an encoder's
latent file (.pt produced by GBM_MAE/scripts/encode_dataset.py) with the
expression matrices written by `python -m multimodal matrix`.

Defaults bake in these conventions (override via kwargs):
  - One row per imaging SUBJECT (not sample, not modality)
  - Modality: t1c (post-contrast T1), fallback t1
  - Sample type: Primary Tumor only
  - Multi-sample subjects: TPM averaged across primary samples
  - Genes: protein-coding only, expressed (TPM ≥ 1) in ≥ 10 subjects
  - Target: log1p(TPM)
  - QC: drops samples with protein_coding_fraction < 0.7
  - Split: subject-level, cohort-stratified, fixed seed

`parse_latent_path` is exposed separately so the verify CLI and tests can
exercise it without loading large .pt files.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import torch

from . import REPO_ROOT

log = logging.getLogger(__name__)

MATRICES_DIR = REPO_ROOT / "data" / "multimodal" / "matrices"


# ---- cohort registry -----------------------------------------------------
#
# Each entry is (project_id, regex over a normalized stem, canonicalizer).
# Add a new cohort = one new line. The regex runs on a path stem that has
# already had the BIDS `sub-` prefix and the nnUNet `_NNNN` channel suffix
# stripped, so patterns don't need to encode those.

@dataclass(frozen=True)
class CohortSpec:
    project_id: str
    pattern: re.Pattern
    canonicalize: Callable[[re.Match], str]


# `_` is a regex word char (`\w`), so `\b` doesn't fire next to underscores —
# which kills boundary detection in BIDS-style filenames like `..._t1_...`.
# Use explicit "not-followed-by-alphanumeric" assertions instead.
_NO_ALNUM_BEFORE = r"(?<![A-Z0-9])"
_NO_ALNUM_AFTER  = r"(?![A-Z0-9])"

REGISTRY: list[CohortSpec] = [
    # CPTAC: C3L00016 → C3L-00016 (also accept already-dashed)
    CohortSpec(
        "CPTAC-3",
        re.compile(_NO_ALNUM_BEFORE + r"C3([LN])-?(\d+)" + _NO_ALNUM_AFTER,
                   re.IGNORECASE),
        lambda m: f"C3{m.group(1).upper()}-{m.group(2)}",
    ),
    # TCGA: accept dashed (TCGA-06-1084) or dashless (TCGA020003, TCGAHTA61C).
    # TCGA submitter IDs are exactly 4 (prefix) + 2 (TSS) + 4 (participant).
    CohortSpec(
        "TCGA",
        re.compile(_NO_ALNUM_BEFORE + r"TCGA-?([A-Z0-9]{2})-?([A-Z0-9]{4})" + _NO_ALNUM_AFTER,
                   re.IGNORECASE),
        lambda m: f"TCGA-{m.group(1).upper()}-{m.group(2).upper()}",
    ),
    # UPENN: accept UPENN-GBM-00001 or UPENNGBM00001
    CohortSpec(
        "UPENN-GBM",
        re.compile(_NO_ALNUM_BEFORE + r"UPENN-?GBM-?(\d+)" + _NO_ALNUM_AFTER,
                   re.IGNORECASE),
        lambda m: f"UPENN-GBM-{m.group(1)}",
    ),
    # UCSF-PDGM
    CohortSpec(
        "UCSF-PDGM",
        re.compile(_NO_ALNUM_BEFORE + r"UCSF-?PDGM-?(\d+)" + _NO_ALNUM_AFTER,
                   re.IGNORECASE),
        lambda m: f"UCSF-PDGM-{m.group(1)}",
    ),
]

# Modality token in the filename (e.g. `..._t1c_preprocessed.nii.gz`).
# Intentionally does NOT infer modality from the BraTS `_0000` channel
# suffix — that mapping is encoder/dataset-specific and not safe to bake in.
MODALITY_RE = re.compile(r"_(t1c|t1gd|t1ce|t1|t2|flair)[_.]", re.IGNORECASE)

_SUB_PREFIX_RE = re.compile(r"^sub-", re.IGNORECASE)
_CHANNEL_SUFFIX_RE = re.compile(r"_\d{4}$")          # nnUNet channel index, e.g. _0000


def _strip_extensions(name: str) -> str:
    """foo.nii.gz → foo  (strip any number of suffixes)."""
    while True:
        new = Path(name).stem
        if new == name:
            return name
        name = new


def _normalize_stem(path: str) -> str:
    """Filename without extensions, BIDS `sub-` prefix, or nnUNet channel suffix."""
    stem = _strip_extensions(Path(path).name)
    stem = _SUB_PREFIX_RE.sub("", stem)
    stem = _CHANNEL_SUFFIX_RE.sub("", stem)
    return stem


def parse_latent_path(path: str) -> tuple[str | None, str | None, str | None]:
    """Extract (project_id, canonical submitter_id, modality) from a volume path.

    All three may be None. Modality detection is filename-token-based only
    (e.g. `..._t1c_...`); channel-suffix inference is intentionally NOT done.
    """
    mod_match = MODALITY_RE.search(path)
    modality = mod_match.group(1).lower() if mod_match else None
    if modality in ("t1gd", "t1ce"):
        modality = "t1c"

    stem = _normalize_stem(path)
    for spec in REGISTRY:
        # Try the stripped stem first (clean), then the full path (cohort name
        # in a parent directory still matches).
        m = spec.pattern.search(stem) or spec.pattern.search(path)
        if m:
            return spec.project_id, spec.canonicalize(m), modality
    return None, None, modality


# ---- dataset assembly ----------------------------------------------------

@dataclass
class PairedDataset:
    X: np.ndarray            # (N, D) imaging latents
    y: np.ndarray            # (N, G) log1p(TPM) expression
    subject_ids: np.ndarray  # (N,)   canonical submitter IDs
    cohort: np.ndarray       # (N,)   final project_id (TCGA-GBM / TCGA-LGG / CPTAC-3)
    split: np.ndarray        # (N,)   'train' / 'val' / 'test'
    gene_ids: np.ndarray     # (G,)   versioned Ensembl IDs
    gene_names: np.ndarray   # (G,)   HGNC symbols
    modality: str            # picked modality (mode across rows)


def load_paired_dataset(
    latents_path: str | Path,
    matrices_dir: str | Path | None = None,
    modality: str = "t1c",
    fallback_modality: str | None = "t1",
    sample_type: str = "Primary Tumor",
    gene_types: tuple[str, ...] = ("protein_coding",),
    min_subjects_expressed: int = 10,
    min_tpm: float = 1.0,
    pcf_min: float = 0.7,
    cohorts: tuple[str, ...] | None = None,
    split_seed: int = 0,
    val_frac: float = 0.15,
    test_frac: float = 0.15,
) -> PairedDataset:
    """Load + join + filter + split. See module docstring for conventions.

    Modality handling: prefers exact filename match on `modality`, then
    `fallback_modality`, then a path whose modality wasn't in the filename
    (trusts the caller's `modality` choice for single-modality latent files
    like `TCGA_T1C_V2/features.pt`).
    """
    matrices_dir = Path(matrices_dir) if matrices_dir else MATRICES_DIR

    # 1. Latents ------------------------------------------------------------
    blob = torch.load(latents_path, map_location="cpu", weights_only=False)
    if "features" in blob:
        latents = blob["features"]
    elif "latents" in blob:
        latents = blob["latents"]
    else:
        raise KeyError(
            f"{latents_path}: neither 'features' nor 'latents' key in blob. "
            f"Keys present: {list(blob.keys())}"
        )
    if latents.ndim == 3:           # (N, T, D) → mean-pool tokens
        latents = latents.mean(dim=1)
    latents = latents.numpy().astype(np.float32)
    paths = blob["file_paths"]

    parsed = pd.DataFrame(
        [parse_latent_path(p) for p in paths],
        columns=["project_hint", "submitter_id", "modality"],
    )
    parsed["latent_idx"] = np.arange(len(parsed))
    parsed = parsed.dropna(subset=["submitter_id"]).reset_index(drop=True)

    # 2. Expression tables --------------------------------------------------
    sample_meta = pd.read_parquet(matrices_dir / "sample_meta.parquet")
    gene_meta   = pd.read_parquet(matrices_dir / "gene_meta.parquet")
    tpm         = pd.read_parquet(matrices_dir / "tpm.parquet")

    sample_meta = sample_meta[
        (sample_meta["sample_type"] == sample_type) &
        (sample_meta["protein_coding_fraction"] >= pcf_min)
    ]
    if cohorts:
        sample_meta = sample_meta[sample_meta["project_id"].isin(cohorts)]
    tpm = tpm.loc[sample_meta.index]

    # 3. Subject-level expression (mean over primary samples) ---------------
    sm_keys = sample_meta[["submitter_id", "project_id"]].reset_index(drop=True)
    tpm_subj = (
        tpm.reset_index(drop=True)
           .assign(submitter_id=sm_keys["submitter_id"].values,
                   project_id=sm_keys["project_id"].values)
           .groupby(["project_id", "submitter_id"])
           .mean(numeric_only=True)
           .reset_index()
    )

    # 4. One latent per subject, preferring `modality` then fallback then None
    def _pick_one(sub: pd.DataFrame) -> pd.Series:
        mod_col = sub["modality"]
        for m in (modality, fallback_modality):
            if m is None:
                continue
            mask = mod_col.notna() & (mod_col.str.lower() == m.lower())
            hit = sub[mask]
            if len(hit):
                return hit.iloc[0]
        # Accept rows without explicit modality token (single-modality .pt files)
        null = sub[mod_col.isna()]
        if len(null):
            return null.iloc[0]
        return sub.iloc[0]

    picked_rows = [_pick_one(g) for _, g in parsed.groupby("submitter_id")]
    picked = pd.DataFrame(picked_rows).reset_index(drop=True)
    chosen_mod = (picked["modality"].dropna().mode().iloc[0]
                  if picked["modality"].notna().any() else modality)

    # 5. Join: expression subjects ⨝ latent subjects ------------------------
    joined = tpm_subj.merge(picked[["submitter_id", "latent_idx", "modality"]],
                            on="submitter_id", how="inner")
    if joined.empty:
        log.warning("No paired subjects after join. Encode TCGA / CPTAC first, "
                    "or check `python -m multimodal verify --latents <path>`.")
        return PairedDataset(
            X=np.zeros((0, latents.shape[1]), dtype=np.float32),
            y=np.zeros((0, 0), dtype=np.float32),
            subject_ids=np.array([]), cohort=np.array([]),
            split=np.array([]), gene_ids=np.array([]), gene_names=np.array([]),
            modality=modality,
        )

    # 6. Gene filter --------------------------------------------------------
    gene_cols = [c for c in joined.columns if c.startswith("ENSG")]
    tpm_arr = joined[gene_cols].values.astype(np.float32)
    type_keep = gene_meta.set_index("gene_id").loc[gene_cols, "gene_type"].isin(gene_types).values
    expr_keep = (tpm_arr >= min_tpm).sum(axis=0) >= min_subjects_expressed
    gene_keep = type_keep & expr_keep
    y = np.log1p(tpm_arr[:, gene_keep]).astype(np.float32)
    kept_gene_ids = np.array(gene_cols)[gene_keep]
    kept_gene_names = gene_meta.set_index("gene_id").loc[kept_gene_ids, "gene_name"].values

    # 7. Imaging features ---------------------------------------------------
    X = latents[joined["latent_idx"].values]

    # 8. Subject-level cohort-stratified split ------------------------------
    rng = np.random.default_rng(split_seed)
    split = np.empty(len(joined), dtype="<U5")
    for project, idx in joined.groupby("project_id").groups.items():
        idx = np.array(idx)
        rng.shuffle(idx)
        n = len(idx)
        n_test = max(1, int(round(n * test_frac)))
        n_val  = max(1, int(round(n * val_frac)))
        split[idx[:n_test]] = "test"
        split[idx[n_test:n_test + n_val]] = "val"
        split[idx[n_test + n_val:]] = "train"

    log.info(f"Paired dataset: N={len(joined)} subjects, G={int(gene_keep.sum())} genes, "
             f"D={X.shape[1]} latent-dim, modality='{chosen_mod}'")
    log.info(f"  per-cohort: {joined['project_id'].value_counts().to_dict()}")
    log.info(f"  split:      {pd.Series(split).value_counts().to_dict()}")

    return PairedDataset(
        X=X, y=y,
        subject_ids=joined["submitter_id"].values,
        cohort=joined["project_id"].values,
        split=split,
        gene_ids=kept_gene_ids,
        gene_names=kept_gene_names,
        modality=chosen_mod,
    )


# ---- diagnostics --------------------------------------------------------

def verify_latents(latents_path: str | Path) -> dict:
    """Run parse_latent_path on every entry in a .pt file. Returns a report.

    Report dict keys: total, matched, unmatched, per_cohort, per_modality,
                      unmatched_examples (first 5).
    """
    blob = torch.load(latents_path, map_location="cpu", weights_only=False)
    paths = blob.get("file_paths") or blob.get("paths")
    if paths is None:
        raise KeyError(f"{latents_path}: no 'file_paths' or 'paths' key.")
    parsed = [parse_latent_path(p) for p in paths]
    matched = [p for p, r in zip(paths, parsed) if r[0]]
    unmatched = [p for p, r in zip(paths, parsed) if not r[0]]
    from collections import Counter
    return {
        "total": len(paths),
        "matched": len(matched),
        "unmatched": len(unmatched),
        "per_cohort": dict(Counter(r[0] for r in parsed if r[0])),
        "per_modality": dict(Counter(r[2] for r in parsed if r[2])),
        "unmatched_examples": unmatched[:5],
    }
