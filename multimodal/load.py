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
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from . import REPO_ROOT

log = logging.getLogger(__name__)

MATRICES_DIR = REPO_ROOT / "data" / "multimodal" / "matrices"

# Patterns for extracting (project_id, submitter_id) from a preprocessed-volume path.
MODALITY_RE = re.compile(r"_(t1c|t1gd|t1ce|t1|t2|flair)[_.]", re.I)
# CPTAC files come as `sub-C3L00016_0000.nii.gz` with no explicit modality
# token; the user confirmed these are T1-post (T1c). Same for the TCGA dump
# the colleague provided. Treat any `_0000` suffix on a tumor cohort as T1c.
BRATS_CHANNEL_TO_MODALITY = {"_0000": "t1c"}


def parse_latent_path(path: str) -> tuple[str | None, str | None, str | None]:
    """Return (project_hint, submitter_id, modality) from a volume file path."""
    m = MODALITY_RE.search(path)
    modality = m.group(1).lower() if m else None
    if modality in ("t1gd", "t1ce"):
        modality = "t1c"
    # BraTS-style _0000 channel as modality fallback (CPTAC + TCGA colleague dump)
    if modality is None:
        for k, v in BRATS_CHANNEL_TO_MODALITY.items():
            if k in path:
                modality = v
                break

    # TCGA: project_hint is generic "TCGA" — final TCGA-GBM vs TCGA-LGG comes
    # from sample_meta after the join.
    m = re.search(r"(TCGA-[A-Z0-9]{2}-[A-Z0-9]{4})", path)
    if m:
        return "TCGA", m.group(1), modality
    # CPTAC: sub-C3L00016 → C3L-00016
    m = re.search(r"C3([LN])-?(\d+)", path)
    if m:
        return "CPTAC-3", f"C3{m.group(1)}-{m.group(2)}", modality
    # UPENN: accept both UPENN-GBM-00001 and UPENNGBM00001
    m = re.search(r"UPENN-?GBM-?(\d+)", path)
    if m:
        return "UPENN-GBM", f"UPENN-GBM-{m.group(1)}", modality
    # UCSF: accept UCSF-PDGM-NNNN
    m = re.search(r"UCSF-?PDGM-?(\d+)", path)
    if m:
        return "UCSF-PDGM", f"UCSF-PDGM-{m.group(1)}", modality
    return None, None, modality


@dataclass
class PairedDataset:
    X: np.ndarray            # (N, D) imaging latents
    y: np.ndarray            # (N, G) log1p(TPM) expression
    subject_ids: np.ndarray  # (N,)   submitter_id strings
    cohort: np.ndarray       # (N,)   project_id strings
    split: np.ndarray        # (N,)   'train' / 'val' / 'test'
    gene_ids: np.ndarray     # (G,)   versioned Ensembl IDs
    gene_names: np.ndarray   # (G,)
    modality: str            # which modality each X row corresponds to


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
    """Load + join + filter + split. See module docstring for conventions."""
    matrices_dir = Path(matrices_dir) if matrices_dir else MATRICES_DIR

    # 1. Latents ------------------------------------------------------------
    blob = torch.load(latents_path, map_location="cpu", weights_only=False)
    latents = blob["latents"]
    if latents.ndim == 3:           # (N, T, D) → mean-pool tokens
        latents = latents.mean(dim=1)
    latents = latents.numpy().astype(np.float32)
    paths = blob["file_paths"]

    parsed = pd.DataFrame(
        [parse_latent_path(p) for p in paths],
        columns=["project_hint", "submitter_id", "modality"],
    )
    parsed["latent_idx"] = np.arange(len(parsed))

    # 2. Expression tables --------------------------------------------------
    sample_meta = pd.read_parquet(matrices_dir / "sample_meta.parquet")
    gene_meta   = pd.read_parquet(matrices_dir / "gene_meta.parquet")
    tpm         = pd.read_parquet(matrices_dir / "tpm.parquet")

    # Sample QC + sample_type filter
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

    # 4. Pick one latent per subject (preferred modality, fallback) ---------
    def _pick_one(sub_df: pd.DataFrame) -> pd.Series:
        for mod in (modality, fallback_modality):
            if mod is None:
                continue
            hit = sub_df[sub_df["modality"] == mod]
            if len(hit):
                return hit.iloc[0]
        return sub_df.iloc[0]  # last resort: any modality

    parsed = parsed.dropna(subset=["submitter_id"]).reset_index(drop=True)
    picked_rows = []
    for sid, sub in parsed.groupby("submitter_id"):
        picked_rows.append(_pick_one(sub))
    picked = pd.DataFrame(picked_rows).reset_index(drop=True)
    chosen_mod = picked["modality"].mode().iloc[0] if len(picked) else modality

    # 5. Join: expression subjects ⨝ latent subjects ------------------------
    joined = tpm_subj.merge(picked[["submitter_id", "latent_idx", "modality"]],
                            on="submitter_id", how="inner")
    if joined.empty:
        log.warning("No paired subjects after join — latent file probably lacks "
                    "TCGA / CPTAC paths. Encode those cohorts first.")
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
        X=X,
        y=y,
        subject_ids=joined["submitter_id"].values,
        cohort=joined["project_id"].values,
        split=split,
        gene_ids=kept_gene_ids,
        gene_names=kept_gene_names,
        modality=chosen_mod,
    )
