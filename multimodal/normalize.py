"""Gene filtering + normalization helpers for the expression matrix.

We keep the *stored* matrix as raw counts + TPM. Normalization (log1p, z-score,
TMM, etc.) is applied at load time — but the basic helpers live here so a
single source-of-truth keeps the modeling code consistent.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def filter_genes_by_type(genes: pd.DataFrame, gene_types: list[str]) -> pd.DataFrame:
    """Keep only rows whose `gene_type` is in `gene_types`."""
    return genes[genes["gene_type"].isin(gene_types)].copy()


def filter_genes_by_expression(
    tpm: pd.DataFrame,
    min_subjects_expressed: int = 10,
    min_tpm: float = 1.0,
) -> np.ndarray:
    """Return a boolean column mask: True for genes expressed in ≥ N subjects.

    Args:
        tpm: rows = subjects, columns = genes
    """
    expressed = (tpm >= min_tpm).sum(axis=0)
    return (expressed >= min_subjects_expressed).values


def log1p_tpm(tpm: pd.DataFrame) -> pd.DataFrame:
    """log1p-transform a TPM matrix. Preserves index/columns."""
    return np.log1p(tpm)


def library_size_normalize_counts(counts: pd.DataFrame, target_sum: float = 1e6) -> pd.DataFrame:
    """Convert raw counts to CPM-like (counts-per-million). Use for non-length-corrected scaling."""
    lib = counts.sum(axis=1).replace(0, np.nan)
    return counts.div(lib, axis=0) * target_sum
