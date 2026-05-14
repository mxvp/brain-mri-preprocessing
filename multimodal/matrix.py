"""Build the subject × gene expression matrix from per-sample TSVs.

The output is a single parquet (or HDF5) file with three logical tables:
    counts       sample_index × gene_id   (int)
    tpm          sample_index × gene_id   (float32)
    gene_meta    gene_id, gene_name, gene_type
    sample_meta  sample_index → cohort, submitter_id, sample_submitter,
                 sample_type, file_id, qc fields

The sample (not subject) is the row identifier — a subject can contribute
multiple samples (primary + recurrent). The pairing module decides how to
collapse to one-vector-per-imaging-subject downstream.
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from .parse import discover_tsvs, parse_star_counts, extract_qc_metrics

log = logging.getLogger(__name__)


def build_expression_matrix(
    manifests: dict[str, pd.DataFrame],
    raw_root: Path,
    sample_filter: str = "all_tumor",
) -> dict[str, pd.DataFrame]:
    """Build counts/tpm/gene_meta/sample_meta tables from the manifests + TSVs.

    Args:
        manifests: mapping of project_id → manifest DataFrame (from gdc.search_rnaseq)
        raw_root:  root of downloaded TSVs (matrices under {raw_root}/<project>/<file_id>/...)
        sample_filter: 'primary_only' | 'all_tumor' | 'all'

    Returns:
        dict with keys 'counts', 'tpm', 'gene_meta', 'sample_meta'.
    """
    if sample_filter == "primary_only":
        keep_types = {"Primary Tumor"}
    elif sample_filter == "all_tumor":
        keep_types = {"Primary Tumor", "Recurrent Tumor"}
    elif sample_filter == "all":
        keep_types = None
    else:
        raise ValueError(f"unknown sample_filter: {sample_filter}")

    sample_meta_rows = []
    per_sample_counts: list[pd.Series] = []
    per_sample_tpm:    list[pd.Series] = []
    gene_meta_ref: pd.DataFrame | None = None

    for project, manifest in manifests.items():
        # Only keep STAR-Counts gene expression rows.
        m = manifest[
            (manifest["workflow_type"] == "STAR - Counts") &
            (manifest["data_type"] == "Gene Expression Quantification")
        ].copy()
        if keep_types is not None:
            m = m[m["sample_type"].isin(keep_types)]
        if m.empty:
            log.info(f"{project}: no rows after sample filter")
            continue

        discovered = discover_tsvs(raw_root / project)
        m = m.merge(discovered, on="file_id", how="inner")
        log.info(f"{project}: {len(m)} samples to ingest after disk join")

        for _, row in m.iterrows():
            tsv_path = Path(row["tsv_path"])
            df = parse_star_counts(tsv_path)
            qc = extract_qc_metrics(tsv_path)

            # Use the dotless Ensembl ID as gene key for stability.
            df = df.set_index("gene_id_clean")
            if gene_meta_ref is None:
                gene_meta_ref = df[["gene_id", "gene_name", "gene_type"]].copy()

            sample_id = row["sample_submitter"] or row["file_id"]
            per_sample_counts.append(df["unstranded"].rename(sample_id))
            per_sample_tpm.append(df["tpm_unstranded"].rename(sample_id))

            sample_meta_rows.append({
                "sample_id":        sample_id,
                "project_id":       row["project_id"],
                "submitter_id":     row["submitter_id"],
                "sample_submitter": row["sample_submitter"],
                "sample_type":      row["sample_type"],
                "file_id":          row["file_id"],
                **qc,
            })

    if not per_sample_counts:
        return {"counts": pd.DataFrame(), "tpm": pd.DataFrame(),
                "gene_meta": pd.DataFrame(), "sample_meta": pd.DataFrame()}

    counts = pd.concat(per_sample_counts, axis=1).T
    tpm    = pd.concat(per_sample_tpm, axis=1).T
    counts.index.name = "sample_id"
    tpm.index.name = "sample_id"
    sample_meta = pd.DataFrame(sample_meta_rows).set_index("sample_id")

    return {
        "counts":      counts.astype(np.int32),
        "tpm":         tpm.astype(np.float32),
        "gene_meta":   gene_meta_ref.reset_index().rename(columns={"gene_id_clean": "gene_id_short"}),
        "sample_meta": sample_meta,
    }


def write_parquet(tables: dict[str, pd.DataFrame], out_dir: Path) -> None:
    """Write each table to its own parquet under `out_dir/`."""
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        path = out_dir / f"{name}.parquet"
        df.to_parquet(path)
        log.info(f"  wrote {path}  shape={df.shape}")
