"""Build the curated expression product: counts/tpm matrices, sample/gene metadata, pairs.

Encapsulates everything from raw STAR-Counts TSV → subject×gene matrix:
  - TSV parsing + per-sample QC extraction
  - Aggregation into counts/tpm/gene_meta/sample_meta tables
  - log1p, library-size norm, biotype/expression filters (helpers)
  - Imaging↔expression pairing + coverage summary
"""
from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


# ---- TSV parsing ---------------------------------------------------------

EXPRESSION_COLS = [
    "gene_id", "gene_name", "gene_type",
    "unstranded", "tpm_unstranded", "fpkm_unstranded",
]
NON_GENE_ROWS = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}


def parse_star_counts(tsv_path: Path) -> pd.DataFrame:
    """Read one STAR-Counts TSV and return the gene-level expression frame."""
    df = pd.read_csv(tsv_path, sep="\t", comment="#")
    df = df[~df["gene_id"].isin(NON_GENE_ROWS)].copy()
    df["gene_id_clean"] = df["gene_id"].str.split(".").str[0]
    return df[EXPRESSION_COLS + ["gene_id_clean"]].reset_index(drop=True)


def extract_qc_metrics(tsv_path: Path) -> dict:
    """STAR alignment-level QC + a few derived stats per sample."""
    df = pd.read_csv(tsv_path, sep="\t", comment="#")
    qc_rows = df[df["gene_id"].isin(NON_GENE_ROWS)].set_index("gene_id")
    qc = {f"{r.lower()}_reads": float(qc_rows.loc[r, "unstranded"])
          for r in NON_GENE_ROWS if r in qc_rows.index}

    genes = df[~df["gene_id"].isin(NON_GENE_ROWS)]
    qc["mapped_reads"] = float(genes["unstranded"].sum())
    qc["mapped_to_protein_coding"] = float(
        genes.loc[genes["gene_type"] == "protein_coding", "unstranded"].sum()
    )
    qc["protein_coding_fraction"] = (
        qc["mapped_to_protein_coding"] / qc["mapped_reads"]
        if qc["mapped_reads"] > 0 else 0.0
    )
    qc["n_genes_expressed_tpm1"] = int((genes["tpm_unstranded"] >= 1.0).sum())
    return qc


def discover_tsvs(root: Path) -> pd.DataFrame:
    """Walk a download directory and list every STAR-Counts TSV. Returns [file_id, tsv_path]."""
    rows = [{"file_id": tsv.parent.name, "tsv_path": str(tsv)}
            for tsv in root.glob("*/*.rna_seq.augmented_star_gene_counts.tsv")]
    return pd.DataFrame(rows)


# ---- matrix build --------------------------------------------------------

def build_expression_matrix(
    manifests: dict[str, pd.DataFrame],
    raw_root: Path,
    sample_filter: str = "all_tumor",
) -> dict[str, pd.DataFrame]:
    """Build counts/tpm/gene_meta/sample_meta from manifests + downloaded TSVs."""
    if sample_filter == "primary_only":
        keep = {"Primary Tumor"}
    elif sample_filter == "all_tumor":
        keep = {"Primary Tumor", "Recurrent Tumor"}
    elif sample_filter == "all":
        keep = None
    else:
        raise ValueError(f"unknown sample_filter: {sample_filter}")

    sample_rows, count_cols, tpm_cols = [], [], []
    gene_meta_ref: pd.DataFrame | None = None

    for project, manifest in manifests.items():
        m = manifest[
            (manifest["workflow_type"] == "STAR - Counts") &
            (manifest["data_type"] == "Gene Expression Quantification")
        ].copy()
        if keep is not None:
            m = m[m["sample_type"].isin(keep)]
        if m.empty:
            log.info(f"{project}: no rows after sample filter")
            continue

        m = m.merge(discover_tsvs(raw_root / project), on="file_id", how="inner")
        log.info(f"{project}: {len(m)} samples to ingest")

        for _, row in m.iterrows():
            tsv = Path(row["tsv_path"])
            genes = parse_star_counts(tsv).set_index("gene_id_clean")
            if gene_meta_ref is None:
                gene_meta_ref = genes[["gene_id", "gene_name", "gene_type"]].copy()

            sample_id = row["sample_submitter"] or row["file_id"]
            count_cols.append(genes["unstranded"].rename(sample_id))
            tpm_cols.append(genes["tpm_unstranded"].rename(sample_id))
            sample_rows.append({
                "sample_id":        sample_id,
                "project_id":       row["project_id"],
                "submitter_id":     row["submitter_id"],
                "sample_submitter": row["sample_submitter"],
                "sample_type":      row["sample_type"],
                "file_id":          row["file_id"],
                **extract_qc_metrics(tsv),
            })

    if not count_cols:
        empty = pd.DataFrame()
        return {"counts": empty, "tpm": empty, "gene_meta": empty, "sample_meta": empty}

    counts = pd.concat(count_cols, axis=1).T
    tpm = pd.concat(tpm_cols, axis=1).T
    counts.index.name = tpm.index.name = "sample_id"
    sample_meta = pd.DataFrame(sample_rows).set_index("sample_id")

    return {
        "counts":      counts.astype(np.int32),
        "tpm":         tpm.astype(np.float32),
        "gene_meta":   gene_meta_ref.reset_index().rename(columns={"gene_id_clean": "gene_id_short"}),
        "sample_meta": sample_meta,
    }


def write_parquet(tables: dict[str, pd.DataFrame], out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in tables.items():
        if df.empty:
            continue
        path = out_dir / f"{name}.parquet"
        df.to_parquet(path)
        log.info(f"  wrote {path}  shape={df.shape}")


# ---- normalization helpers (applied at load time, not baked in) ----------

def filter_genes_by_type(genes: pd.DataFrame, gene_types: list[str]) -> pd.DataFrame:
    return genes[genes["gene_type"].isin(gene_types)].copy()


def filter_genes_by_expression(tpm: pd.DataFrame, min_subjects: int = 10, min_tpm: float = 1.0) -> np.ndarray:
    """Bool mask over genes: True if expressed (tpm ≥ min_tpm) in ≥ N subjects."""
    return ((tpm >= min_tpm).sum(axis=0) >= min_subjects).values


def log1p_tpm(tpm: pd.DataFrame) -> pd.DataFrame:
    return np.log1p(tpm)


def library_size_normalize(counts: pd.DataFrame, target_sum: float = 1e6) -> pd.DataFrame:
    """Raw counts → CPM (no length correction)."""
    lib = counts.sum(axis=1).replace(0, np.nan)
    return counts.div(lib, axis=0) * target_sum


# ---- pairing -------------------------------------------------------------

def make_pairs(imaging: pd.DataFrame, sample_meta: pd.DataFrame) -> pd.DataFrame:
    """Inner-join imaging subjects with expression samples on submitter_id."""
    s = sample_meta.reset_index()
    df = imaging.merge(s, on="submitter_id", how="inner",
                       suffixes=("_imaging", "_expression"))
    cols = ["cohort", "imaging_id", "submitter_id", "sample_id", "sample_submitter",
            "sample_type", "file_id", "project_id"]
    return df[[c for c in cols if c in df.columns]].drop_duplicates()


def coverage_summary(imaging: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cohort, sub in imaging.groupby("cohort"):
        paired = pairs[pairs["cohort"] == cohort]
        rows.append({
            "cohort":             cohort,
            "imaging_subjects":   sub["imaging_id"].nunique(),
            "paired_subjects":    paired["imaging_id"].nunique(),
            "paired_samples":     len(paired),
            "has_public":         bool(sub["has_public_expression"].iloc[0]),
        })
    return pd.DataFrame(rows).sort_values("paired_subjects", ascending=False)
