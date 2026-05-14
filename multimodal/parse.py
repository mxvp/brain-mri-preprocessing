"""Parse GDC STAR-Counts gene expression TSVs into clean per-sample frames.

GDC STAR-Counts files look like:
    # gene-model: GENCODE v36
    gene_id  gene_name  gene_type  unstranded  stranded_first  stranded_second  tpm_unstranded  fpkm_unstranded  fpkm_uq_unstranded
    N_unmapped  -  -  100  100  100  -  -  -
    N_multimapping  -  -  ...
    N_noFeature  -  -  ...
    N_ambiguous  -  -  ...
    ENSG00000000003.15  TSPAN6  protein_coding  3417  1716  1701  20.74  9.81  9.65
    ...

The first 4 non-gene rows (N_unmapped etc.) are STAR alignment QC and must be
skipped. The header starts at the line beginning with `gene_id`.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# Columns we keep from each TSV.
EXPRESSION_COLS = [
    "gene_id", "gene_name", "gene_type",
    "unstranded", "tpm_unstranded", "fpkm_unstranded",
]

# Rows that aren't genes — STAR alignment-level QC metrics.
NON_GENE_ROWS = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}


def parse_star_counts(tsv_path: Path) -> pd.DataFrame:
    """Read one STAR-Counts TSV and return the gene-level expression frame."""
    df = pd.read_csv(tsv_path, sep="\t", comment="#")
    df = df[~df["gene_id"].isin(NON_GENE_ROWS)].copy()
    df["gene_id_clean"] = df["gene_id"].str.split(".").str[0]
    keep = EXPRESSION_COLS + ["gene_id_clean"]
    return df[keep].reset_index(drop=True)


def extract_qc_metrics(tsv_path: Path) -> dict:
    """Pull the STAR alignment-level QC counts (N_unmapped etc.) plus a few derived stats."""
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
    """Walk a download directory and list every `*.rna_seq.augmented_star_gene_counts.tsv`.

    Returns DataFrame with [file_id, tsv_path] for cross-joining with manifests.
    Each TSV lives under {root}/{file_id}/{file_name}.
    """
    rows = []
    for tsv in root.glob("*/*.rna_seq.augmented_star_gene_counts.tsv"):
        rows.append({"file_id": tsv.parent.name, "tsv_path": str(tsv)})
    return pd.DataFrame(rows)
