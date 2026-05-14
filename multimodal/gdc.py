"""GDC REST API client (read-only, public endpoints, no auth required).

We use only two endpoints:
    POST /files   — search for files by filters
    POST /data    — bulk download (we delegate to gdc-client instead)
"""
from __future__ import annotations

import json
import time
from typing import Any

import pandas as pd
import requests

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"

# Fields we always want returned for each RNA-seq file hit.
DEFAULT_FIELDS = [
    "file_id", "file_name", "data_type", "data_category", "data_format",
    "experimental_strategy", "analysis.workflow_type", "access",
    "cases.submitter_id", "cases.case_id", "cases.disease_type",
    "cases.samples.sample_type", "cases.samples.submitter_id",
    "md5sum", "file_size",
]


def search_rnaseq(
    project_id: str,
    submitter_ids: list[str],
    fields: list[str] | None = None,
    page_size: int = 500,
    sleep_between_pages: float = 0.2,
) -> pd.DataFrame:
    """Return all RNA-seq files for the given cases in a project.

    The returned frame is flattened so each (file × case × sample) combo is one
    row. For STAR-Counts gene expression each file corresponds to one sample.
    """
    fields = fields or DEFAULT_FIELDS
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "experimental_strategy",     "value": ["RNA-Seq"]}},
            {"op": "in", "content": {"field": "cases.submitter_id",         "value": submitter_ids}},
        ],
    }
    hits: list[dict] = []
    frm = 0
    while True:
        payload = {
            "filters": json.dumps(filters),
            "fields": ",".join(fields),
            "format": "JSON",
            "size": page_size,
            "from": frm,
        }
        r = requests.post(GDC_FILES_URL, json=payload, timeout=60)
        r.raise_for_status()
        data = r.json()["data"]
        hits.extend(data["hits"])
        pagination = data["pagination"]
        if frm + page_size >= pagination["total"]:
            break
        frm += page_size
        time.sleep(sleep_between_pages)
    return _flatten(hits, project_id)


def _flatten(hits: list[dict], project_id: str) -> pd.DataFrame:
    rows = []
    for h in hits:
        cases = h.get("cases") or [{}]
        for c in cases:
            samples = c.get("samples") or [{}]
            for s in samples:
                rows.append({
                    "project_id":       project_id,
                    "submitter_id":     c.get("submitter_id"),
                    "case_id":          c.get("case_id"),
                    "sample_submitter": s.get("submitter_id"),
                    "sample_type":      s.get("sample_type"),
                    "file_id":          h["file_id"],
                    "file_name":        h["file_name"],
                    "data_type":        h.get("data_type"),
                    "data_category":    h.get("data_category"),
                    "workflow_type":    (h.get("analysis") or {}).get("workflow_type"),
                    "access":           h.get("access"),
                    "md5sum":           h.get("md5sum"),
                    "file_size":        h.get("file_size"),
                })
    return pd.DataFrame(rows)


def filter_star_counts_expression(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only the canonical STAR-Counts gene-expression TSVs.

    These are the files with rows = genes and columns = counts/TPM/FPKM.
    Drops fusion calls, splice junctions, and BAMs.
    """
    return df[
        (df["workflow_type"] == "STAR - Counts") &
        (df["data_type"] == "Gene Expression Quantification")
    ].copy()


def write_gdc_manifest(df: pd.DataFrame, path) -> int:
    """Write a `gdc-client`-compatible TSV manifest. Returns # rows written."""
    manifest = df[["file_id", "file_name", "md5sum", "file_size"]].rename(columns={
        "file_id": "id", "file_name": "filename",
        "md5sum": "md5", "file_size": "size",
    })
    manifest["state"] = "validated"
    manifest.to_csv(path, sep="\t", index=False)
    return len(manifest)
