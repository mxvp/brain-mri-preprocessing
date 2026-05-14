"""GDC REST API client + file downloader.

Public open-data endpoints, no auth needed for the projects we target
(TCGA-GBM, TCGA-LGG, CPTAC-3).

Two operations:
    search_rnaseq(project_id, submitter_ids)  → DataFrame of matching files
    download_files(manifest, dest)            → streams each file_id to disk
"""
from __future__ import annotations

import hashlib
import json
import logging
import time
from pathlib import Path

import pandas as pd
import requests

log = logging.getLogger(__name__)

GDC_FILES_URL = "https://api.gdc.cancer.gov/files"
GDC_DATA_URL = "https://api.gdc.cancer.gov/data"

DEFAULT_FIELDS = [
    "file_id", "file_name", "data_type", "data_category", "data_format",
    "experimental_strategy", "analysis.workflow_type", "access",
    "cases.submitter_id", "cases.case_id", "cases.disease_type",
    "cases.samples.sample_type", "cases.samples.submitter_id",
    "md5sum", "file_size",
]


# ---- search --------------------------------------------------------------

def search_rnaseq(
    project_id: str,
    submitter_ids: list[str],
    fields: list[str] | None = None,
    page_size: int = 500,
    sleep_between_pages: float = 0.2,
) -> pd.DataFrame:
    """All RNA-seq files for the given cases in a project, flat (one row per file×sample)."""
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
        if frm + page_size >= data["pagination"]["total"]:
            break
        frm += page_size
        time.sleep(sleep_between_pages)
    return _flatten(hits, project_id)


def _flatten(hits: list[dict], project_id: str) -> pd.DataFrame:
    rows = []
    for h in hits:
        for c in h.get("cases") or [{}]:
            for s in c.get("samples") or [{}]:
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
    """The canonical GDC RNA-seq gene-expression TSV (drop BAMs, fusions, splice junctions)."""
    return df[
        (df["workflow_type"] == "STAR - Counts") &
        (df["data_type"] == "Gene Expression Quantification")
    ].copy()


def write_gdc_manifest(df: pd.DataFrame, path) -> int:
    """Write a `gdc-client`-compatible TSV manifest."""
    m = df[["file_id", "file_name", "md5sum", "file_size"]].rename(columns={
        "file_id": "id", "file_name": "filename",
        "md5sum": "md5", "file_size": "size",
    })
    m["state"] = "validated"
    m.to_csv(path, sep="\t", index=False)
    return len(m)


# ---- download ------------------------------------------------------------

def _md5(path: Path, chunk: int = 1 << 16) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for buf in iter(lambda: f.read(chunk), b""):
            h.update(buf)
    return h.hexdigest()


def download_files(
    manifest: pd.DataFrame,
    dest: Path,
    verify_md5: bool = True,
    overwrite: bool = False,
    session: requests.Session | None = None,
) -> dict[str, int]:
    """Stream each manifest row into dest/{file_id}/{filename}. Idempotent."""
    dest.mkdir(parents=True, exist_ok=True)
    sess = session or requests.Session()
    stats = {"downloaded": 0, "skipped": 0, "failed": 0}

    for _, row in manifest.iterrows():
        file_id = row["file_id"]
        file_name = row["file_name"]
        target_dir = dest / file_id
        target = target_dir / file_name

        if target.exists() and not overwrite:
            if verify_md5 and pd.notna(row.get("md5sum")):
                if _md5(target) == row["md5sum"]:
                    stats["skipped"] += 1
                    continue
                log.warning(f"md5 mismatch on {target.name}, re-downloading")
            else:
                stats["skipped"] += 1
                continue

        target_dir.mkdir(parents=True, exist_ok=True)
        try:
            r = sess.get(f"{GDC_DATA_URL}/{file_id}", stream=True, timeout=120)
            r.raise_for_status()
            tmp = target.with_suffix(target.suffix + ".part")
            with open(tmp, "wb") as f:
                for chunk in r.iter_content(chunk_size=1 << 16):
                    f.write(chunk)
            if verify_md5 and pd.notna(row.get("md5sum")):
                got = _md5(tmp)
                if got != row["md5sum"]:
                    tmp.unlink()
                    raise RuntimeError(f"md5 mismatch: expected {row['md5sum']}, got {got}")
            tmp.replace(target)
            stats["downloaded"] += 1
        except Exception as e:
            log.error(f"failed {file_id} ({file_name}): {e}")
            stats["failed"] += 1
    return stats
