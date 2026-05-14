"""Download GDC files via the public REST API (no auth needed for open data).

We avoid the `gdc-client` binary dependency. Each file is fetched by file_id
from the `/data` endpoint, streamed to disk, MD5-verified, and skipped on
re-runs (idempotent).
"""
from __future__ import annotations

import hashlib
import logging
from pathlib import Path

import pandas as pd
import requests

log = logging.getLogger(__name__)

GDC_DATA_URL = "https://api.gdc.cancer.gov/data"


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
    """Download every file in the manifest into `dest/{file_id}/{filename}`.

    Manifest columns required: file_id, file_name, md5sum.
    Returns counts of {downloaded, skipped, failed}.
    """
    dest.mkdir(parents=True, exist_ok=True)
    sess = session or requests.Session()
    stats = {"downloaded": 0, "skipped": 0, "failed": 0}

    for _, row in manifest.iterrows():
        file_id = row["file_id"]
        file_name = row["file_name"]
        target_dir = dest / file_id
        target = target_dir / file_name

        if target.exists() and not overwrite:
            if verify_md5 and "md5sum" in row and pd.notna(row["md5sum"]):
                if _md5(target) == row["md5sum"]:
                    stats["skipped"] += 1
                    continue
                log.warning(f"MD5 mismatch on existing {target.name}, redownloading")
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
            log.error(f"Failed {file_id} ({file_name}): {e}")
            stats["failed"] += 1
    return stats
