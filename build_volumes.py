"""Build one big volumes.csv — one row per preprocessed .nii.gz file on disk,
joined with clinical metadata from each dataset's master CSV.

Usage:

    # On Sherlock:
    find <sanitized>/projects/brain-mri-preprocessing/data/preprocessed \
         <sanitized>/projects/GBM_MAE/data/UPENN/normalize \
         /oak/stanford/groups/ogevaert/data/brain_mri_tumor_project/UCSF-PDGM-v3 \
         -name "*.nii.gz" > files.txt

    # Locally:
    python build_volumes.py \
        --files files.txt \
        --metadata-root data/metadata \
        --output data/metadata/volumes.csv

Output schema (one row per file):
    file_path          — absolute Sherlock path
    dataset            — canonical dataset name (matches master)
    subject_id         — canonical subject ID
    session_id         — e.g. BAS1 / MR1 / ses-1 / NA
    modality           — atomic: t1 / t2 / t1c / t1gd / flair
    + all clinical cols from the dataset's master CSV (age, sex, dx, mmse, ...)
"""

import argparse
import logging
import re
from pathlib import Path

import pandas as pd


log = logging.getLogger(__name__)


# =============================================================================
# Parent dir name → which master CSV to load and which parser to apply
# =============================================================================

DATASET_MAP = {
    # Preprocessed subdir name → (master CSV relative path, parser function name)
    # NOTE: adni_old deliberately excluded — it's superseded by adni_full and
    # the files would duplicate subjects with different naming.
    "adni_full": ("ADNI/adni_master.csv", "parse_adni"),
    "ppmi":      ("PPMI/ppmi_master.csv", "parse_ppmi"),
    "oasis1":    ("OASIS1/oasis1_master.csv", "parse_oasis1"),
    "oasis2":    ("OASIS2/oasis2_master.csv", "parse_oasis2"),
    "ixi":       ("IXI/ixi_master.csv", "parse_ixi"),
    "schizo":    ("SCHIZO/schizo_master.csv", "parse_schizo"),
    "stanford":  (None, "parse_stanford"),  # Stanford has no master (no metadata pulled)
    "hcp_ya":    ("HCP_YA/hcp_ya_master.csv", "parse_hcp_ya"),
    "bgsp":      ("BGSP/bgsp_master.csv", "parse_bgsp"),
    "abide":     ("ABIDE/abide_master.csv", "parse_abide"),
    "abide2":    ("ABIDE2/abide2_master.csv", "parse_abide2"),
    "adhd200":   ("ADHD200/adhd200_master.csv", "parse_adhd200"),
    "corr":      ("CORR/corr_master.csv", "parse_corr"),
    "fcon1000":  ("FCON1000/fcon1000_master.csv", "parse_fcon1000"),
    "hbn":       ("HBN/hbn_master.csv", "parse_hbn"),
    "nki":       ("NKI/nki_master.csv", "parse_nki"),
    "nki2":      ("NKI2/nki2_master.csv", "parse_nki2"),
    "upenn":     ("UPENN/upenn_master.csv", "parse_upenn"),
    "ucsf":      ("UCSF/ucsf_master.csv", "parse_ucsf"),
}


# Canonical dataset names (applied in output): rename parent-dir keys to
# stable lowercase identifiers the downstream pipeline expects.
DATASET_CANONICAL_NAME = {
    "adni_full":    "adni",
    "hcp_ya":       "hcp_ya",
}


# Collapse aliases for the harmonized `dx` column.
# CN (ADNI/OASIS cognitively normal) == HC (healthy control) for conditioning.
DX_ALIASES = {
    "CN": "HC",
}

# Collapse aliases for the `modality` column. t1gd (Stanford) == t1c (UCSF):
# both are post-gadolinium contrast-enhanced T1.
MODALITY_ALIASES = {
    "t1gd": "t1c",
}


def _normalize_apoe(val) -> str | float:
    """Normalize APOE genotype to standard E{n}/E{n} format.
    ADNI stores as '34.0' (concatenated allele digits), PPMI as 'E3/E4'."""
    if pd.isna(val):
        return pd.NA
    s = str(val).strip()
    # Already in E_/E_ format
    if s.startswith("E"):
        return s
    # Numeric format: '34.0' or '34' → 'E3/E4'
    try:
        digits = s.replace(".", "").replace("0", "")  # '34.0' → '34'
        if len(digits) == 2 and digits.isdigit():
            return f"E{digits[0]}/E{digits[1]}"
    except (ValueError, IndexError):
        pass
    return pd.NA


# =============================================================================
# Per-dataset file parsers
# Each returns {subject_key: <col>, value: <...>, session_id: ..., modality: ...}
# where subject_key tells us which master column to join on.
# Returns None if the filename doesn't match (e.g. log files, wrong format).
# =============================================================================


# Common modality extraction from filename token
MODALITY_TOKEN = re.compile(r"_(t1gd|t1c|t1|t2|flair)_preprocessed\.nii\.gz$")


def _modality(name: str) -> str | None:
    m = MODALITY_TOKEN.search(name)
    return m.group(1) if m else None


def parse_adni(filename: str) -> dict | None:
    """ADNI_{subject}_{YYYYMMDD}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"ADNI_(.+?)_(\d{8})_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    subject, yyyymmdd, mod = m.groups()
    return {
        "lookup_cols": {"subject_id": subject,
                         "scan_date": f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:]}"},
        "session_id": pd.NA,
        "modality": mod,
    }


def parse_adni_old(filename: str) -> dict | None:
    """ADNI_{subject}_{mod}_preprocessed.nii.gz (old pipeline, no date)"""
    m = re.match(r"ADNI_(.+?)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    subject, mod = m.groups()
    # adni_old: just subject, no date — pick the earliest-dated row in master
    return {
        "lookup_cols": {"subject_id": subject},
        "session_id": "baseline",
        "modality": mod,
    }


def parse_ppmi(filename: str) -> dict | None:
    """PPMI_{PATNO}_{YYYYMMDD}_{mod}_preprocessed.nii.gz

    PPMI master is per-subject (baseline clinical), so we just match by PATNO.
    The filename's scan date is preserved in build_volumes output via `scan_date`
    column override so per-scan longitudinal info isn't lost.
    """
    m = re.match(r"PPMI_(\d+)_(\d{8})_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    patno, yyyymmdd, mod = m.groups()
    scan_date = f"{yyyymmdd[:4]}-{yyyymmdd[4:6]}-{yyyymmdd[6:]}"
    return {
        "lookup_cols": {"subject_id": patno},
        "session_id": pd.NA,
        "modality": mod,
        "override": {"scan_date": scan_date},  # preserve actual scan date in output
    }


def parse_oasis1(filename: str) -> dict | None:
    """OAS1_{NNNN}_MR{n}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"(OAS1_\d+_MR\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    return {"lookup_cols": {"scan_id": m.group(1)}, "session_id": pd.NA, "modality": m.group(2)}


def parse_oasis2(filename: str) -> dict | None:
    """OAS2_{NNNN}_MR{n}_{mod}_preprocessed.nii.gz (session baked into scan_id)"""
    m = re.match(r"(OAS2_\d+_MR\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    scan_id = m.group(1)
    session = scan_id.rsplit("_", 1)[1]  # MR1, MR2, MR3
    return {"lookup_cols": {"scan_id": scan_id}, "session_id": session, "modality": m.group(2)}


def parse_ixi(filename: str) -> dict | None:
    """IXI{NNN}-{Site}-{scan}_{mod}_preprocessed.nii.gz — master key is IXI{NNN}"""
    m = re.match(r"(IXI\d+)-[^_]+_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    return {"lookup_cols": {"subject_id": m.group(1)}, "session_id": pd.NA, "modality": m.group(2)}


def parse_schizo(filename: str) -> dict | None:
    """SCHIZO_{eid}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"SCHIZO_([A-Z0-9]+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    return {"lookup_cols": {"subject_id": m.group(1)}, "session_id": pd.NA, "modality": m.group(2)}


def parse_stanford(filename: str) -> dict | None:
    """Stanford_{Patient-NN}_{mod}_preprocessed.nii.gz — no master CSV"""
    m = re.match(r"Stanford_(Patient-\d+)_(t1gd|flair|t1|t2)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    return {
        "lookup_cols": {"subject_id": m.group(1)},
        "session_id": pd.NA,
        "modality": m.group(2),
        "placeholder": True,  # signals no master join, emit raw row
    }


def parse_hcp_ya(filename: str) -> dict | None:
    """HCP_{subject}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"HCP_(\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    return {"lookup_cols": {"subject_id": m.group(1)}, "session_id": pd.NA, "modality": m.group(2)}


def parse_bgsp(filename: str) -> dict | None:
    """BGSP_sub-{NNNN}_ses-{NN}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"BGSP_(sub-\d+)_ses-(\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    sub, ses, mod = m.groups()
    return {"lookup_cols": {"subject_id": sub}, "session_id": f"ses-{ses}", "modality": mod}


def parse_adhd200(filename: str) -> dict | None:
    """sub-{NNNNNNN}_{mod}_preprocessed.nii.gz — master has ScanDir ID as
    int (no leading zeros)."""
    m = re.match(r"(sub-\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    sub_num = m.group(1).replace("sub-", "").lstrip("0") or "0"
    return {"lookup_cols": {"subject_id": sub_num},
            "session_id": pd.NA, "modality": m.group(2)}


def parse_corr(filename: str) -> dict | None:
    """sub-{NNNNNNN}_{mod}_preprocessed.nii.gz — master CoRR has SUBID as int-string"""
    m = re.match(r"(sub-\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    # Master CoRR has subject_id = str(SUBID) with NO leading zeros (e.g. "3001")
    # File has zero-padded 7-digit (e.g. "sub-0003001") → strip "sub-" and leading zeros
    sub_num = m.group(1).replace("sub-", "").lstrip("0") or "0"
    return {"lookup_cols": {"subject_id": sub_num},
            "session_id": pd.NA, "modality": m.group(2)}


def parse_fcon1000(filename: str) -> dict | None:
    """FCON_sub-{NNNNN}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"FCON_(sub-\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    sub = m.group(1)
    # Master FCON1000 has subject_id as-is from per-site TSV (could be "sub-00031" or "31")
    return {"lookup_cols": {"subject_id": sub},
            "session_id": pd.NA, "modality": m.group(2),
            "fallback_lookup_cols": {"subject_id": sub.replace("sub-", "").lstrip("0") or "0"}}


def parse_abide(filename: str) -> dict | None:
    """[ABIDE_[_cleaned_]]sub-{NNNNNNN}_{mod}_preprocessed.nii.gz"""
    # Strip optional intermediate prefixes (ABIDE_, _cleaned_, _com_aligned_)
    m = re.search(r"(sub-\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    sub = m.group(1)
    # Master ABIDE has subject_id = str(int SUB_ID) (e.g. "51472")
    sub_num = sub.replace("sub-", "").lstrip("0") or "0"
    return {"lookup_cols": {"subject_id": sub_num},
            "session_id": pd.NA, "modality": m.group(2)}


def parse_abide2(filename: str) -> dict | None:
    """ABIDE[_ intermediate_]sub-{NNNNN}_ses-{N}[_run-{N}|_acq-...]_{mod}_preprocessed.nii.gz"""
    m = re.search(r"(sub-\d+)(?:_ses-(\d+))?.*_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    sub, ses, mod = m.group(1), m.group(2), m.group(3)
    sub_num = sub.replace("sub-", "").lstrip("0") or "0"
    return {"lookup_cols": {"subject_id": sub_num},
            "session_id": f"ses-{ses}" if ses else pd.NA,
            "modality": mod}


def parse_hbn(filename: str) -> dict | None:
    """sub-NDAR{XXXXXXXX}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"(sub-NDAR\w+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    return {"lookup_cols": {"subject_id": m.group(1)},
            "session_id": pd.NA, "modality": m.group(2)}


def parse_nki(filename: str) -> dict | None:
    """sub-{A00NNNNNNN}_ses-{BAS1|BAS2|FLU1|...}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"(sub-A\d+)_ses-(\w+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    sub, ses, mod = m.groups()
    sid = sub.replace("sub-", "")
    return {"lookup_cols": {"subject_id": sid}, "session_id": ses, "modality": mod}


def parse_nki2(filename: str) -> dict | None:
    """sub-{A00NNNNNNN}_ses-MRI{N}_{mod}_preprocessed.nii.gz"""
    m = re.match(r"(sub-A\d+)_ses-(MRI\d+)_(t1|t2|flair)_preprocessed\.nii\.gz$", filename)
    if not m:
        return None
    sub, ses, mod = m.groups()
    sid = sub.replace("sub-", "")
    return {"lookup_cols": {"subject_id": sid}, "session_id": ses, "modality": mod}


def parse_upenn(filename: str) -> dict | None:
    """UPenn_UPENN-GBM-{NNNNN}_{NN}_{modality}_preprocessed.nii.gz
    Example: UPenn_UPENN-GBM-00001_11_t1_preprocessed.nii.gz"""
    m = re.match(
        r"UPenn_UPENN-GBM-(\d+)_(\d+)_(t1|t2|t1gd|t1c|flair)_preprocessed\.nii\.gz$",
        filename, re.IGNORECASE,
    )
    if not m:
        return None
    num, sess, mod = m.group(1), m.group(2), m.group(3).lower()
    return {"lookup_cols": {"subject_id": f"UPENN-GBM-{num}_{sess}"},
            "session_id": pd.NA, "modality": mod}


def parse_ucsf(filename: str) -> dict | None:
    """UCSF-PDGM-{NNNN}_{modality}_preprocessed.nii.gz
    Example: UCSF-PDGM-0004_t1_preprocessed.nii.gz
    Files have 4-digit zero-padding (e.g. 0465); master uses 3-digit (e.g. 465)."""
    m = re.match(
        r"UCSF-PDGM-(\d+)_(t1|t1c|t2|flair)_preprocessed\.nii\.gz$",
        filename, re.IGNORECASE,
    )
    if not m:
        return None
    num, mod = m.group(1), m.group(2).lower()
    master_sid = f"UCSF-PDGM-{int(num):03d}"
    return {"lookup_cols": {"subject_id": master_sid}, "session_id": pd.NA, "modality": mod}


PARSERS = {
    "parse_adni": parse_adni, "parse_adni_old": parse_adni_old, "parse_ppmi": parse_ppmi,
    "parse_oasis1": parse_oasis1, "parse_oasis2": parse_oasis2, "parse_ixi": parse_ixi,
    "parse_schizo": parse_schizo, "parse_stanford": parse_stanford, "parse_hcp_ya": parse_hcp_ya,
    "parse_bgsp": parse_bgsp, "parse_adhd200": parse_adhd200, "parse_corr": parse_corr,
    "parse_fcon1000": parse_fcon1000, "parse_abide": parse_abide, "parse_abide2": parse_abide2,
    "parse_hbn": parse_hbn, "parse_nki": parse_nki, "parse_nki2": parse_nki2,
    "parse_upenn": parse_upenn, "parse_ucsf": parse_ucsf,
}


# =============================================================================
# Main build
# =============================================================================


def _dataset_key_from_path(path: Path) -> str | None:
    """Walk from immediate parent outward, return first match in DATASET_MAP.
    Stops after 4 levels to avoid matching ancestor names like /oak/stanford/ → 'stanford'."""
    for i, parent in enumerate(path.parents):
        if i >= 4:
            break
        if parent.name in DATASET_MAP:
            return parent.name
    return None


def _load_master(metadata_root: Path, relpath: str) -> pd.DataFrame:
    p = metadata_root / relpath
    if not p.exists():
        raise FileNotFoundError(f"Master not found: {p}")
    df = pd.read_csv(p, low_memory=False)
    # Force all join columns to strings for robust matching
    for col in ("subject_id", "scan_id", "scan_date"):
        if col in df.columns:
            df[col] = df[col].astype(str)
    return df


def build_volumes(files_txt: Path, metadata_root: Path, output_csv: Path):
    paths = [Path(line.strip()) for line in files_txt.read_text().splitlines() if line.strip()]
    log.info(f"Loaded {len(paths)} file paths from {files_txt}")

    # Group paths by dataset key
    by_dataset: dict[str, list[Path]] = {}
    skipped = 0
    for p in paths:
        ds = _dataset_key_from_path(p)
        if ds is None:
            skipped += 1
            continue
        by_dataset.setdefault(ds, []).append(p)
    log.info(f"{skipped} files skipped (no recognized dataset in path)")

    # Cache masters (some datasets share a master — e.g. adni_full + adni_old)
    master_cache: dict[str, pd.DataFrame] = {}
    all_rows = []
    no_match = []

    for ds_key, ds_paths in sorted(by_dataset.items()):
        relpath, parser_name = DATASET_MAP[ds_key]
        parser = PARSERS[parser_name]
        if relpath is None:
            master = None  # Stanford — no clinical master
        else:
            if relpath not in master_cache:
                try:
                    master_cache[relpath] = _load_master(metadata_root, relpath)
                except FileNotFoundError as e:
                    log.warning(f"[{ds_key}] {e}, emitting rows without clinical data")
                    master_cache[relpath] = None
            master = master_cache[relpath]

        matched = 0
        for p in ds_paths:
            parsed = parser(p.name)
            if parsed is None:
                no_match.append(str(p))
                continue

            row = {
                "file_path": str(p),
                "dataset": ds_key,
                "modality": parsed["modality"],
                "session_id": parsed["session_id"],
            }

            if master is None:
                # No clinical master for this dataset — emit the row with
                # whatever parsed lookup info we have, dx stays NaN unless
                # the parser hardcodes a default.
                row.update(parsed["lookup_cols"])
                if ds_key == "stanford":
                    # Visual inspection of representative samples showed
                    # imaging characteristics consistent with WHO grade IV
                    # glioblastoma (large, contrast-enhancing, necrotic,
                    # mass effect). No per-subject pathology available.
                    row["dx"] = "GBM"
                    row["dx_detail"] = "Stanford tumor cohort (weak GBM label from cohort-level imaging review)"
                    row["tumor_grade"] = 4  # weak cohort-level label
                    # IDH/MGMT stay NaN (no pathology)
                elif ds_key == "bgsp":
                    row["dx"] = "HC"
                    row["dx_detail"] = "BGSP healthy (clinical pending DUA)"
                all_rows.append(row)
                matched += 1
                continue

            # Look up in master
            mask = pd.Series(True, index=master.index)
            for col, val in parsed["lookup_cols"].items():
                if col not in master.columns:
                    mask = pd.Series(False, index=master.index)
                    break
                mask &= (master[col] == str(val))
            hits = master[mask]

            # Fallback lookup (used by FCON1000 for sub-XXX vs numeric mismatch)
            if hits.empty and "fallback_lookup_cols" in parsed:
                mask2 = pd.Series(True, index=master.index)
                for col, val in parsed["fallback_lookup_cols"].items():
                    mask2 &= (master[col] == str(val))
                hits = master[mask2]

            if hits.empty:
                # Still emit the row so pretraining can see the file — just
                # without clinical fields. Label-dependent tasks will filter on NaN.
                for k, v in parsed["lookup_cols"].items():
                    row[k] = v
                row["dx"] = pd.NA
                all_rows.append(row)
                no_match.append(str(p))
                continue

            # Inherit master clinical columns (dedupe if multiple matches — take first)
            master_row = hits.iloc[0].to_dict()
            # File-side fields override master-side (modality, session, scan_id)
            master_row.update(row)
            # Parser-specified overrides (e.g. PPMI: actual scan_date from filename)
            if "override" in parsed:
                master_row.update(parsed["override"])
            all_rows.append(master_row)
            matched += 1

        log.info(f"[{ds_key}] {matched}/{len(ds_paths)} files matched to clinical data")

    df = pd.DataFrame(all_rows)

    # Canonicalize dataset names (normalize → upenn, UCSF-PDGM-v3 → ucsf, etc.)
    df["dataset"] = df["dataset"].replace(DATASET_CANONICAL_NAME)

    # Collapse dx aliases (CN → HC)
    df["dx"] = df["dx"].replace(DX_ALIASES)

    # Collapse modality aliases (t1gd → t1c)
    df["modality"] = df["modality"].replace(MODALITY_ALIASES)

    # Normalize APOE genotype format (ADNI '34.0' → 'E3/E4')
    if "apoe" in df.columns:
        df["apoe"] = df["apoe"].apply(_normalize_apoe)

    # Put the core columns first
    lead = ["file_path", "dataset", "subject_id", "session_id", "modality", "scan_date",
            "age_at_scan", "sex", "dx", "dx_detail", "site"]
    ordered = [c for c in lead if c in df.columns] + [c for c in df.columns if c not in lead]
    df = df[ordered]

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_csv, index=False)
    log.info(f"Wrote {len(df):,} volume rows to {output_csv}")

    # Report no-match paths (limited)
    if no_match:
        log.warning(f"{len(no_match)} files had no clinical match. Sample:")
        for p in no_match[:5]:
            log.warning(f"  {p}")

    # Brief summary
    log.info(f"Datasets: {df['dataset'].nunique()}, rows: {len(df):,}, "
             f"subjects: {df['subject_id'].nunique():,}")
    log.info("Modality breakdown:")
    for mod, n in df["modality"].value_counts().items():
        log.info(f"  {mod}: {n:,}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Build volumes.csv: one row per preprocessed .nii.gz file"
    )
    parser.add_argument("--files", type=Path, required=True,
                        help="Text file with one path per line (from find on Sherlock)")
    parser.add_argument("--metadata-root", type=Path, default=Path("data/metadata"),
                        help="Root of per-dataset *_master.csv files")
    parser.add_argument("--output", "-o", type=Path, required=True)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")

    build_volumes(args.files, args.metadata_root, args.output)


if __name__ == "__main__":
    main()
