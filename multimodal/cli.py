"""Command-line entry point for the multimodal pipeline.

    python -m multimodal <command> [options]

Commands (run in this order, or just `all`):
    inventory   Build the imaging-subject inventory CSV from cohort configs.
    query       Hit GDC for paired RNA-seq → per-project manifests.
    download    Pull all manifest TSVs into data/multimodal/raw_tsvs/<project>/.
    matrix      Parse TSVs → counts / tpm / gene_meta / sample_meta parquet.
    pairs       Imaging-subjects × sample_meta join → pairs.csv.
    all         Run inventory → query → download → matrix → pairs.

All output paths are taken from configs/curate.yaml. Re-runs are idempotent.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

from . import REPO_ROOT, CONFIGS_DIR
from . import cohorts as cohorts_mod
from . import gdc
from . import download as download_mod
from . import matrix as matrix_mod
from . import pairs as pairs_mod

log = logging.getLogger("multimodal")


# ----- config helpers ------------------------------------------------------

def load_curate(path: Path | None = None) -> dict:
    path = path or (CONFIGS_DIR / "curate.yaml")
    with open(path) as f:
        return yaml.safe_load(f)


def _path(p: str | Path) -> Path:
    p = Path(p)
    return p if p.is_absolute() else (REPO_ROOT / p).resolve()


# ----- commands ------------------------------------------------------------

def cmd_inventory(curate: dict) -> Path:
    """Build the imaging-subject inventory."""
    cfg = cohorts_mod.load_cohorts_config()
    df = cohorts_mod.load_imaging_subjects(cfg)
    out = _path(curate["paths"]["matrix_root"]) / "imaging_subjects.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    log.info(f"inventory → {out}")
    log.info("\n" + df.groupby("cohort").size().to_string())
    return out


def cmd_query(curate: dict) -> dict[str, Path]:
    """Query GDC for RNA-seq paired with imaging subjects; write manifests."""
    cfg = cohorts_mod.load_cohorts_config()
    inv = cohorts_mod.load_imaging_subjects(cfg)
    man_root = _path(curate["paths"]["manifest_root"])
    man_root.mkdir(parents=True, exist_ok=True)

    out: dict[str, Path] = {}
    for key, spec in cfg["cohorts"].items():
        if key not in curate["include_cohorts"]:
            continue
        if not spec["expression"].get("has_public"):
            log.info(f"{key}: skipping (no public RNA-seq)")
            continue
        project = spec["expression"]["project_id"]
        sub_ids = inv.loc[inv["cohort"] == key, "submitter_id"].dropna().unique().tolist()
        if not sub_ids:
            log.warning(f"{key}: no imaging subjects to query")
            continue
        log.info(f"{key}: querying GDC ({project}) for {len(sub_ids)} subjects")
        df = gdc.search_rnaseq(project, sub_ids)
        # Save the full search results and a STAR-Counts-only download manifest.
        full = man_root / f"{key}_gdc_full.csv"
        df.to_csv(full, index=False)
        gene = gdc.filter_star_counts_expression(df)
        man = man_root / f"{key}_manifest.tsv"
        n = gdc.write_gdc_manifest(gene, man)
        size_gb = gene["file_size"].fillna(0).sum() / 1e9
        log.info(f"  {key}: {n} files, {gene['submitter_id'].nunique()} subjects, ~{size_gb:.2f} GB → {man}")
        out[key] = man
    return out


def cmd_download(curate: dict) -> dict[str, dict[str, int]]:
    """Download all manifest files into data/multimodal/raw_tsvs/<project>/."""
    cfg = cohorts_mod.load_cohorts_config()
    raw_root = _path(curate["paths"]["raw_tsv_root"])
    man_root = _path(curate["paths"]["manifest_root"])
    stats: dict[str, dict[str, int]] = {}
    for key in curate["include_cohorts"]:
        spec = cfg["cohorts"].get(key)
        if not spec or not spec["expression"].get("has_public"):
            continue
        project = spec["expression"]["project_id"]
        man_path = man_root / f"{key}_manifest.tsv"
        if not man_path.exists():
            log.warning(f"{key}: no manifest at {man_path}, run `query` first")
            continue
        manifest = pd.read_csv(man_path, sep="\t").rename(columns={
            "id": "file_id", "filename": "file_name", "md5": "md5sum", "size": "file_size",
        })
        dest = raw_root / project
        log.info(f"{key}: downloading {len(manifest)} files → {dest}")
        s = download_mod.download_files(manifest, dest)
        log.info(f"  {key}: {s}")
        stats[key] = s
    return stats


def cmd_matrix(curate: dict) -> Path:
    """Build the expression matrix from downloaded TSVs."""
    cfg = cohorts_mod.load_cohorts_config()
    raw_root = _path(curate["paths"]["raw_tsv_root"])
    man_root = _path(curate["paths"]["manifest_root"])
    out_root = _path(curate["paths"]["matrix_root"])
    manifests: dict[str, pd.DataFrame] = {}
    for key in curate["include_cohorts"]:
        spec = cfg["cohorts"].get(key)
        if not spec or not spec["expression"].get("has_public"):
            continue
        full = man_root / f"{key}_gdc_full.csv"
        if not full.exists():
            log.warning(f"{key}: no manifest at {full}, run `query` first")
            continue
        manifests[spec["expression"]["project_id"]] = pd.read_csv(full)
    tables = matrix_mod.build_expression_matrix(
        manifests, raw_root, sample_filter=curate["samples"]
    )
    matrix_mod.write_parquet(tables, out_root)
    log.info(f"matrix tables: " + ", ".join(f"{k} {v.shape}" for k, v in tables.items()))
    return out_root


def cmd_pairs(curate: dict) -> Path:
    """Join imaging inventory with expression sample_meta to make pairs.csv."""
    out_root = _path(curate["paths"]["matrix_root"])
    sm_path = out_root / "sample_meta.parquet"
    if not sm_path.exists():
        raise FileNotFoundError(f"missing {sm_path} — run `matrix` first")
    sample_meta = pd.read_parquet(sm_path)
    inv = cohorts_mod.load_imaging_subjects()
    pairs = pairs_mod.make_pairs(inv, sample_meta)
    pairs_path = out_root / "pairs.csv"
    pairs.to_csv(pairs_path, index=False)
    coverage = pairs_mod.coverage_summary(inv, pairs)
    coverage.to_csv(out_root / "coverage.csv", index=False)
    log.info(f"pairs → {pairs_path}")
    log.info("\n" + coverage.to_string(index=False))
    return pairs_path


def cmd_all(curate: dict) -> None:
    cmd_inventory(curate)
    cmd_query(curate)
    cmd_download(curate)
    cmd_matrix(curate)
    cmd_pairs(curate)


# ----- argparse glue -------------------------------------------------------

COMMANDS = {
    "inventory": cmd_inventory,
    "query":     cmd_query,
    "download":  cmd_download,
    "matrix":    cmd_matrix,
    "pairs":     cmd_pairs,
    "all":       cmd_all,
}


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(prog="python -m multimodal")
    p.add_argument("command", choices=list(COMMANDS))
    p.add_argument("--curate-config", type=Path, default=None,
                   help="path to curate.yaml (default: configs/curate.yaml)")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    curate = load_curate(args.curate_config)
    COMMANDS[args.command](curate)


if __name__ == "__main__":
    main()
