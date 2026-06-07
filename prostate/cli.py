"""CLI entry point: `python -m prostate <command>`.

Commands:
    inventory   Scan dcm2niix output, group into per-patient (T2/ADC/DWI) records, print summary.
    preprocess  Run the full pipeline for all complete patients.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import yaml

from . import REPO_ROOT, CONFIGS_DIR
from . import patients as patients_mod
from . import preprocess as preprocess_mod

log = logging.getLogger("prostate")


def _load_cfg(path: Path | None) -> dict:
    path = path or (CONFIGS_DIR / "preprocess.yaml")
    with open(path) as f:
        return yaml.safe_load(f)


def _path(p: str | Path) -> Path:
    p = Path(p)
    return p if p.is_absolute() else (REPO_ROOT / p).resolve()


def cmd_inventory(cfg: dict) -> None:
    nifti_root = _path(cfg["paths"]["nifti_root"])
    records = patients_mod.load_patients(nifti_root, cfg["series"])
    summary = patients_mod.summarize(records.values())
    log.info(f"patients: {summary['total']}")
    log.info(f"  complete (T2+ADC+DWI): {summary['complete']}")
    log.info(f"  missing T2:  {summary['missing_t2']}")
    log.info(f"  missing ADC: {summary['missing_adc']}")
    log.info(f"  missing DWI: {summary['missing_dwi']}")
    log.info(f"  missing ≥2 modalities: {summary['missing_2_or_more']}")


def cmd_preprocess(cfg: dict, limit: int | None = None) -> None:
    nifti_root = _path(cfg["paths"]["nifti_root"])
    output_root = _path(cfg["paths"]["output_root"])
    records = patients_mod.load_patients(nifti_root, cfg["series"])
    complete = [r for r in records.values() if r.is_complete()]
    log.info(f"Preprocessing {len(complete)} complete patients → {output_root}")
    if limit:
        complete = complete[:limit]
        log.info(f"  (limit={limit})")
    n_ok = n_skip = n_fail = 0
    for i, rec in enumerate(complete, 1):
        try:
            out = preprocess_mod.preprocess_patient(rec, cfg, output_root)
            if out is None:
                n_skip += 1
            else:
                n_ok += 1
        except Exception:
            log.exception(f"{rec.patient_id}: failed")
            n_fail += 1
        if i % 10 == 0 or i == len(complete):
            log.info(f"  progress: {i}/{len(complete)}  ok={n_ok} skip={n_skip} fail={n_fail}")
    log.info(f"done: ok={n_ok}  skip={n_skip}  fail={n_fail}")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="python -m prostate")
    p.add_argument("command", choices=["inventory", "preprocess"])
    p.add_argument("--config", type=Path, default=None,
                   help="path to preprocess.yaml (default: configs/preprocess.yaml)")
    p.add_argument("--limit", type=int, default=None,
                   help="preprocess only the first N complete patients (testing)")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    cfg = _load_cfg(args.config)
    if args.command == "inventory":
        cmd_inventory(cfg)
    else:
        cmd_preprocess(cfg, limit=args.limit)


if __name__ == "__main__":
    main()
