"""CLI entry point: `python -m prostate <command>`.

Commands:
    inventory   Scan dcm2niix output, group into per-patient (T2/ADC/DWI) records, print summary.
    organize    Materialize a flat, VAE-loadable dir of single-channel volumes,
                one per (patient, modality). Symlinks by default.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import yaml

from . import REPO_ROOT, CONFIGS_DIR
from . import patients as patients_mod
from . import organize as organize_mod

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


def cmd_organize(cfg: dict, limit: int | None = None) -> None:
    nifti_root = _path(cfg["paths"]["nifti_root"])
    output_root = _path(cfg["paths"]["output_root"])
    modalities = cfg["output"]["modalities"]
    template = cfg["output"]["filename_template"]
    copy = bool(cfg["output"].get("copy", False))

    records = patients_mod.load_patients(nifti_root, cfg["series"])
    # Include any patient with at least one of the target modalities — drop
    # only the truly empty records.
    eligible = [r for r in records.values()
                if any(getattr(r, m, None) is not None for m in modalities)]
    log.info(f"Organizing {len(eligible)} patients (of {len(records)} total) → {output_root}  "
             f"({'copy' if copy else 'symlink'})")
    if limit:
        eligible = eligible[:limit]
        log.info(f"  (limit={limit})")

    n_ok = n_files = n_fail = 0
    per_modality = {m: 0 for m in modalities}
    for i, rec in enumerate(eligible, 1):
        try:
            out = organize_mod.organize_patient(rec, output_root, modalities, template, copy=copy)
            n_ok += 1
            n_files += len(out)
            for p in out:
                for m in modalities:
                    if p.name.endswith(f"_{m}.nii.gz"):
                        per_modality[m] += 1
                        break
        except Exception:
            log.exception(f"{rec.patient_id}: failed")
            n_fail += 1
        if i % 50 == 0 or i == len(eligible):
            log.info(f"  progress: {i}/{len(eligible)}  patients_ok={n_ok} fail={n_fail} files={n_files}")
    log.info(f"done: {n_ok} patients, {n_files} files, {n_fail} failed")
    for m, c in per_modality.items():
        log.info(f"  {m}: {c} files")
    log.info(f"output dir: {output_root}")


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="python -m prostate")
    p.add_argument("command", choices=["inventory", "organize"])
    p.add_argument("--config", type=Path, default=None,
                   help="path to preprocess.yaml (default: configs/preprocess.yaml)")
    p.add_argument("--limit", type=int, default=None,
                   help="process only the first N complete patients")
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
        cmd_organize(cfg, limit=args.limit)


if __name__ == "__main__":
    main()
