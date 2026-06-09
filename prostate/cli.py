"""CLI entry: `python -m prostate <command>`.

Commands:
    discover    Walk every enabled source, print volume counts per
                (source, modality, view). No disk writes.
    organize    Materialize the harmonized layout under cfg.paths.output_root,
                emit manifest.csv with train/val splits.
"""
from __future__ import annotations

import argparse
import logging
from collections import Counter
from pathlib import Path

import yaml

from . import REPO_ROOT, CONFIGS_DIR
from . import organize as organize_mod
from . import preprocess as preprocess_mod
from . import qc as qc_mod
from . import sources as sources_mod

log = logging.getLogger("prostate")


def _load_cfg(path: Path | None) -> dict:
    path = path or (CONFIGS_DIR / "preprocess.yaml")
    with open(path) as f:
        return yaml.safe_load(f)


def cmd_discover(cfg: dict) -> None:
    grand_total = 0
    for source_name, source_cfg in cfg["sources"].items():
        if not source_cfg.get("enabled", True):
            continue
        root = Path(source_cfg["root"]).expanduser()
        if not root.exists():
            log.warning(f"{source_name}: root {root} not found")
            continue
        records = list(sources_mod.WALKERS[source_name](root))
        log.info(f"\n{source_name}: {len(records)} volumes  ({root})")
        counts = Counter((r.modality, r.view) for r in records)
        for (mod, view), n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
            log.info(f"  {mod}_{view:<10s} {n}")
        grand_total += len(records)
    log.info(f"\nTotal across enabled sources: {grand_total}")


def cmd_organize(cfg: dict, limit_per_source: int | None) -> None:
    organize_mod.organize_all(cfg, limit_per_source=limit_per_source)


def cmd_qc(cfg: dict, no_canonicalize: bool) -> None:
    output_root = Path(cfg["paths"]["output_root"]).expanduser().resolve()
    qc_cfg = cfg.get("qc", {})
    qc_mod.run(
        output_root=output_root,
        target_orientation=qc_cfg.get("target_orientation", "LPS"),
        canonicalize=not no_canonicalize,
        max_workers=int(cfg.get("workers", 8)),
    )


def cmd_preprocess(cfg: dict) -> None:
    organized_root = Path(cfg["paths"]["output_root"]).expanduser().resolve()
    pp_cfg = cfg.get("preprocess", {})
    pp_root = Path(pp_cfg.get(
        "output_root",
        str(organized_root.parent / (organized_root.name + "_preprocessed"))
    )).expanduser().resolve()
    spacing = tuple(pp_cfg.get("spacing", [0.5, 0.5, 3.0]))
    matrix  = tuple(pp_cfg.get("matrix_size", [160, 160, 20]))
    filters = pp_cfg.get("filters", {})
    preprocess_mod.run(
        organized_root=organized_root,
        preprocessed_root=pp_root,
        spacing_xyz=spacing,
        matrix_size_xyz=matrix,
        max_workers=int(cfg.get("workers", 8)),
        filters=filters,
    )


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="python -m prostate")
    p.add_argument("command", choices=["discover", "organize", "qc", "preprocess"])
    p.add_argument("--config", type=Path, default=None,
                   help="path to preprocess.yaml (default: configs/preprocess.yaml)")
    p.add_argument("--limit-per-source", type=int, default=None,
                   help="cap each source to N records (smoke testing)")
    p.add_argument("--no-canonicalize", action="store_true",
                   help="qc: skip the orientation rewrite, only add manifest columns")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)s  %(name)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    cfg = _load_cfg(args.config)
    if args.command == "discover":
        cmd_discover(cfg)
    elif args.command == "organize":
        cmd_organize(cfg, args.limit_per_source)
    elif args.command == "qc":
        cmd_qc(cfg, args.no_canonicalize)
    else:  # preprocess
        cmd_preprocess(cfg)


if __name__ == "__main__":
    main()
