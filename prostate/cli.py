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


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="python -m prostate")
    p.add_argument("command", choices=["discover", "organize"])
    p.add_argument("--config", type=Path, default=None,
                   help="path to preprocess.yaml (default: configs/preprocess.yaml)")
    p.add_argument("--limit-per-source", type=int, default=None,
                   help="cap each source to N records (smoke testing)")
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
    else:
        cmd_organize(cfg, args.limit_per_source)


if __name__ == "__main__":
    main()
