"""Multimodal MRI ↔ RNA-seq curation pipeline.

This package builds reproducible subject×gene expression matrices paired with
the brain-mri-preprocessing imaging cohorts. The pipeline is config-driven
(see configs/cohorts.yaml and configs/curate.yaml) and decomposed into
idempotent steps that can be run individually or end-to-end.

Quick start:
    python -m multimodal all           # full pipeline
    python -m multimodal inventory     # imaging subjects
    python -m multimodal query         # GDC search → manifests
    python -m multimodal download      # pull TSVs
    python -m multimodal matrix        # build subject×gene matrix
    python -m multimodal pairs         # imaging↔expression join
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
PACKAGE_ROOT = Path(__file__).resolve().parent
CONFIGS_DIR = PACKAGE_ROOT / "configs"
DATA_ROOT = REPO_ROOT / "data" / "multimodal"

__all__ = ["REPO_ROOT", "PACKAGE_ROOT", "CONFIGS_DIR", "DATA_ROOT"]
