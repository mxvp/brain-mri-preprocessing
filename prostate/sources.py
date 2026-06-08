"""Per-source walkers that emit a stream of `VolumeRecord` for every usable
MRI volume across the four prostate datasets.

A single common record format means the downstream organizer / splitter doesn't
care which source a file came from. New cohort = add a walker, register it
in WALKERS.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator


@dataclass
class VolumeRecord:
    source: str          # 'picai' | 'tcia_biopsy' | 'prostate158' | 'promise12'
    patient_id: str      # canonical per-source patient identifier
    modality: str        # 't2' | 'adc' | 'dwi'
    view: str            # 'axial' | 'sagittal' | 'coronal'
    src_path: Path
    src_format: str      # 'nii.gz' | 'mha' | 'mhd'


# ---- PI-CAI -----------------------------------------------------------------
#
# Layout:
#     picai/images/<patient>/<patient>_<study>_<suffix>.mha
#
# Suffixes:
#     t2w  → axial T2     adc → axial ADC      hbv → axial high-b DWI
#     sag  → sagittal T2  cor → coronal T2
#
# Multi-study patients exist (~0.5% of cases). Treat each (patient, study)
# pair as a separate sample — same as how the brain side treats multiple
# scans of the same subject.

_PICAI_SUFFIX_TO_MODALITY_VIEW = {
    "t2w": ("t2", "axial"),
    "sag": ("t2", "sagittal"),
    "cor": ("t2", "coronal"),
    "adc": ("adc", "axial"),
    "hbv": ("dwi", "axial"),
}


def walk_picai(root: Path) -> Iterator[VolumeRecord]:
    pat = re.compile(r"^(\d+)_(\d+)_([a-z0-9]+)\.mha$", re.IGNORECASE)
    for patient_dir in sorted(p for p in root.iterdir()
                              if p.is_dir() and not p.name.startswith(".")):
        for f in sorted(patient_dir.glob("*.mha")):
            m = pat.match(f.name)
            if not m:
                continue
            pid, study, suffix = m.groups()
            mv = _PICAI_SUFFIX_TO_MODALITY_VIEW.get(suffix.lower())
            if mv is None:
                continue
            modality, view = mv
            yield VolumeRecord(
                source="picai",
                patient_id=f"{pid}_{study}",
                modality=modality,
                view=view,
                src_path=f,
                src_format="mha",
            )


# ---- TCIA Prostate-MRI-US-Biopsy ------------------------------------------
#
# Layout (after dcm2niix):
#     <root>/<patient>_<series-description>_<series-num>.nii.gz
#     <root>/<patient>_<series-description>_<series-num>.json   (sidecar)
#
# Classify by JSON SeriesDescription. All sequences are axial (no sag/cor in
# this collection). Multi-acquisition patients keep the longest stack.

_TCIA_BIOPSY_PATTERNS = {
    "t2":  [r"t2.*prostate", r"t2_spc", r"ax_t2", r"t2_tra"],
    "adc": [r"adc"],
    "dwi": [r"calc_bval", r"high_b", r"tracew", r"b\d{3,4}\b"],
}
_TCIA_PATIENT_RE = re.compile(r"(Prostate-MRI-US-Biopsy-\d+)")


def _tcia_classify(desc: str) -> str | None:
    desc = (desc or "").lower()
    for mod, patterns in _TCIA_BIOPSY_PATTERNS.items():
        if any(re.search(p, desc, re.IGNORECASE) for p in patterns):
            return mod
    return None


def walk_tcia_biopsy(root: Path) -> Iterator[VolumeRecord]:
    # Same dedupe-by-longest logic as the old patients.load_patients: per
    # (patient, modality) keep the acquisition with the most slices.
    import nibabel as nib

    best: dict[tuple[str, str], tuple[int, Path]] = {}
    for jpath in sorted(root.glob("*.json")):
        try:
            meta = json.loads(jpath.read_text())
        except Exception:
            continue
        nii = jpath.with_suffix("").with_suffix(".nii.gz")
        if not nii.exists():
            continue
        pat_match = _TCIA_PATIENT_RE.search(jpath.name)
        if not pat_match:
            continue
        pid = pat_match.group(1)
        modality = _tcia_classify(meta.get("SeriesDescription", ""))
        if modality is None:
            continue
        try:
            n_z = int(nib.load(nii).shape[2])
        except Exception:
            continue
        key = (pid, modality)
        if key not in best or n_z > best[key][0]:
            best[key] = (n_z, nii)

    for (pid, modality), (_, path) in best.items():
        yield VolumeRecord(
            source="tcia_biopsy",
            patient_id=pid,
            modality=modality,
            view="axial",
            src_path=path,
            src_format="nii.gz",
        )


# ---- Prostate158 ----------------------------------------------------------
#
# Layout:
#     prostate158_train/train/<patient>/{t2,adc,dwi}.nii.gz
#     prostate158_test/test/<patient>/{t2,adc,dwi}.nii.gz
#
# All axial. Patient IDs are 3-digit strings; collisions across train/test are
# rare but possible — prefix with the split name to be safe.

_PROSTATE158_MODALITIES = ("t2", "adc", "dwi")


def walk_prostate158(root: Path) -> Iterator[VolumeRecord]:
    splits = {
        "train": root / "prostate158_train" / "train",
        "test":  root / "prostate158_test" / "test",
    }
    for split_name, split_dir in splits.items():
        if not split_dir.is_dir():
            continue
        for patient_dir in sorted(p for p in split_dir.iterdir()
                                  if p.is_dir() and not p.name.startswith(".")):
            for mod in _PROSTATE158_MODALITIES:
                f = patient_dir / f"{mod}.nii.gz"
                if not f.exists():
                    continue
                yield VolumeRecord(
                    source="prostate158",
                    patient_id=f"{split_name}_{patient_dir.name}",
                    modality=mod,
                    view="axial",
                    src_path=f,
                    src_format="nii.gz",
                )


# ---- PROMISE12 ------------------------------------------------------------
#
# Layout:
#     TrainingDataPart{1,2,3}/Case##.{mhd,raw}              (image)
#     TrainingDataPart{1,2,3}/Case##_segmentation.{mhd,raw} (label — skip)
#     TestData/Case##.{mhd,raw}                              (image)
#
# T2 only by design (it's a prostate segmentation challenge).
# Prefix patient ID with the part name to disambiguate across parts.

_PROMISE12_PARTS = ("TrainingDataPart1", "TrainingDataPart2",
                    "TrainingDataPart3", "TestData")


def walk_promise12(root: Path) -> Iterator[VolumeRecord]:
    pat = re.compile(r"^(Case\d+)\.mhd$")
    for part in _PROMISE12_PARTS:
        part_dir = root / part
        if not part_dir.is_dir():
            continue
        for f in sorted(part_dir.glob("Case*.mhd")):
            if "_segmentation" in f.name:
                continue
            m = pat.match(f.name)
            if not m:
                continue
            yield VolumeRecord(
                source="promise12",
                patient_id=f"{part}_{m.group(1)}",
                modality="t2",
                view="axial",
                src_path=f,
                src_format="mhd",
            )


# ---- registry --------------------------------------------------------------

WALKERS = {
    "picai":       walk_picai,
    "tcia_biopsy": walk_tcia_biopsy,
    "prostate158": walk_prostate158,
    "promise12":   walk_promise12,
}
