"""Smoke tests for the per-source walkers — verify each finds and classifies
the right files in a fake tree mirroring each dataset's layout.
"""
import json
import tempfile
from pathlib import Path

import nibabel as nib
import numpy as np
import pytest

from prostate.sources import (
    walk_picai, walk_tcia_biopsy, walk_prostate158, walk_promise12,
)


def _empty_nifti(path: Path, slices: int = 20):
    path.parent.mkdir(parents=True, exist_ok=True)
    nib.save(nib.Nifti1Image(np.zeros((4, 4, slices), dtype=np.int16),
                             np.eye(4)), path)


def _empty_mha(path: Path, slices: int = 20):
    path.parent.mkdir(parents=True, exist_ok=True)
    import SimpleITK as sitk
    img = sitk.GetImageFromArray(np.zeros((slices, 4, 4), dtype=np.int16))
    sitk.WriteImage(img, str(path))


def _empty_mhd(path: Path, slices: int = 20):
    _empty_mha(path, slices=slices)  # SimpleITK writes mhd+raw from .mhd


# ---- PI-CAI -------------------------------------------------------------

def test_picai_finds_all_5_views():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        for suffix in ("t2w", "sag", "cor", "adc", "hbv"):
            _empty_mha(root / "10000" / f"10000_1000000_{suffix}.mha")
        records = list(walk_picai(root))
    by_key = {(r.modality, r.view) for r in records}
    assert by_key == {("t2", "axial"), ("t2", "sagittal"), ("t2", "coronal"),
                      ("adc", "axial"), ("dwi", "axial")}
    assert all(r.source == "picai" for r in records)
    assert all(r.patient_id == "10000_1000000" for r in records)


def test_picai_multi_study_patient():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        _empty_mha(root / "10000" / "10000_1000000_t2w.mha")
        _empty_mha(root / "10000" / "10000_1000001_t2w.mha")
        records = list(walk_picai(root))
    pids = sorted(r.patient_id for r in records)
    assert pids == ["10000_1000000", "10000_1000001"]


# ---- TCIA Biopsy --------------------------------------------------------

def _tcia_pair(root: Path, base: str, description: str, slices: int = 20):
    _empty_nifti(root / f"{base}.nii.gz", slices=slices)
    (root / f"{base}.json").write_text(json.dumps({
        "Modality": "MR", "SeriesDescription": description,
    }))


def test_tcia_biopsy_picks_longest_per_modality():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        _tcia_pair(root, "Prostate-MRI-US-Biopsy-0001_t2_short_3",
                   "t2_spc_rst_axial obl_Prostate", slices=10)
        _tcia_pair(root, "Prostate-MRI-US-Biopsy-0001_t2_long_5",
                   "t2_spc_rst_axial obl_Prostate", slices=30)
        _tcia_pair(root, "Prostate-MRI-US-Biopsy-0001_ADC_6",
                   "ep2d-advdiff-3Scan-4bval_spair_std_ADC", slices=20)
        records = list(walk_tcia_biopsy(root))
    by_mod = {r.modality: r for r in records}
    assert "t2" in by_mod and "long_5" in by_mod["t2"].src_path.name
    assert "adc" in by_mod
    assert all(r.view == "axial" for r in records)


# ---- Prostate158 --------------------------------------------------------

def test_prostate158_walks_train_and_test():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        for split, n_patients in [("prostate158_train/train", 2),
                                  ("prostate158_test/test", 1)]:
            for i in range(n_patients):
                pdir = root / split / f"{i:03d}"
                for mod in ("t2", "adc", "dwi"):
                    _empty_nifti(pdir / f"{mod}.nii.gz")
        records = list(walk_prostate158(root))
    assert len(records) == 3 * 3   # 3 patients × 3 modalities
    sources_split = {r.patient_id.split("_")[0] for r in records}
    assert sources_split == {"train", "test"}


# ---- PROMISE12 ----------------------------------------------------------

def test_promise12_skips_segmentation_files():
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        for part in ("TrainingDataPart1", "TrainingDataPart2",
                     "TrainingDataPart3", "TestData"):
            _empty_mhd(root / part / "Case00.mhd")
            _empty_mhd(root / part / "Case00_segmentation.mhd")  # should skip
        records = list(walk_promise12(root))
    assert len(records) == 4
    assert all(r.modality == "t2" and r.view == "axial" for r in records)
    assert all("segmentation" not in r.src_path.name for r in records)
