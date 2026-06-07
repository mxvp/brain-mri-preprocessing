"""Smoke tests for the patient discovery / classification step."""
import json
import tempfile
from pathlib import Path

import numpy as np
import pytest

import nibabel as nib

from prostate.patients import load_patients, summarize


SERIES_PATTERNS = {
    "t2":  ["t2.*prostate", "t2_spc", "ax_t2"],
    "adc": ["adc"],
    "dwi": ["calc_bval", "high_b", r"b\d{3,4}\b"],
}


def _write_pair(dirpath: Path, base: str, description: str, slices: int = 20):
    """Write a (.nii.gz, .json) pair mimicking dcm2niix output."""
    img = nib.Nifti1Image(np.zeros((10, 10, slices), dtype=np.int16), np.eye(4))
    nib.save(img, dirpath / f"{base}.nii.gz")
    (dirpath / f"{base}.json").write_text(json.dumps({
        "Modality": "MR", "SeriesDescription": description,
    }))


def test_classify_t2_adc_dwi():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        _write_pair(td, "Prostate-MRI-US-Biopsy-0001_t2_spc_axial_obl_Prostate_11",
                    "t2_spc_rst_axial obl_Prostate")
        _write_pair(td, "Prostate-MRI-US-Biopsy-0001_ADC_5",
                    "ep2d-advdiff-3Scan-4bval_spair_std_ADC")
        _write_pair(td, "Prostate-MRI-US-Biopsy-0001_DWI_7",
                    "ep2d-advdiff-3Scan-4bval_spair_std_CALC_BVAL")
        recs = load_patients(td, SERIES_PATTERNS)
    assert len(recs) == 1
    rec = recs["Prostate-MRI-US-Biopsy-0001"]
    assert rec.is_complete()
    assert "_t2_" in rec.t2.name
    assert "ADC" in rec.adc.name
    assert "DWI" in rec.dwi.name


def test_multiple_candidates_pick_more_slices():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        _write_pair(td, "Prostate-MRI-US-Biopsy-0001_t2_short_3",
                    "t2_spc_rst_axial obl_Prostate", slices=10)
        _write_pair(td, "Prostate-MRI-US-Biopsy-0001_t2_long_5",
                    "t2_spc_rst_axial obl_Prostate", slices=30)
        recs = load_patients(td, SERIES_PATTERNS)
    assert "t2_long_5" in recs["Prostate-MRI-US-Biopsy-0001"].t2.name


def test_incomplete_patient_flagged():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        _write_pair(td, "Prostate-MRI-US-Biopsy-0001_t2_X_1",
                    "t2_spc_rst_axial obl_Prostate")
        recs = load_patients(td, SERIES_PATTERNS)
    s = summarize(recs.values())
    assert s["complete"] == 0 and s["total"] == 1
    assert s["missing_adc"] == 1 and s["missing_dwi"] == 1
