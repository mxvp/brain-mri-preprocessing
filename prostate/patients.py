"""Group flat dcm2niix output into per-patient (T2 / ADC / DWI) records.

dcm2niix produces files like:
    Prostate-MRI-US-Biopsy-0001_t2_spc_rst_axial_obl_Prostate_11.nii.gz
    Prostate-MRI-US-Biopsy-0001_..._ADC_DFC_MIX_5.nii.gz
    Prostate-MRI-US-Biopsy-0001_..._CALC_BVAL_DFC_MIX_7.nii.gz

Each `.nii.gz` has a matching `.json` sidecar with `SeriesDescription`. We
classify each volume by description using regex patterns from the config,
then build a per-patient record:

    {
      "patient_id": "Prostate-MRI-US-Biopsy-0001",
      "t2":  Path(...),     # required; patient skipped if missing
      "adc": Path(...) | None,
      "dwi": Path(...) | None,
    }

When a patient has multiple candidates for a series (e.g. two T2 acquisitions),
the one with the most slices wins — that's the closest thing to a "primary"
acquisition the dcm2niix output gives us.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

PATIENT_RE = re.compile(r"(Prostate-MRI-US-Biopsy-\d+)")


@dataclass
class PatientRecord:
    patient_id: str
    t2:  Path | None = None
    adc: Path | None = None
    dwi: Path | None = None
    extras: dict[str, list[Path]] = field(default_factory=dict)

    def is_complete(self) -> bool:
        return self.t2 is not None and self.adc is not None and self.dwi is not None

    def missing(self) -> list[str]:
        return [k for k in ("t2", "adc", "dwi") if getattr(self, k) is None]


def _classify(description: str, series_patterns: dict[str, list[str]]) -> str | None:
    """Return 't2' | 'adc' | 'dwi' | None based on the first matching regex."""
    if not description:
        return None
    desc = description.lower()
    for series, patterns in series_patterns.items():
        for pat in patterns:
            if re.search(pat, desc, re.IGNORECASE):
                return series
    return None


def _slices_in_nifti(nifti_path: Path) -> int:
    """Return the z-dim of the volume. Cheap header read, no pixel decode."""
    import nibabel as nib
    return int(nib.load(nifti_path).shape[2])


def load_patients(
    nifti_root: Path,
    series_patterns: dict[str, list[str]],
) -> dict[str, PatientRecord]:
    """Scan a flat dcm2niix output and return {patient_id: PatientRecord}."""
    nifti_root = Path(nifti_root)
    records: dict[str, PatientRecord] = {}

    # Read JSON sidecars to get SeriesDescription. Faster than peeking inside
    # the NIfTI header.
    for json_path in sorted(nifti_root.glob("*.json")):
        try:
            meta = json.loads(json_path.read_text())
        except Exception:
            continue
        nifti = json_path.with_suffix("").with_suffix(".nii.gz")
        if not nifti.exists():
            continue

        m = PATIENT_RE.search(json_path.name)
        if not m:
            continue
        pid = m.group(1)
        rec = records.setdefault(pid, PatientRecord(patient_id=pid))

        series = _classify(meta.get("SeriesDescription", ""), series_patterns)
        if series is None:
            rec.extras.setdefault("unclassified", []).append(nifti)
            continue

        # Multiple candidates? Keep the one with the most slices.
        current = getattr(rec, series)
        if current is None or _slices_in_nifti(nifti) > _slices_in_nifti(current):
            setattr(rec, series, nifti)

    return records


def summarize(records: Iterable[PatientRecord]) -> dict[str, int]:
    """Tally completeness across a set of patients."""
    records = list(records)
    counts = {"total": len(records), "complete": 0, "missing_t2": 0,
              "missing_adc": 0, "missing_dwi": 0, "missing_2_or_more": 0}
    for r in records:
        if r.is_complete():
            counts["complete"] += 1
        else:
            missing = r.missing()
            for k in missing:
                counts[f"missing_{k}"] += 1
            if len(missing) >= 2:
                counts["missing_2_or_more"] += 1
    return counts
