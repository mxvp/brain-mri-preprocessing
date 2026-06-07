"""Organize the dcm2niix output into a VAE-loadable directory layout.

No preprocessing. For each patient: pick one T2 / one ADC / one DWI (the
acquisition with the most slices wins, see patients.py) and link/copy each
into `output_root/<patient_id>_<modality>.nii.gz`. The VAE loader globs the
flat dir and treats each file as one training sample.

Why so minimal:
  - VAE generative training learns the data manifold; bias fields / scanner
    intensity quirks are part of that variability, not noise to remove.
  - Multi-modal stacking would force registration. We're keeping each
    modality as its own training sample, so DWI/ADC don't need to share
    T2's voxel grid.
  - VAE loader resizes + z-scores at load time. On-disk shape and intensity
    don't matter as long as background is reasonably ≈ 0 (true here).

Links by default (cheap, doesn't duplicate disk). Pass --copy if you need
self-contained outputs (e.g. for moving to another machine).
"""
from __future__ import annotations

import logging
import shutil
from pathlib import Path

from .patients import PatientRecord

log = logging.getLogger(__name__)


def organize_patient(
    record: PatientRecord,
    output_dir: Path,
    modalities: list[str],
    template: str,
    copy: bool = False,
) -> list[Path]:
    """Materialize whichever modality files this patient has in `output_dir`.

    Patients with only some of `modalities` are still processed — each present
    modality produces a link/copy; missing ones are silently skipped. Returns
    an empty list for patients with none of the requested modalities.
    Idempotent.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for mod in modalities:
        src = getattr(record, mod, None)
        if src is None:
            continue
        dst = output_dir / template.format(patient_id=record.patient_id, modality=mod)
        if dst.exists() or dst.is_symlink():
            written.append(dst)
            continue
        if copy:
            shutil.copy2(src, dst)
        else:
            dst.symlink_to(src.resolve())
        written.append(dst)
    return written
