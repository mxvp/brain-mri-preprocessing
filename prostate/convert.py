"""File-format conversion: MetaImage (`.mha`, `.mhd`) → NIfTI (`.nii.gz`).

NIfTI is the brain-side convention, MONAI's LoadImaged reads it cleanly, and a
single extension simplifies downstream globs. Idempotent — pre-existing
outputs are left alone.
"""
from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)


def materialize_nifti(
    src_path: Path,
    src_format: str,
    dst_path: Path,
    symlink_if_already_nifti: bool = True,
) -> None:
    """Make sure `dst_path` is a NIfTI on disk pointing at the same data.

    Behavior by source format:
      'nii.gz' → symlink dst → src (no disk duplication)
      'mha' or 'mhd' → SimpleITK roundtrip; new .nii.gz at dst
    """
    if dst_path.exists() or dst_path.is_symlink():
        return
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    if src_format == "nii.gz":
        if symlink_if_already_nifti:
            dst_path.symlink_to(src_path.resolve())
        else:
            import shutil
            shutil.copy2(src_path, dst_path)
        return
    import SimpleITK as sitk
    img = sitk.ReadImage(str(src_path))
    sitk.WriteImage(img, str(dst_path))


def n_slices(nifti_path: Path) -> int:
    """Slice count along the z axis. Cheap — header-only read."""
    import nibabel as nib
    return int(nib.load(nifti_path).shape[2])
