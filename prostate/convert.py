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
    use_symlinks: bool = False,
) -> None:
    """Make sure `dst_path` is a real NIfTI file on disk.

    Default behavior is a real copy (or SimpleITK roundtrip), so the output
    directory is portable — you can rsync it to another machine without the
    targets going stale. Pass `use_symlinks=True` to symlink already-NIfTI
    sources for a transient / local-only run.
    """
    if dst_path.exists() or dst_path.is_symlink():
        return
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    if src_format == "nii.gz":
        if use_symlinks:
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
