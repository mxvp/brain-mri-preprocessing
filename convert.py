"""Convert various neuroimaging formats to NIfTI (.nii.gz).

Supported: Analyze (.hdr/.img), MGZ (.mgz), DICOM (directories), MINC (.mnc).

Usage:
    python convert.py input.hdr output.nii.gz
    python convert.py input.mgz output.nii.gz
    python convert.py dicom_dir/ output.nii.gz
    python convert.py input_dir/ output_dir/ --batch
"""

import argparse
import logging
import subprocess
from pathlib import Path

import nibabel as nib
import numpy as np

log = logging.getLogger(__name__)


def _to_nifti(data: np.ndarray, affine: np.ndarray, output_path: Path):
    data = data.squeeze().astype(np.float32)
    nii = nib.Nifti1Image(data, affine)
    nib.save(nii, output_path)


def convert_analyze(input_path: Path, output_path: Path):
    img = nib.load(input_path)
    _to_nifti(img.get_fdata(), img.affine, output_path)
    log.info(f"Analyze: {input_path.name} -> {output_path.name}")


def convert_mgz(input_path: Path, output_path: Path):
    img = nib.load(input_path)
    _to_nifti(img.get_fdata(), img.affine, output_path)
    log.info(f"MGZ: {input_path.name} -> {output_path.name}")


def convert_minc(input_path: Path, output_path: Path):
    img = nib.load(input_path)
    _to_nifti(img.get_fdata(), img.affine, output_path)
    log.info(f"MINC: {input_path.name} -> {output_path.name}")


def convert_dicom(input_dir: Path, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["dcm2niix", "-z", "y", "-f", output_path.stem, "-o", str(output_path.parent), str(input_dir)],
            check=True,
            capture_output=True,
            text=True,
        )
        log.info(f"DICOM: {input_dir.name} -> {output_path.name}")
    except FileNotFoundError:
        raise RuntimeError("dcm2niix not found. Install it or load the module (module load dcm2niix).")


FORMAT_MAP = {
    ".hdr": convert_analyze,
    ".img": lambda inp, out: convert_analyze(inp.with_suffix(".hdr"), out),
    ".mgz": convert_mgz,
    ".mgh": convert_mgz,
    ".mnc": convert_minc,
}


def convert(input_path: Path, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if input_path.is_dir():
        # Assume DICOM directory
        dcm_files = list(input_path.glob("*.dcm")) + list(input_path.glob("*.DCM"))
        if not dcm_files:
            # Some DICOM dirs have no extension
            dcm_files = [f for f in input_path.iterdir() if f.is_file()]
        if not dcm_files:
            raise ValueError(f"No files found in {input_path}")
        convert_dicom(input_path, output_path)
        return

    suffix = input_path.suffix.lower()
    if suffix in (".gz",):
        # Already nifti
        log.info(f"Already NIfTI: {input_path.name}, skipping")
        return
    if suffix in (".nii",):
        # Uncompressed nifti, just compress
        img = nib.load(input_path)
        nib.save(img, output_path)
        log.info(f"Compressed: {input_path.name} -> {output_path.name}")
        return

    converter = FORMAT_MAP.get(suffix)
    if converter is None:
        raise ValueError(f"Unsupported format: {suffix} ({input_path.name})")
    converter(input_path, output_path)


def convert_batch(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    extensions = {".hdr", ".mgz", ".mgh", ".mnc"}
    files = sorted(f for f in input_dir.rglob("*") if f.suffix.lower() in extensions)

    # Also check for DICOM directories (dirs containing .dcm files)
    dicom_dirs = sorted(
        d for d in input_dir.rglob("*") if d.is_dir() and list(d.glob("*.dcm"))
    )

    total = len(files) + len(dicom_dirs)
    if total == 0:
        log.warning(f"No convertible files found in {input_dir}")
        return

    log.info(f"Found {len(files)} files + {len(dicom_dirs)} DICOM dirs in {input_dir}")
    failed = []
    for i, f in enumerate(files):
        out = output_dir / f"{f.stem}.nii.gz"
        try:
            convert(f, out)
        except Exception:
            log.exception(f"Failed: {f.name}")
            failed.append(f.name)

    for i, d in enumerate(dicom_dirs):
        out = output_dir / f"{d.name}.nii.gz"
        try:
            convert(d, out)
        except Exception:
            log.exception(f"Failed: {d.name}")
            failed.append(d.name)

    if failed:
        log.warning(f"{len(failed)}/{total} failed: {failed}")
    else:
        log.info(f"All {total} converted successfully")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    parser = argparse.ArgumentParser(description="Convert neuroimaging formats to NIfTI")
    parser.add_argument("input", type=Path, help="Input file/dir or directory (with --batch)")
    parser.add_argument("output", type=Path, help="Output .nii.gz file or directory (with --batch)")
    parser.add_argument("--batch", action="store_true", help="Convert all supported files recursively")
    args = parser.parse_args()

    if args.batch:
        convert_batch(args.input, args.output)
    else:
        convert(args.input, args.output)
