"""Convert Analyze format (.hdr/.img) to NIfTI (.nii.gz).

Needed for OASIS-1/2 raw data which ships in Analyze format.

Usage:
    python convert_analyze.py input.hdr output.nii.gz
    python convert_analyze.py input_dir/ output_dir/ --batch
"""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np


def convert_volume(input_path: Path, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    img = nib.load(input_path)
    data = img.get_fdata()

    # Drop trailing singleton dimensions (OASIS has shape 256x256x128x1)
    data = data.squeeze()

    nii = nib.Nifti1Image(data.astype(np.float32), img.affine, img.header)
    nib.save(nii, output_path)
    print(f"Converted: {input_path.name} -> {output_path.name} (shape: {data.shape})")


def convert_batch(input_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(input_dir.rglob("*.hdr"))
    if not files:
        print(f"No .hdr files found in {input_dir}")
        return

    print(f"Found {len(files)} Analyze volumes in {input_dir}")
    for i, f in enumerate(files):
        out = output_dir / f"{f.stem}.nii.gz"
        print(f"[{i+1}/{len(files)}] ", end="")
        try:
            convert_volume(f, out)
        except Exception as e:
            print(f"  FAILED: {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Analyze (.hdr/.img) to NIfTI (.nii.gz)")
    parser.add_argument("input", type=Path, help="Input .hdr file or directory (with --batch)")
    parser.add_argument("output", type=Path, help="Output .nii.gz file or directory (with --batch)")
    parser.add_argument("--batch", action="store_true", help="Convert all .hdr files in input dir recursively")
    args = parser.parse_args()

    if args.batch:
        convert_batch(args.input, args.output, )
    else:
        convert_volume(args.input, args.output)
