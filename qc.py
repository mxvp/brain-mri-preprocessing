"""Quick visual QC: render mid-slice PNGs for a directory of preprocessed volumes.

Outputs a grid of axial mid-slices for fast visual inspection.

Usage:
    python qc.py input_dir/ output.png
    python qc.py input_dir/ output.png --cols 10
"""

import argparse
from pathlib import Path

import nibabel as nib
import numpy as np


def render_midslice(nii_path: Path) -> np.ndarray:
    data = nib.load(nii_path).get_fdata().squeeze()
    mid = data.shape[2] // 2
    slc = data[:, :, mid].T
    # Normalize to 0-255 for display
    nz = slc[slc != 0]
    if len(nz) == 0:
        return np.zeros_like(slc, dtype=np.uint8)
    lo, hi = np.percentile(nz, [1, 99])
    slc = np.clip((slc - lo) / (hi - lo + 1e-8), 0, 1)
    return (slc * 255).astype(np.uint8)


def make_grid(images: list, cols: int = 10, pad: int = 2) -> np.ndarray:
    if not images:
        return np.zeros((1, 1), dtype=np.uint8)

    h, w = images[0].shape[:2]
    rows = (len(images) + cols - 1) // cols

    # Pad list to fill grid
    while len(images) < rows * cols:
        images.append(np.zeros((h, w), dtype=np.uint8))

    grid = np.zeros((rows * (h + pad) - pad, cols * (w + pad) - pad), dtype=np.uint8)
    for idx, img in enumerate(images):
        r, c = divmod(idx, cols)
        y, x = r * (h + pad), c * (w + pad)
        grid[y : y + h, x : x + w] = img[:h, :w]

    return grid


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="QC grid of mid-axial slices")
    parser.add_argument("input_dir", type=Path, help="Directory of .nii.gz files")
    parser.add_argument("output", type=Path, help="Output PNG path")
    parser.add_argument("--cols", type=int, default=10, help="Columns in grid (default: 10)")
    args = parser.parse_args()

    files = sorted(args.input_dir.glob("*.nii*"))
    print(f"Found {len(files)} volumes")

    slices = []
    for f in files:
        try:
            slices.append(render_midslice(f))
        except Exception as e:
            print(f"  Skipping {f.name}: {e}")

    grid = make_grid(slices, cols=args.cols)

    # Save as PNG using raw bytes (avoid matplotlib/PIL dependency)
    import struct
    import zlib

    def write_png(path: Path, data: np.ndarray):
        h, w = data.shape
        raw = b""
        for row in range(h):
            raw += b"\x00" + data[row].tobytes()

        def chunk(ctype, cdata):
            c = ctype + cdata
            return struct.pack(">I", len(cdata)) + c + struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)

        sig = b"\x89PNG\r\n\x1a\n"
        ihdr = struct.pack(">IIBBBBB", w, h, 8, 0, 0, 0, 0)
        with open(path, "wb") as f:
            f.write(sig)
            f.write(chunk(b"IHDR", ihdr))
            f.write(chunk(b"IDAT", zlib.compress(raw)))
            f.write(chunk(b"IEND", b""))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_png(args.output, grid)
    print(f"Saved QC grid: {args.output} ({grid.shape[1]}x{grid.shape[0]})")
