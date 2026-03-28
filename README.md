# brain-mri-preprocessing

Standardized preprocessing pipeline for brain MRI: N4 bias correction, skull-stripping, and atlas registration.

Takes raw volumes from any dataset and produces skull-stripped, SRI24-registered volumes.

## Repo structure
```
datasets.py            # Dataset registry — per-dataset T1 discovery + format conversion
prepare.py             # Find and prepare T1w files from a dataset
preprocess.py          # Single-volume preprocessing (wraps brainles-preprocessing)
preprocess_slurm.sh    # SLURM array job for batch processing on HPC
convert.py             # Standalone format converter (Analyze, MGZ, MINC, DICOM)
qc.py                  # Render mid-slice PNG grid for visual QC
pyproject.toml         # Dependencies
```

## Supported datasets

| Dataset | Format | Selection |
|---------|--------|-----------|
| IXI | .nii.gz | `*-T1.nii.gz` |
| OASIS-1 | Analyze (.hdr/.img) | `RAW/mpr-1` (converted to NIfTI) |
| OASIS-2 | Analyze (.hdr/.img) | `RAW/mpr-1` (converted to NIfTI) |
| OASIS-3 | BIDS .nii.gz | `*T1w.nii.gz` |
| PPMI | DICOM | `*T1*weighted*` / `*MPRAGE*` dirs (converted via dcm2niix) |
| ADNI | .nii | `T1.nii` per subject |
| SCHIZO | .nii.gz | `T1.nii.gz` per subject (pre-skull-stripped, raw unavailable) |
| Stanford | .nii.gz | `T1Gd.nii.gz`, `FLAIR.nii.gz` per subject |
| TCGA | .nii.gz | `*_t1.nii.gz`, `*_t1Gd.nii.gz`, `*_t2.nii.gz`, `*_flair.nii.gz` |

To add a new dataset: subclass `Dataset` in `datasets.py`, implement `prepare()`, add to `REGISTRY`.

## Preprocessing pipeline

1. **N4 bias field correction** — remove scanner intensity inhomogeneity
2. **Skull-stripping** — remove non-brain tissue (HD-BET)
3. **Affine registration to SRI24** — standard atlas space

Powered by [brainles-preprocessing](https://github.com/BrainLesion/preprocessing).

## Output spec
| Property | Value |
|----------|-------|
| Shape | ~240 x 240 x 155 (SRI24 atlas space) |
| Voxel size | 1mm isotropic |
| Orientation | LPS |
| Intensity | Raw (not normalized) |
| Non-brain voxels | 0 |
| Format | .nii.gz |

## Usage

```bash
# Step 1: Prepare — find T1 files, convert formats if needed
python prepare.py ixi data/IXI staging/IXI
python prepare.py ppmi data/PPMI/PPMI staging/PPMI
python prepare.py oasis1 data/OASIS1 staging/OASIS1
python prepare.py --list  # show all supported datasets

# Step 2a: Preprocess locally
python preprocess.py --filelist staging/IXI/files.txt --output preprocessed/IXI
python preprocess.py --filelist staging/IXI/files.txt --output preprocessed/IXI --device cpu

# Step 2b: Preprocess on HPC (SLURM)
N=$(wc -l < staging/IXI/files.txt)
sbatch --array=1-${N}%20 preprocess_slurm.sh staging/IXI/files.txt preprocessed/IXI/

# Single volume (no prepare step needed)
python preprocess.py input.nii.gz output.nii.gz

# QC: visual inspection grid
python qc.py preprocessed/IXI/ qc_ixi.png
```

## Notes
- Z-score normalization is NOT done at preprocessing time — the downstream loader handles it.
- DICOM conversion (PPMI) requires dcm2niix (`module load dcm2niix` on HPC).
- Runtime: ~1-2 min/volume with GPU, ~3-5 min/volume on CPU.
