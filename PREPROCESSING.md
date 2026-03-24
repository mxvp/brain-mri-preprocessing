# MRI-COLLECT

Brain MRI preprocessing pipeline for GBM-MAE. Takes raw volumes from any dataset and produces
skull-stripped, atlas-registered volumes ready for the GBM-MAE data loader.

## Repo structure
```
preprocess.py          # Single-volume preprocessing (wraps brainles-preprocessing)
preprocess_slurm.sh    # SLURM array job for batch processing on Sherlock
convert.py             # Convert any format to NIfTI (Analyze, MGZ, MINC, DICOM)
qc.py                  # Render mid-slice PNG grid for visual QC
pyproject.toml         # Dependencies
data/                  # Sample volumes for testing
```

## Preprocessing target

Match existing UPenn-GBM / UCSF-PDGM data so volumes are ready for the GBM-MAE data pipeline.

## What the GBM-MAE loader handles at load time
- Resize/pad/crop to 160x192x160
- Z-score normalization (nonzero, channel-wise)
- Augmentations (flips, intensity shifts)

## What must be done BEFORE the loader

### 1. N4 bias field correction
Remove scanner-induced intensity inhomogeneity.

### 2. Skull-stripping
Remove non-brain tissue (skull, scalp, eyes).

### 3. Affine registration to SRI24 atlas
Register to SRI24 atlas space. Result: ~240x240x155, 1mm isotropic voxels, LPS orientation.

### 4. Save as .nii.gz

## Tool: brainles-preprocessing

First choice: `pip install brainles-preprocessing` — does all 3 steps in one pipeline.
- Uses ANTs for N4 + registration, HD-BET or SynthStrip for skull-stripping
- Supports SRI24 atlas out of the box (use `Atlas.SRI24`, not the default `Atlas.BRATS_SRI24`)
- Works on any brain MRI — no lesion-specific logic despite the name
- Validated on healthy brains (HD-BET Dice >96.9% on LPBA40/NFBS/CC-359)
- GitHub: https://github.com/BrainLesion/preprocessing
- Docs: https://brainles-preprocessing.readthedocs.io

Fallback: build our own pipeline wrapping ANTs + HD-BET/SynthStrip directly.

## Output spec
| Property | Value |
|----------|-------|
| Shape | ~240 x 240 x 155 (SRI atlas space) |
| Voxel size | 1mm isotropic |
| Orientation | LPS |
| Intensity | Raw (not normalized) |
| Non-brain voxels | 0 |
| Format | .nii.gz |

## Usage

```bash
# Single volume
python preprocess.py raw.nii.gz output.nii.gz

# Batch (all .nii* in directory)
python preprocess.py input_dir/ output_dir/ --batch

# CPU only (no GPU)
python preprocess.py raw.nii.gz output.nii.gz --device cpu

# Convert non-NIfTI formats first (Analyze, MGZ, MINC, DICOM)
python convert.py input.hdr output.nii.gz
python convert.py dicom_dir/ output.nii.gz
python convert.py input_dir/ output_dir/ --batch

# QC: render mid-slice grid
python qc.py output/ qc_grid.png

# Sherlock batch job
find /path/to/raw/ -name "*.nii.gz" | sort > input_files.txt
sbatch preprocess_slurm.sh input_files.txt /path/to/output/
```

## Notes
- Z-score normalization is NOT needed at preprocessing time. The loader does it.
- Some existing UPenn data has pre-baked z-normalization — this is harmless (double z-score is ~idempotent) but not required.
- Current training data is T1c (post-contrast). New datasets may be T1w (no contrast) — different intensity distribution but same preprocessing pipeline.
- Runtime estimate: ~3-8 min per volume (N4: 1-2 min, skull-strip: 10-30s, registration: 2-5 min).
- Input to preprocess.py must be NIfTI (.nii or .nii.gz). Use convert.py first for Analyze, MGZ, MINC, or DICOM.
- DICOM conversion requires dcm2niix (check: `which dcm2niix` or `module load dcm2niix`).
