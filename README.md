# brain-mri-preprocessing

Standardized preprocessing pipeline for brain MRI: N4 bias correction, skull-stripping, and atlas registration to SRI24.

Supports per-subject multi-modality preprocessing where T1 drives atlas registration and other modalities (T2, FLAIR) are co-registered through T1.

## Repo structure
```
datasets.py            # Dataset registry — per-dataset file discovery, conversion, subject grouping
prepare.py             # Find and prepare files, output manifest.json
preprocess.py          # Per-subject preprocessing (wraps brainles-preprocessing)
diagnosis.py           # Analyze NIfTI files: format, atlas, modality, skull-strip status
slurm_scripts/         # SLURM job scripts for HPC
  preprocess.sh        # Template: one dataset per job
  run_all.sh           # Submit all datasets
convert.py             # Standalone format converter (Analyze, MGZ, MINC, DICOM)
qc.py                  # Render mid-slice PNG grid for visual QC
STATUS.md              # Per-dataset status tracker
pyproject.toml         # Dependencies
```

## Supported datasets

### Needs preprocessing (in registry)

| Dataset | Format | Modalities | Notes |
|---------|--------|-----------|-------|
| IXI | .nii.gz | T1, T2 | Anisotropic, some T2 are 2D (filtered) |
| OASIS-1 | Analyze | T1 | Uses T88 (Talairach) version, re-registered to SRI24 |
| OASIS-2 | Analyze | T1 | Raw axis permutation (AP,SI,LR)->(LR,AP,SI) + SRI24 registration |
| OASIS-3 | BIDS .nii.gz | T1 | |
| PPMI | DICOM | T1, T2, FLAIR | 2D sequences filtered, dcm2niix conversion |
| ADNI | .nii | T1 | |
| SCHIZO | .nii.gz | T1, T2 | Pre-skull-stripped: uses SS atlas + CoM alignment, skips HD-BET |
| Stanford | .nii.gz | T1Gd, FLAIR | Negative values clipped |

### Already in SRI24 (no preprocessing needed)

| Dataset | Modalities | Volumes | Notes |
|---------|-----------|---------|-------|
| UPENN | T1, T1GD, T2, FLAIR | 2,684 | SRI24, skull-stripped |
| TCGA | T1, T1Gd, T2, FLAIR | 668 | SRI24, skull-stripped |
| UCSF | T1, T1c, T2, FLAIR | 2,004 | SRI24, skull-stripped, bias-corrected |

To add a new dataset: subclass `Dataset` in `datasets.py`, implement `prepare()`, add to `REGISTRY`.

## Pipeline

1. **prepare.py** — finds files per dataset, converts formats, groups modalities by subject, outputs `manifest.json`
2. **preprocess.py** — reads manifest, runs per-subject: N4 + skull-strip + affine registration to SRI24
   - T1/T1Gd = center modality (drives atlas registration)
   - T2/FLAIR = moving modalities (co-registered through T1)
   - `pre_registered: true` — skip atlas registration (OASIS, already registered during prepare)
   - `pre_skull_stripped: true` — use skull-stripped atlas + CoM alignment, skip HD-BET (SCHIZO)
   - Negative input values clipped to 0 automatically
3. **qc.py** — visual inspection of outputs

## Output spec
| Property | Value |
|----------|-------|
| Shape | 240 x 240 x 155 (SRI24 atlas space) |
| Voxel size | 1mm isotropic |
| Intensity | Raw (not normalized) |
| Non-brain voxels | 0 |
| Format | .nii.gz |

## Usage

```bash
# Step 1: Prepare (find files, convert, group by subject)
python prepare.py ixi data/IXI staging/IXI
python prepare.py ppmi data/PPMI staging/PPMI
python prepare.py oasis1 data/OASIS/OASIS1 staging/OASIS1
python prepare.py --list

# Step 2a: Preprocess locally
python preprocess.py --manifest staging/IXI/manifest.json --output preprocessed/IXI

# Step 2b: Preprocess on HPC (SLURM)
sbatch --job-name=preproc-ixi slurm_scripts/preprocess.sh staging/IXI/manifest.json preprocessed/IXI/

# Or submit all at once
bash slurm_scripts/run_all.sh

# Single volume (quick test)
python preprocess.py input.nii.gz --output output.nii.gz --device cpu

# Diagnose a dataset
python diagnosis.py data/IXI -r

# QC
python qc.py preprocessed/IXI/ qc_ixi.png
```
