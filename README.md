# brain-mri-preprocessing

Standardized preprocessing pipeline for brain MRI: N4 bias correction, skull-stripping, and affine registration to the SRI24 atlas.

Supports per-subject multi-modality preprocessing where the center modality (T1 or T1c) drives atlas registration and other modalities (T2, FLAIR) are co-registered through it.

A parallel prostate MRI pipeline lives in `prostate/` — see `prostate/README.md`.

## Repo structure
```
datasets.py            # Dataset registry — per-dataset file discovery, conversion, subject grouping
prepare.py             # Find and prepare files, output manifest.json
preprocess.py          # Per-subject preprocessing (wraps brainles-preprocessing)
_cgga_parallel.py      # Shard-based parallel driver for preprocess.py (splits a manifest across N processes)
convert.py             # Standalone format converter (Analyze, MGZ, MINC, DICOM)
diagnosis.py           # Analyze NIfTI files: format, atlas, modality, skull-strip status
qc.py                  # Render mid-slice PNG grid for visual QC
build_volumes.py       # Legacy helpers for older cohort dumps
helpers.py             # Shared helpers
metadata.py            # Clinical / demographic metadata parsing (per-cohort)

slurm_scripts/         # SLURM job scripts for HPC
  preprocess.sh        # Template: one dataset per job (GPU)
  preprocess_cpu.sh    # CPU-only variant (for pre-skull-stripped cohorts, HD-BET skipped)
  run_all.sh           # Submit all datasets

prostate/              # Prostate MRI pipeline — separate CLI `python -m prostate`

STATUS.md              # Per-dataset preprocessing status tracker
pyproject.toml         # Dependencies (uv-managed)
```

## Supported datasets

### Registered (`REGISTRY` in `datasets.py`)

Datasets you can invoke via `python prepare.py <name> <input> <output>`:

| Name | Cohorts | Format | Modalities | Notes |
|---|---|---|---|---|
| `ixi` | IXI | .nii.gz | T1, T2 | Anisotropic; 2D T2 sequences filtered |
| `ppmi` | PPMI | DICOM | T1, T2, FLAIR | 2D sequences filtered, dcm2niix conversion |
| `adni` | ADNI | .nii | T1 | |
| `schizo` | COBRE | .nii.gz | T1, T2 | Pre-skull-stripped: SS atlas + CoM align, skips HD-BET |
| `cgga` | Chinese Glioma Genome Atlas | .nii.gz | CE (T1c) | Pre-BET-stripped; raw clinical (LAS, ~6 mm slices); CE-only manifest (T1/T2 dropped for the current downstream) |
| `stanford` | Stanford tumor | .nii.gz | T1Gd, FLAIR | Negative intensities clipped to 0 |
| `bids_defaced` | ABIDE I/II, NKI, CORR, FCON1000 | BIDS .nii.gz | T1 | Defaced |
| `fmriprep_mni` | ADHD200, CORR (fmriprep outputs) | .nii.gz | T1 | Already MNI-registered; re-registered to SRI24 |
| `bgsp` | Brain Genomics Superstruct Project | .nii.gz | T1 | |
| `hcp` | HCP-Young Adult | .nii.gz | T1 | |
| `upenn` | UPenn-GBM | BIDS .nii.gz | T1, T1c, T2, FLAIR | Baseline (`_11`) sessions only; post-op (`_21`) excluded |
| `ucsf` | UCSF-PDGM | .nii.gz | T1, T1c, T2, FLAIR | |

**Not currently registered** — classes exist in `datasets.py` but commented out:
- `OASIS1`, `OASIS2`, `OASIS3` — legacy, need re-audit
- `TCGA` — needs the same audit as UPENN/UCSF got before re-enabling
- `FCON1000` — covered by `bids_defaced`

### Preprocessed elsewhere (feed straight into an encoder)

Cohorts that arrive already-SRI24-preprocessed from collaborators — no registry entry needed; drop the files into `data/preprocessed/<name>/` and encode directly:

| Cohort | Volumes | Provenance |
|---|---:|---|
| Rohan's TCGA / CPTAC / UPENN dumps | varies | Colleague sends already SRI24-aligned, LPS, skull-stripped |
| Ivy GAP (`ivy`) | 31 | nnU-Net `imagesTr/` layout; encoder-ready except `sub-W36` (severely truncated) |

To add a new dataset: subclass `Dataset` in `datasets.py`, implement `prepare()`, add to `REGISTRY`.

## Pipeline

1. **`prepare.py <name> <input_dir> <output_dir>`** — finds files per dataset, converts formats, groups modalities by subject, writes `manifest.json`.
2. **`preprocess.py --manifest <path> --output <dir>`** — reads the manifest, runs per subject: N4 bias correction → skull-strip (HD-BET, GPU) → affine registration to SRI24 → LPS reorient.
   - Center modality (T1 / T1c) drives atlas registration.
   - Moving modalities (T2 / FLAIR / T1c when T1 is center) are co-registered through the center.
   - `pre_registered: true` → skip atlas registration (subject arrived already in SRI24).
   - `pre_skull_stripped: true` → use skull-stripped atlas + CoM alignment, skip HD-BET (SCHIZO, CGGA).
   - `com_align: true` → CoM-align to atlas before registration (implied by `pre_skull_stripped`).
   - Negative input intensities clipped to 0 automatically.
   - **Every output is reoriented to LPS** as a post-step — `brainles_preprocessing` sometimes inherits the input's orientation labels rather than the SRI24 atlas's, so a LAS input (e.g. CGGA) would otherwise write LAS output and be A↔P-flipped vs the encoder's training data.
3. **`qc.py`** — mid-slice PNG grid for visual sanity.

## Parallel preprocessing

For datasets big enough that per-subject serial execution is painful, use `_cgga_parallel.py` (name is legacy, works for any cohort):

```bash
python _cgga_parallel.py <manifest.json> <output_dir> <n_workers>
```

Shards the manifest stride-wise across N subprocesses, each running its own `preprocess.py`. `brainles_preprocessing` uses per-call tempdirs so concurrent runs don't collide. On a 32-core Sherlock node with 12 workers, ~50 s/subject → ~30–45 min for a few hundred subjects (CE-only single-modality).

## Output spec

| Property | Value |
|---|---|
| Shape | 240 × 240 × 155 (SRI24 atlas space) |
| Voxel size | 1 mm isotropic |
| Orientation | LPS |
| Intensity | Raw (not normalized) |
| Non-brain voxels | 0 |
| Format | .nii.gz |

## Usage

```bash
# Step 1: Prepare (find files, convert, group by subject)
python prepare.py ixi data/IXI staging/IXI
python prepare.py ppmi data/PPMI staging/PPMI
python prepare.py cgga data/CGGA staging/CGGA
python prepare.py --list

# Step 2a: Preprocess locally (single process)
python preprocess.py --manifest staging/IXI/manifest.json --output preprocessed/IXI

# Step 2b: Preprocess in parallel (many subprocesses)
python _cgga_parallel.py staging/CGGA/manifest.json data/preprocessed/cgga 12

# Step 2c: Preprocess on HPC (SLURM)
sbatch --job-name=preproc-ixi slurm_scripts/preprocess.sh staging/IXI/manifest.json preprocessed/IXI/

# Or submit all at once
bash slurm_scripts/run_all.sh

# Single volume (quick test)
python preprocess.py input.nii.gz --output output.nii.gz --device cpu

# Diagnose a raw dataset (shape, orient, skull-strip state, modality guess)
python diagnosis.py data/IXI -r

# QC render
python qc.py preprocessed/IXI/ qc_ixi.png
```
