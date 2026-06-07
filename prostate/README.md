# prostate

Per-patient preprocessing for multiparametric prostate MRI, following the PI-CAI / Radboud nnU-Net community baseline.

Pipeline (per patient, in order):

1. **N4 bias correction** on T2w (and ADC)
2. **Prostate gland segmentation** on T2w (stub box for now → swap in Radboud nnU-Net later)
3. **Rigid registration** of ADC + DWI to T2w (SimpleITK Mattes MI)
4. **Crop** around the gland bbox + margin
5. **Resample** all volumes to a common voxel spacing (PI-CAI default: 0.5 × 0.5 × 3.0 mm)
6. **Per-modality normalization** (T2/DWI: z-score over nonzero voxels; ADC: clip 0.5–99.5 percentile then z-score)
7. **Stack as 4 channels**: `[T2w, DWI, ADC, mask]` → one NIfTI per patient

Tuned via `configs/preprocess.yaml`. All steps idempotent — re-runs skip existing outputs.

## Run

```bash
# What's in our dcm2niix output? Counts T2/ADC/DWI per patient.
python -m prostate inventory

# Full pipeline for every complete patient. Use --limit N for smoke tests.
python -m prostate preprocess --limit 5
python -m prostate preprocess
```

Outputs land in `data/prostate/preprocessed/<patient_id>.nii.gz`.

## Stub vs real prostate segmentation

`configs/preprocess.yaml::segmentation.backend` is `stub` by default — uses a 50%-of-FOV center box as the prostate mask. **Do not train on this.** It exists so the pipeline runs end-to-end without the 1.5 GB Radboud nnU-Net checkpoint.

To switch to the real segmenter:
1. Grab the checkpoint from <https://github.com/DIAGNijmegen/AbdomenMRUS-prostate-segmentation>
2. Set `segmentation.backend: nnunet` + `segmentation.nnunet_checkpoint: <path>`
3. Implement the actual nnU-Net call inside `segment._nnunet_mask` (interface documented in the docstring)

The PI-CAI paper reports Dice ≈ 0.96 internal / 0.82 external for that model. TotalSegmentator's prostate channel is Dice ≈ 0.15 (under-segments badly) — don't use it.

## Adding a cohort

Adding e.g. PI-CAI / PROSTATEx:
1. Run dcm2niix into a flat output dir
2. Point `paths.nifti_root` at it
3. If the SeriesDescription patterns differ, extend `series:` in `preprocess.yaml`

## Files

- `patients.py` — discover dcm2niix output, classify by `SeriesDescription`, build per-patient records
- `segment.py` — prostate gland mask (stub / nnU-Net)
- `registration.py` — rigid SimpleITK DWI/ADC → T2w
- `normalize.py` — per-modality intensity normalization
- `preprocess.py` — the per-patient pipeline orchestrator
- `cli.py` + `__main__.py` — `python -m prostate <command>`
- `configs/preprocess.yaml` — all science knobs
- `tests/` — unit tests
