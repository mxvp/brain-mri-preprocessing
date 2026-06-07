# prostate

Organize multiparametric prostate MRI for VAE training. No preprocessing.

For VAE generative pretraining the raw dcm2niix output is enough — the VAE loader handles resize and normalization at load time, and bias-field / intensity quirks are part of the data variability the model should learn. So this module just:

1. Groups the flat dcm2niix output into per-patient T2 / ADC / DWI records (regex on `SeriesDescription` from each `.json` sidecar)
2. Picks the best acquisition per modality (most slices wins; one T2, one ADC, one DWI per patient)
3. Materializes them as symlinks (or copies) into a flat dir the VAE loader can glob

```
data/prostate/preprocessed/
  Prostate-MRI-US-Biopsy-0001_t2.nii.gz   -> ../../.../nifti/Prostate-MRI-US-Biopsy-0001_..._t2_spc_..._11.nii.gz
  Prostate-MRI-US-Biopsy-0001_adc.nii.gz  -> ../../.../nifti/Prostate-MRI-US-Biopsy-0001_..._ADC_..._5.nii.gz
  Prostate-MRI-US-Biopsy-0001_dwi.nii.gz  -> ../../.../nifti/Prostate-MRI-US-Biopsy-0001_..._CALC_BVAL_..._7.nii.gz
  ...
```

## Run

```bash
# what's in the dcm2niix output? counts T2/ADC/DWI per patient
python -m prostate inventory

# materialize the VAE-loadable flat dir (symlinks; ~5 sec for 732 patients)
python -m prostate organize

# or self-contained copies (~7 GB on disk, slower)
python -m prostate organize  # then edit configs/preprocess.yaml: output.copy: true
```

## Wiring into the VAE

The flat dir is a `type: simple` dataset entry. In `VAE/configs/vae.yaml`:

```yaml
datasets: [..., prostate_t2, prostate_adc, prostate_dwi]
dataset_paths:
  prostate_t2:  {type: simple, img_dir: /oak/.../brain-mri-preprocessing/data/prostate/preprocessed_t2}
  prostate_adc: {type: simple, img_dir: /oak/.../brain-mri-preprocessing/data/prostate/preprocessed_adc}
  prostate_dwi: {type: simple, img_dir: /oak/.../brain-mri-preprocessing/data/prostate/preprocessed_dwi}
```

(If you want per-modality dirs instead of flat-with-prefix, change `output_root` in the config to point at three different dirs and run `organize` three times — or just route by suffix at load time.)

## What this dataset has

The TCIA Prostate-MRI-US-Biopsy collection is a **biopsy-confirmed prostate-cancer cohort**, not healthy prostates. Beyond imaging it has Gleason scores per biopsy core, PI-RADS-like Likert scores per lesion, PSA, and prostate volume measurements (separate `.xlsx` files on TCIA, not included here). Useful later if you want to condition the DiT on cancer grade / risk.

## Files

- `patients.py` — discover dcm2niix output, classify by `SeriesDescription`, pick best per modality
- `organize.py` — materialize the VAE-loadable layout (symlinks or copies)
- `cli.py` + `__main__.py` — `python -m prostate <command>`
- `configs/preprocess.yaml` — paths + series regex + output naming
- `tests/` — unit tests
