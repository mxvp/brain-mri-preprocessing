# prostate

Harmonize four multi-source prostate MRI datasets into a flat, VAE-loadable directory layout with patient-grouped train/val splits.

No preprocessing — the VAE handles resize and normalization at load time. Bias-field / intensity variability is part of the distribution we want the model to learn.

## What this does

1. **Walks each enabled source** with cohort-specific patient ID + modality extraction.
2. **Converts MetaImage (`.mha`, `.mhd`) to NIfTI** via SimpleITK (~0.3 s/file, parallel).
3. **Symlinks already-NIfTI sources** (TCIA-Biopsy, Prostate158) — no disk duplication.
4. **Drops series with < 5 slices** (localizers, derived 1-slice products).
5. **Patient-grouped, source-stratified train/val split** — a patient's volumes (across modalities + views) always land in the same split, never some in each.
6. **Emits `manifest.csv`** — one row per volume with source, patient_id, modality, view, path, src_format, n_slices, split.

## Sources

| Source | Format on disk | Patient ID convention |
|---|---|---|
| PI-CAI | `.mha` | `<patient>_<study>` (multi-study per patient supported) |
| TCIA Prostate-MRI-US-Biopsy | `.nii.gz` (dcm2niix) | `Prostate-MRI-US-Biopsy-NNNN` (per-modality longest acquisition kept) |
| Prostate158 | `.nii.gz` | `<split>_<NNN>` (split prefix disambiguates train/test ID collisions) |
| PROMISE12 | `.mhd` + `.raw` | `<part>_CaseNN` (part prefix disambiguates across TrainingDataPart{1,2,3}/TestData) |

Adding a new cohort = one walker in `sources.py` registered into `WALKERS`.

## Output layout

```
<output_root>/
  t2_axial/<source>_<patient>.nii.gz         # 2,532 — combined axial T2 pool
  t2_sagittal/<source>_<patient>.nii.gz      # 1,498 — PI-CAI only
  t2_coronal/<source>_<patient>.nii.gz       # 1,497 — PI-CAI only
  adc/<source>_<patient>.nii.gz              # 2,487
  dwi/<source>_<patient>.nii.gz              # 2,425
  manifest.csv
```

## Run

```bash
# what each source contributes — no disk writes
python -m prostate discover

# materialize the harmonized layout + manifest
python -m prostate organize

# smoke / debug: cap each source to N records
python -m prostate organize --limit-per-source 100
```

Tunables in `configs/preprocess.yaml`:
- `paths.output_root`
- `sources.<name>.{enabled, root, min_slices}`
- `split.{val_fraction, seed}`
- `workers` (parallel MetaImage conversions; default 8)

## Wiring into the VAE

Each output subdir is a `type: simple` dataset entry. In `VAE/configs/vae_prostate.yaml`:

```yaml
datasets: [prostate_t2_axial]
dataset_paths:
  prostate_t2_axial:
    type: simple
    img_dir: /oak/.../prostate_organized/t2_axial
```

To train on multi-view T2:

```yaml
datasets: [prostate_t2_axial, prostate_t2_sagittal, prostate_t2_coronal]
dataset_paths:
  prostate_t2_axial:    {type: simple, img_dir: .../t2_axial}
  prostate_t2_sagittal: {type: simple, img_dir: .../t2_sagittal}
  prostate_t2_coronal:  {type: simple, img_dir: .../t2_coronal}
```

For held-out evaluation, filter the manifest by `split == 'val'` rather than relying on directory listings.

## Files

- `sources.py` — per-source walkers + `WALKERS` registry
- `convert.py` — SimpleITK MetaImage → NIfTI (with symlink fast-path for already-NIfTI inputs)
- `splits.py` — patient-grouped, source-stratified split + leakage verifier
- `organize.py` — orchestrator
- `cli.py` + `__main__.py` — `python -m prostate <discover|organize>`
- `configs/preprocess.yaml` — paths, sources, split config
- `tests/` — unit tests (11 currently)
