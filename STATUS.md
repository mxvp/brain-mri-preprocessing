# Dataset Status

Target: SRI24 (240x240x155, 1mm iso, LPS, skull-stripped, raw intensity)

## Ready (no preprocessing needed)


| Dataset | Subjects | Modalities          | Volumes | Notes                           |
| ------- | -------- | ------------------- | ------- | ------------------------------- |
| UPENN   | 671      | T1, T1GD, T2, FLAIR | 2,684   | SRI24, stripped, LPS            |
| TCGA    | 167      | T1, T1Gd, T2, FLAIR | 668     | SRI24, stripped, LPS            |
| UCSF    | 501      | T1, T1c, T2, FLAIR  | 2,004   | SRI24, stripped, bias-corrected |


## Needs full pipeline (N4 + skull-strip + SRI24 registration)


| Dataset  | Subjects | Modalities    | Volumes | Format       | Issues                               |
| -------- | -------- | ------------- | ------- | ------------ | ------------------------------------ |
| IXI      | 581      | T1, T2        | ~1,159  | NIfTI        | Anisotropic, some T2 are 2D (filter) |
| ADNI     | 818      | T1            | 818     | NIfTI (.nii) | Rotated affine (IPL), anisotropic    |
| Stanford | 80       | T1Gd, FLAIR   | 160     | NIfTI        | Has negatives (clip), some FLAIR 2D  |
| PPMI     | ~1,500   | T1, T2, FLAIR | TBD     | DICOM        | dcm2niix first, many 2D to filter    |
| OASIS-1  | 611      | T1            | 611     | Analyze      | Use T88 version, re-register to SRI24 |
| OASIS-2  | 373      | T1            | 373     | Analyze      | Raw permute (AP,SI,LR)->(LR,AP,SI)   |


## Needs registration only (already skull-stripped)


| Dataset | Subjects | Modalities | Volumes | Issue                                |
| ------- | -------- | ---------- | ------- | ------------------------------------ |
| SCHIZO  | 335      | T1, T2     | 670     | Skip HD-BET, use SS atlas + CoM init |


## Not usable

- `ADNI/processed/`, `SCHIZO/processed/norm*`, `TCGA/processed/` — MNI152 atlas, z-score normalized. Wrong space.

## Pipeline fixes — ALL RESOLVED

- OASIS-1: T88 → SRI24 re-registration (solved)
- OASIS-2: Raw axis permutation + SRI24 registration (solved)
- SCHIZO: SS atlas + CoM + skip HD-BET (solved)
- Multi-modal: T1 as center, others as moving (solved)

## Totals


| Status         | Volumes                   |
| -------------- | ------------------------- |
| Ready          | ~5,356                    |
| Needs pipeline | ~3,494+                   |
| Needs reg only | 670                       |
| **Total**      | **~9,520+** (before PPMI) |


