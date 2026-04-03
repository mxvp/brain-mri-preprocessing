# Dataset Status

Target: SRI24 (240x240x155, 1mm iso, LPS, skull-stripped, raw intensity)

## Preprocessing done

| Dataset  | Volumes | Modalities          | Notes                           |
| -------- | ------- | ------------------- | ------------------------------- |
| UPENN    | 2,684   | T1, T1GD, T2, FLAIR | Already in SRI24                |
| TCGA     | 668     | T1, T1Gd, T2, FLAIR | Already in SRI24                |
| UCSF     | 2,004   | T1, T1c, T2, FLAIR  | Already in SRI24, bias-corrected|
| ADNI     | 815     | T1                   | Full pipeline                   |
| Stanford | 158     | T1Gd, FLAIR          | Full pipeline                   |
| SCHIZO   | 670     | T1, T2               | Reg only (pre-stripped)         |
| OASIS-1  | 436     | T1                   | T88 → SRI24                    |
| OASIS-2  | 373     | T1                   | Raw permute → SRI24            |
| IXI      | ~795    | T1, T2               | Full pipeline, done             |
| PPMI     | ~3,908  | T1, T2, FLAIR        | Full pipeline, done             |

**Subtotal preprocessed: ~12,511 volumes**

## Downloading / incoming

| Dataset   | Est. volumes | Modalities | Status                              |
| --------- | ------------ | ---------- | ----------------------------------- |
| ADNI full | ~10,913      | T1         | Downloading on Mac (250GB), NIfTI+DICOM via LONI |
| ABIDE I   | ~1,100       | T1         | Downloading on Sherlock (S3, no auth) |
| ABIDE II  | ~1,100       | T1         | Downloading on Sherlock (S3, no auth) |

## Pending access

| Dataset   | Est. volumes | Status                              |
| --------- | ------------ | ----------------------------------- |
| OASIS-3   | ~2,842       | DUA approved, NITRC access pending (emailed) |
| OASIS-4   | ~600+        | Access requested                    |
| HCP       | ~1,800       | Not yet applied                     |
| BraTS 2023| ~5,000       | Available on Synapse, not started   |
| ADNI (full)| ~10K more   | Access granted, downloading         |

## Not usable

- `ADNI/processed/`, `SCHIZO/processed/norm*`, `TCGA/processed/` — MNI152 atlas, wrong space.

## Pipeline fixes — ALL RESOLVED

- OASIS-1: T88 → SRI24 re-registration
- OASIS-2: Raw axis permutation + SRI24 registration
- SCHIZO: SS atlas + CoM + skip HD-BET
- Multi-modal: T1 as center, others as moving

## Totals

| Status          | Volumes     |
| --------------- | ----------- |
| Preprocessed    | ~12,511     |
| Downloading     | ~13,113     |
| Pending access  | ~10,242+    |
| **Projected**   | **~35,866+**|
