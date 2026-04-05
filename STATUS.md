# Dataset Status

Target: SRI24 (240x240x155, 1mm iso, LPS, skull-stripped, raw intensity)

## Preprocessing done

| Dataset  | Volumes | Modalities          | Notes                           |
| -------- | ------- | ------------------- | ------------------------------- |
| UPENN    | 2,684   | T1, T1GD, T2, FLAIR | Already in SRI24                |
| TCGA     | 668     | T1, T1Gd, T2, FLAIR | Already in SRI24                |
| UCSF     | 2,004   | T1, T1c, T2, FLAIR  | Already in SRI24, bias-corrected|
| BraTS    | ~5,000  | T1, T1c, T2, FLAIR  | Already in SRI24                |
| ADNI old | 815     | T1                   | Full pipeline                   |
| Stanford | 158     | T1Gd, FLAIR          | Full pipeline                   |
| SCHIZO   | 670     | T1, T2               | Reg only (pre-stripped)         |
| OASIS-1  | 436     | T1                   | T88 → SRI24                    |
| OASIS-2  | 373     | T1                   | Raw permute → SRI24            |
| IXI      | ~795    | T1, T2               | Full pipeline                   |
| PPMI     | ~3,908  | T1, T2, FLAIR        | Full pipeline                   |

**Subtotal done: ~17,511 volumes**

## Preprocessing running (weekend batch)

| Dataset    | Est. volumes | GPUs | Status                                    |
| ---------- | ------------ | ---- | ----------------------------------------- |
| ADNI full  | ~10,000      | 8    | Preparing + submitting                    |
| ABIDE I    | 1,102        | 3    | Running (SS atlas + CoM)                  |
| ABIDE II   | 1,427        | 4    | Running (SS atlas + CoM)                  |
| HCP-YA     | 2,226        | 4    | Running (pre-stripped, ~742 remaining)    |
| NKI        | 2,455        | 4    | Running (SS atlas + CoM)                  |

**Subtotal running: ~17,210 volumes**

## Downloaded, needs preprocessing

| Dataset    | Est. volumes | Notes                                     |
| ---------- | ------------ | ----------------------------------------- |
| CORR       | ~1,600       | Downloaded (301GB, healthy adults)         |
| FCON1000   | ~1,000       | Downloaded (13GB, healthy adults)          |
| ADHD200    | ~970         | Downloaded (212GB, pediatric, low priority)|

## Pending access / download

| Dataset   | Est. volumes | Status                              |
| --------- | ------------ | ----------------------------------- |
| OASIS-3   | ~2,842       | DUA approved, NITRC access pending  |
| OASIS-4   | ~600+        | Access requested                    |
| ABCD      | ~12,000      | Need PI sign-off for NDA (pediatric)|
| HCP Aging | ~5,578       | Need external SSD for 2TB Aspera    |
| HBN       | ~2,000       | Killed download (pediatric, 1TB+, low priority) |

## Not usable

- `ADNI/processed/`, `SCHIZO/processed/norm*`, `TCGA/processed/` — MNI152 atlas, wrong space.

## Pipeline fixes — summary

- OASIS-1: T88 → SRI24 re-registration
- OASIS-2: Raw axis permutation + SRI24 registration
- SCHIZO: SS atlas + CoM + skip HD-BET (`pre_skull_stripped`)
- ABIDE/NKI: SS atlas + CoM for defaced inputs (`use_ss_atlas` + `com_align`)
- HCP: SS atlas + CoM, pre-stripped (`pre_skull_stripped`)
- ADNI full: dedup by subject+date, prefer best NIfTI variant, DICOM fallback
- Multi-modal: T1 as center, others as moving

## Totals

| Status             | Volumes      |
| ------------------ | ------------ |
| Done               | ~17,511      |
| Running            | ~17,210      |
| Downloaded         | ~3,570       |
| Pending            | ~23,020+     |
| **Available soon** | **~38,291**  |
| **Projected**      | **~61,311+** |
