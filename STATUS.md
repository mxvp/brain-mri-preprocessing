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
| IXI      | 1,159   | T1, T2               | Full pipeline                   |
| PPMI     | 3,908   | T1, T2, FLAIR        | Full pipeline                   |
| ABIDE I  | 984     | T1                   | SS atlas + CoM                  |
| ABIDE II | 1,430   | T1                   | SS atlas + CoM                  |
| HCP-YA   | 2,226   | T1, T2               | Pre-stripped → SRI24            |
| NKI      | 2,208   | T1                   | SS atlas + CoM                  |

**Subtotal done: ~24,723 volumes**

## Preprocessing running

| Dataset    | Est. volumes | Compute | Status                                    |
| ---------- | ------------ | ------- | ----------------------------------------- |
| ADNI full  | ~10,000      | GPU     | ~8,141 done, ~2K remaining                |
| CORR       | 546          | CPU     | Running (fmriprep MNI → SRI24)            |
| FCON1000   | 1,250        | CPU     | Running (pre-stripped, native → SRI24)     |
| ADHD200    | 599          | CPU     | Running (fmriprep MNI → SRI24)            |
| BGSP       | 1,636        | CPU     | Running (skull-stripped 1.2mm, 3 jobs)     |
| NKI2       | 222          | CPU     | Running (MNI → SRI24)                     |
| HBN        | 2,136        | CPU     | Running (qsiprep MNI → SRI24, 3 jobs)    |

**Subtotal running: ~20,228 volumes**

## Pending access / download

| Dataset   | Est. volumes | Status                              |
| --------- | ------------ | ----------------------------------- |
| OASIS-3   | ~2,842       | DUA approved, NITRC access pending  |
| OASIS-4   | ~600+        | Access requested                    |
| ABCD      | ~12,000      | Need PI sign-off for NDA (pediatric)|
| HCP Aging | ~5,578       | Need external SSD for 2TB Aspera    |
| HBN (raw) | ~479         | Remaining subjects without qsiprep, needs GPU |
| UK Biobank| ~100,000+    | Need to apply, months for access (100K first scans + 60K repeats) |

## Disease/population coverage

| Category | Datasets |
|----------|----------|
| Healthy | IXI, HCP-YA, NKI, CORR, FCON1000, OASIS controls, ADNI CN, ABIDE controls, PPMI controls |
| Alzheimer's/dementia | ADNI, OASIS-1/2 |
| Brain tumors (GBM/glioma) | UPENN, TCGA, UCSF, Stanford, BraTS |
| Parkinson's | PPMI |
| Autism | ABIDE I/II |
| Schizophrenia | SCHIZO |
| ADHD | ADHD200 |
| Meningioma | BraTS-MEN |
| Pediatric tumors | BraTS-PED |

## Data licenses

| Dataset | License | Commercial OK? |
|---------|---------|----------------|
| IXI | CC BY-SA 3.0 | Yes (with share-alike) |
| UPENN | CC BY 4.0 (TCIA) | Yes |
| TCGA | CC BY 4.0 (TCIA) | Yes |
| UCSF | CC BY 4.0 (TCIA) | Yes |
| BraTS-GLI/MEN/PED | CC BY 4.0 | Yes |
| FCON1000 | CC BY (1000 FCP) | Yes |
| CORR | CC BY (1000 FCP) | Yes |
| HBN | CC BY-NC-SA (older) / CC BY 4.0 (newer) | Partially |
| ABIDE I/II | CC BY-NC-SA | No |
| NKI | CC BY-NC | No |
| ADHD200 | CC BY-NC | No |
| ADNI | Research DUA | No |
| PPMI | Research DUA | No |
| OASIS | Research DUA (non-commercial only) | No |
| HCP | Open Access DUA | No (unclear) |
| SCHIZO | Research DUA | No |
| UK Biobank | Restricted DUA + fees | No (needs commercial license) |

## Not yet explored

| Dataset | Est. subjects | Population | Source |
|---------|---------------|-----------|--------|
| CamCAN | ~650 | Healthy, 18-88yo | camcan-archive.mrc-cbu.cam.ac.uk |
| AOMIC ID1000 | ~1,000 | Healthy | OpenNeuro |
| DLBS | ~315 | Healthy, lifespan | OpenNeuro |
| SALD | ~494 | Healthy adults | NITRC |
| SLIM | ~600 | Healthy, longitudinal | NITRC |
| ENIGMA | 10,000+ | Mixed (consortium) | Via working groups |
| OpenNeuro (bulk) | ~114K (OpenMind) | Mixed, 800 datasets | openneuro.org |

## Pipeline fixes — summary

- OASIS-1: T88 → SRI24 re-registration
- OASIS-2: Raw axis permutation + SRI24 registration
- SCHIZO: SS atlas + CoM + skip HD-BET (`pre_skull_stripped`)
- ABIDE/NKI: SS atlas + CoM for defaced inputs (`use_ss_atlas` + `com_align`)
- HCP/FCON1000: SS atlas + CoM, pre-stripped (`pre_skull_stripped`)
- ADHD200/CORR: fmriprep MNI → SRI24 (`pre_skull_stripped`)
- ADNI full: dedup by subject+date, prefer best NIfTI variant, DICOM fallback
- Multi-modal: T1 as center, others as moving

## Totals

| Status             | Volumes      |
| ------------------ | ------------ |
| Done               | ~24,723      |
| Running            | ~12,395      |
| Pending            | ~118,020+    |
| **Available soon** | **~44,951**  |
| **Projected**      | **~162,971+**|
