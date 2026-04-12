# Dataset Status

Target: SRI24 (240x240x155, 1mm iso, LPS, skull-stripped, raw intensity)

## Preprocessing done

| Dataset  | Volumes | Modalities          | Notes                           |
| -------- | ------- | ------------------- | ------------------------------- |
| UPENN    | 2,684   | T1, T1GD, T2, FLAIR | Already in SRI24                |
| TCGA     | 668     | T1, T1Gd, T2, FLAIR | Already in SRI24                |
| UCSF     | 2,004   | T1, T1c, T2, FLAIR  | Already in SRI24, bias-corrected|
| BraTS    | ~5,000  | T1, T1c, T2, FLAIR  | Already in SRI24                |
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
| ADNI     | 10,912  | T1                   | Full pipeline (GPU)             |
| CORR     | 546     | T1                   | fmriprep MNI → SRI24            |
| FCON1000 | 1,250   | T1                   | Pre-stripped, native → SRI24    |
| ADHD200  | 599     | T1                   | fmriprep MNI → SRI24            |
| BGSP     | 1,636   | T1                   | Skull-stripped 1.2mm → SRI24    |
| NKI2     | 222     | T1                   | MNI → SRI24                     |
| HBN      | 2,136   | T1                   | qsiprep MNI → SRI24             |

**Subtotal done: ~45,048 volumes across 200+ sites worldwide** (ADNI ~60, FCON1000 ~35, PPMI ~33, CORR ~30, ABIDE I/II 36, ADHD200 8, HBN 4, IXI 3, BGSP 2, BraTS/TCGA pooled multi-institution, rest single-site)

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
| Healthy | IXI, HCP-YA, NKI, NKI2, CORR, FCON1000, BGSP, OASIS controls, ADNI CN, ABIDE controls, PPMI controls, ADHD200 controls |
| Alzheimer's/dementia | ADNI, OASIS-1/2 |
| Brain tumors (GBM/glioma) | UPENN, TCGA, UCSF, Stanford, BraTS |
| Parkinson's | PPMI |
| Autism | ABIDE I/II |
| Schizophrenia | SCHIZO |
| ADHD | ADHD200 |
| Meningioma | BraTS-MEN |
| Pediatric tumors | BraTS-PED |
| Pediatric psychiatric (mixed) | HBN (ADHD, ASD, anxiety, learning disorders) |

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
| NKI / NKI2 | CC BY-NC | No |
| BGSP | BGSP Data Use Terms (non-commercial) | No |
| Stanford | Internal (in-house) | N/A |
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

## Metadata / phenotypes

Raw clinical sources live in `data/metadata/<DATASET>/` (gitignored). Code that builds harmonized masters from them is in `metadata.py`; code that joins masters to preprocessed files on disk is in `build_volumes.py`.

**Pipeline**:
1. `python metadata.py <dataset> --input data/metadata/<DATASET>/ -o <dataset>_master.csv` — per-dataset harmonized master, 25-col SCHEMA (subject, scan, visit, age, sex, dx, mmse, cdr, moca, mds_updrs_iii, gds, scopa_aut, apoe, scanner, modality, ...).
2. `python metadata.py concat --input data/metadata -o all_datasets_master.csv` — concatenate all masters.
3. On Sherlock: `find data/preprocessed GBM_MAE/data/UPENN/normalize UCSF-PDGM-v3 -name '*.nii.gz' > files.txt`
4. `python build_volumes.py --files files.txt --metadata-root data/metadata -o volumes.csv` — **one row per preprocessed .nii.gz**, joined with clinical.

**`volumes.csv` is the training-ready file.** Columns: `file_path, dataset, subject_id, session_id, modality, scan_date, age_at_scan, sex, dx, dx_detail, site, mmse, cdr_global, moca, mds_updrs_iii, gds, scopa_aut, apoe, ...`.

### Current volumes.csv state
- **34,453 volume rows × 17,231 subjects × 20 datasets**
- **Clinical dx match: 96.5%** (33,246 / 34,453)
- Pretraining: use all 34,453 rows (every row has `file_path`). ControlNet: filter `df[df.dx.notna()]`.

### Per-dataset coverage

| Dataset | Rows | With-dx | Notes |
|---------|-----:|--------:|-------|
| ADNI | 10,910 | 100% | XML per-scan (CN/EMCI/LMCI/AD) + LONI longitudinal (MMSE/CDR/MoCA) |
| PPMI | 3,908 | 100% | Per-subject baseline (cohort, UPDRS-III, MoCA, GDS, SCOPA-AUT, APOE + PD variants) |
| NKI | 2,455 | 64% | 876 orphans need COINS DUA |
| HCP-YA | 2,226 | 100% | BALSA Column Selector export (age bin, sex, MMSE) |
| HBN | 2,135 | 99% | 4 per-site participants.tsv; rich instruments need HBN DUA |
| UCSF | 1,980 | 100% | 4 structural × 495 subjects (bias-corrected only) |
| BGSP | 1,636 | HC only | Dataverse CSV still pending |
| ABIDE-II | 1,430 | 92% | Composite + 22 per-site |
| FCON1000 | 1,197 | 97% | 33 per-site participants.tsv |
| IXI | 1,159 | 97% | xls |
| ABIDE-I | 984 | 100% | Composite + 19 per-site |
| ADNI old | 815 | 100% | Pre-full download batch |
| SCHIZO | 670 | 100% | COBRE+MCIC merged |
| UPENN | 615 | 97% | cBioPortal + local clinical |
| ADHD200 | 598 | 96% | Master TSV + 9 per-site |
| CoRR | 546 | 83% | Aggregated phenotypic |
| OASIS-1 | 436 | 100% | CDR + MMSE |
| OASIS-2 | 373 | 100% | Longitudinal CDR |
| NKI2 | 222 | 100% | Backfilled from NKI1 |
| Stanford | 158 | 100% | Brain tumor cohort |

### Gaps (not fixable without more downloads)
- **NKI full phenotype** — 388 subjects need COINS DUA
- **BGSP demographics** — waiting on Dataverse request
- **HBN rich instruments** — need HBN DUA (WISC, CBCL, KSADS)
- **BraTS 2023 clinical** — only Synapse image manifest, no clinical pulled

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
| Done               | ~45,048      |
| Pending            | ~118,020+    |
| **Projected**      | **~163,068+**|
