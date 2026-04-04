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

| Dataset    | Est. volumes | Modalities | Status                                    |
| ---------- | ------------ | ---------- | ----------------------------------------- |
| ADNI full  | ~10,913      | T1         | Downloading on Mac (250GB), LONI          |
| ABIDE I    | ~2,228       | T1         | Downloaded on Sherlock (S3)               |
| ABIDE II   | ~1,100       | T1         | Downloaded on Sherlock (S3)               |
| BraTS 2023 | ~5,000       | T1,T1c,T2,FLAIR | Downloading on Sherlock (Synapse)    |
| HCP-YA     | ~2,226       | T1, T2     | Downloading on Sherlock (S3)              |
| HCP-Aging  | ~2,789       | T1, T2     | Access granted, S3 bucket TBD            |

## Pending access

| Dataset   | Est. volumes | Status                              |
| --------- | ------------ | ----------------------------------- |
| OASIS-3   | ~2,842       | DUA approved, NITRC access pending  |
| OASIS-4   | ~600+        | Access requested                    |
| ABCD      | ~12,000      | Applying via NDA (pediatric, ages 9-10) |
| HCP Aging | ~5,578       | Access granted, need external SSD for 2TB Aspera download |
| NKI-RS    | ~1,000       | On S3 (fcp-indi), no auth, not yet checked |

## Not usable

- `ADNI/processed/`, `SCHIZO/processed/norm*`, `TCGA/processed/` — MNI152 atlas, wrong space.

## Pipeline fixes — ALL RESOLVED

- OASIS-1: T88 → SRI24 re-registration
- OASIS-2: Raw axis permutation + SRI24 registration
- SCHIZO: SS atlas + CoM + skip HD-BET
- Multi-modal: T1 as center, others as moving

## Download commands reference

### ABIDE (S3, no auth)
```bash
module load aws-cli
aws s3 sync s3://fcp-indi/data/Projects/ABIDE/RawDataBIDS/ data/ABIDE/ --no-sign-request --exclude "*" --include "*T1w*"
aws s3 sync s3://fcp-indi/data/Projects/ABIDE2/RawData/ data/ABIDE2/ --no-sign-request --exclude "*" --include "*T1w*"
```

### HCP Young Adult (S3, needs credentials)
```bash
module load aws-cli
export AWS_ACCESS_KEY_ID="<from connectomedb>"
export AWS_SECRET_ACCESS_KEY="<from connectomedb>"

# Skull-stripped + bias-corrected (~38GB for 1113 subjects)
aws s3 sync s3://hcp-openaccess/HCP_1200/ data/HCP_YA/ --region us-east-1 \
  --exclude "*" \
  --include "*/T1w/T1w_acpc_dc_restore_brain.nii.gz" \
  --include "*/T1w/T2w_acpc_dc_restore_brain.nii.gz"

# With skull (backup, ~190GB) — add if needed:
#  --include "*/T1w/T1w_acpc_dc_restore.nii.gz"
#  --include "*/T1w/T2w_acpc_dc_restore.nii.gz"
```
HCP files: `T1w/` folder per subject. `_restore_brain` = bias-corrected + skull-stripped (FreeSurfer).
`_restore` = bias-corrected, with skull. Both in native ACPC space, need SRI24 registration.
Pre-stripped → use same approach as SCHIZO (SS atlas + CoM, skip HD-BET).

### BraTS 2023 (Synapse, needs token)
```bash
pip install synapseclient
synapse login -p "<your_synapse_personal_access_token>"
synapse get -r syn51156910 --downloadLocation data/BraTS2023/
```
BraTS data is already preprocessed (SRI24, skull-stripped, multi-modal). Likely overlaps with TCGA/UPENN.

### ADNI full (LONI, IP-bound)
Download via LONI IDA web interface to local Mac, then rsync to Sherlock.
Collection: ADSP-PHC ADNI T1 1.0 (10,913 scans, 2,592 subjects, 8,602 NIfTI + 2,311 DICOM).
```bash
caffeinate -i curl -L -C - -o ADNI_T1.zip "https://ida.loni.usc.edu/download/files/ida1/<session-id>/ADSP-PHC%3A%20ADNI%20T1%201.0.zip"
# Then rsync to Sherlock
rsync -ahP ADNI_T1.zip maxvpuyv@login.sherlock.stanford.edu:/oak/.../data/ADNI_full/
```

### HCP Aging / AABC (Aspera, IP-bound token)
ConnectomeDB → AABC Release 2 → Structural Preprocessed (2,788 subjects, ~2TB).
Aspera only, no S3 access. Token is IP-bound — cannot transfer from Sherlock.
```bash
# Server details (from Aspera Connect logs):
#   host: asp-connect1.wustl.edu
#   user: asperaxfer
#   port: 33001
#   source: packages/prerelease/aabc/AABC_FZ1/
#   auth: ASPERA_SCP_TOKEN (session-based, IP-bound)
#   key: aspera_tokenauth_id_rsa
#
# To download on Sherlock, would need:
#   1. Aspera token generated from same IP as Sherlock (not possible via browser)
#   2. Or S3 access (not yet available for AABC)
#   3. Or download 2TB to Mac via Aspera Connect, then rsync (impractical)
#
# TODO: Email hcp-users@humanconnectome.org asking for S3 access to AABC structural data
# Each subject is a zip: HCA{id}_V{visit}_MR_StructuralRecommended.zip
# Contains T1w + T2w + FreeSurfer outputs. We only need T1w/T2w_restore_brain files.
```

## Pipeline fixes — summary

- OASIS-1: T88 → SRI24 re-registration
- OASIS-2: Raw axis permutation + SRI24 registration
- SCHIZO: SS atlas + CoM + skip HD-BET (`pre_skull_stripped`)
- ABIDE: SS atlas + CoM for defaced inputs (`use_ss_atlas` + `com_align`)
- HCP: SS atlas + CoM, pre-stripped (`pre_skull_stripped`)
- Multi-modal: T1 as center, others as moving

## Totals

| Status          | Volumes      |
| --------------- | ------------ |
| Preprocessed    | ~12,511      |
| In progress     | ~15,668      |
| Pending access  | ~22,020+     |
| **Projected**   | **~50,199+** |
