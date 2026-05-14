# multimodal

Pair brain MRI cohorts with bulk RNA-seq for downstream expression-prediction modeling.

Config-driven, idempotent steps. Re-runs don't re-download. Producing a fresh expression matrix is one command: `python -m multimodal all`.

## What it does

| Step | Reads | Writes |
|---|---|---|
| `inventory` | `configs/cohorts.yaml`, local imaging metadata | `data/multimodal/matrices/imaging_subjects.csv` |
| `query`     | GDC REST API | `data/multimodal/manifests/{cohort}_manifest.tsv` (+ `_gdc_full.csv`) |
| `download`  | manifests | `data/multimodal/raw_tsvs/<project>/<file_id>/<file_name>` |
| `matrix`    | raw TSVs | `counts.parquet`, `tpm.parquet`, `gene_meta.parquet`, `sample_meta.parquet` |
| `pairs`     | inventory + sample_meta | `pairs.csv`, `coverage.csv` |

## Cohorts

Defined in `configs/cohorts.yaml`. Currently included:

| Cohort     | Imaging | Public RNA-seq |
|------------|:-:|:-:|
| TCGA-GBM   | ✓ (1110 metadata) | ✓ (GDC `TCGA-GBM`, STAR-Counts, ~287 paired) |
| TCGA-LGG   | ✓ (513) | ✓ (GDC `TCGA-LGG`, ~513 paired) |
| CPTAC-GBM  | ✓ (53 local) | ✓ (GDC `CPTAC-3`, ~40 paired) |
| UPENN-GBM  | ✓ (630) | ✗ (not in any public repo — requires direct ask to Bakas lab) |
| UCSF-PDGM  | ✓ (495) | ✗ (no public bulk RNA-seq) |

## Run

```bash
# From repo root
python -m multimodal all

# Or step-by-step:
python -m multimodal inventory
python -m multimodal query
python -m multimodal download
python -m multimodal matrix
python -m multimodal pairs
```

Tunables live in `configs/curate.yaml`:
- `samples`: `primary_only` / `all_tumor` / `all`
- `output.values`: `counts_and_tpm` / `tpm_only` / `counts_only`
- `gene_filter.gene_types`: which biotypes to keep
- `gene_filter.min_subjects_expressed`, `min_tpm`: per-gene expression filter

## Loading into a model

Single function for downstream training:

```python
from multimodal.load import load_paired_dataset

ds = load_paired_dataset(
    "data/latents/your_encoder.pt",   # from GBM_MAE/scripts/encode_dataset.py
    matrices_dir="data/multimodal/matrices",
)

ds.X            # (N, D)  imaging latents
ds.y            # (N, G)  log1p(TPM) expression
ds.subject_ids  # (N,)    canonical submitter IDs
ds.cohort       # (N,)    TCGA-GBM / TCGA-LGG / CPTAC-3
ds.split        # (N,)    'train' / 'val' / 'test'
ds.gene_ids     # (G,)    versioned Ensembl IDs
ds.gene_names   # (G,)    HGNC symbols
ds.modality     # which imaging modality (default t1c)
```

Defaults that are baked in (kwargs override any of them):

- one row per imaging **subject** (not sample, not modality)
- modality `t1c`, fallback `t1`
- sample type `Primary Tumor`
- multi-sample subjects: TPM averaged
- genes: protein-coding only, expressed (TPM ≥ 1) in ≥ 10 subjects (~17K genes)
- target: `log1p(TPM)`
- QC: drops samples with `protein_coding_fraction < 0.7`
- split: subject-level, cohort-stratified, fixed seed 0

The loader **assumes the latents file contains TCGA / CPTAC paths**. If it doesn't (e.g. you only have UPENN/UCSF latents), it returns an empty dataset with a warning — encode the tumor cohorts first via `scripts/encode_dataset.py` in GBM_MAE pointed at the preprocessed image dirs.

## Output format

`counts.parquet` and `tpm.parquet` are sample × gene matrices (rows are samples, columns are dotless Ensembl IDs). `sample_meta.parquet` carries cohort, subject, sample type, and QC fields per row. `gene_meta.parquet` carries the gene annotation. `pairs.csv` is the join table — one row per (imaging_id, sample_id).

A subject can contribute >1 sample (e.g. primary + recurrent). `pairs.csv` reflects that; the modeling code decides whether to average, pick one, or keep them separate.

## Adding a cohort

1. Add an entry under `cohorts:` in `configs/cohorts.yaml`.
2. If the imaging parser isn't already in `multimodal/cohorts.py`, add one and register it in `PARSERS`.
3. Run `python -m multimodal inventory query` to confirm IDs match the molecular repo.

