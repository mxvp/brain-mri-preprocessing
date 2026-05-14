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

## Output format

`counts.parquet` and `tpm.parquet` are sample × gene matrices (rows are samples, columns are dotless Ensembl IDs). `sample_meta.parquet` carries cohort, subject, sample type, and QC fields per row. `gene_meta.parquet` carries the gene annotation. `pairs.csv` is the join table — one row per (imaging_id, sample_id).

A subject can contribute >1 sample (e.g. primary + recurrent). `pairs.csv` reflects that; the modeling code decides whether to average, pick one, or keep them separate.

## Adding a cohort

1. Add an entry under `cohorts:` in `configs/cohorts.yaml`.
2. If the imaging parser isn't already in `multimodal/cohorts.py`, add one and register it in `PARSERS`.
3. Run `python -m multimodal inventory query` to confirm IDs match the molecular repo.

