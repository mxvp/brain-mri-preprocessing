# Design notes — multimodal

This doc captures the science decisions behind the curation pipeline so future
contributors (and future-you) can challenge or change them without spelunking
the code.

## What problem we're solving

Predict bulk tumor RNA expression from preoperative MRI. Tactically: produce
a clean `subject × gene` matrix paired one-to-one (or one-to-many) with the
existing imaging cohorts so a regression head can be trained on encoder
latents → expression.

## Why these cohorts

Among the cohorts in `brain-mri-preprocessing`, only three have **public** bulk
RNA-seq paired with the imaging subject IDs we already process:

- **TCGA-GBM / TCGA-LGG** — bulk RNA-seq via GDC, harmonized STAR-Counts.
- **CPTAC-GBM** — proteogenomic CPTAC-3 release on GDC.

UPENN-GBM is the largest single-site high-quality MRI cohort (630 subjects)
but its TCIA release contains only IDH1+MGMT marker calls — no bulk RNA-seq.
Closing this gap requires a direct ask to the Bakas/Davatzikos lab.

UCSF-PDGM has only IDH/MGMT/1p19q in its public release.

## Why STAR-Counts gene expression specifically

GDC harmonizes RNA-seq through the [STAR-Counts pipeline](https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/).
Output per sample is one TSV with raw counts plus pre-computed TPM and FPKM.

We ignore:

- **Aligned Reads (BAM)** — gigabytes per sample, redundant given count output.
- **Splice Junction Quantification** — useful for splice-aware modeling but
  out of scope for the first-pass expression head.
- **Transcript Fusion** (Arriba, STAR-Fusion) — categorical events, modeled
  separately if needed.

The legacy TCGA microarray (Affymetrix U133+2) covers ~310 additional TCGA-GBM
subjects but mixes platforms in a way that's painful to harmonize. We exclude
it for the first pass; add later if needed.

## Sample-level handling

A subject can have a primary tumor sample, recurrent samples, and (rarely) a
normal tissue sample. The default policy is `all_tumor`: keep primary +
recurrent. The pairing table reports cohort, sample type, and subject for
each row, so the modeling code can:

- average replicates (deterministic, simple),
- pick the primary (matches "preop MRI ↔ baseline tumor"),
- or keep recurrent samples as a separate timepoint for longitudinal work.

Storing per-sample rows preserves these options; `primary_only` is configurable
if a stricter one-to-one join is needed.

## Why store both counts and TPM

Counts are the raw signal — needed for proper statistical models (DESeq2,
edgeR) and for re-normalization across cohorts. TPM is convenient for fast
loading and most ML use cases. They share an index; cost of storing both is
negligible (~250 MB).

`log1p` is applied at load time rather than baked in, so the same matrix
serves both linear-space (statistics) and log-space (neural-net) consumers.

## Gene filter

Default published matrix: protein-coding only, ≥10 subjects with TPM ≥ 1.
That's ~17–19 K genes, which is what radiogenomic papers typically use as
a starting point for a regression head.

The raw matrix (unfiltered) is preserved on disk so a researcher can re-filter
without re-running the pipeline.

## Reproducibility

- Cohort identity rules live in `configs/cohorts.yaml`, never in code.
- All science knobs live in `configs/curate.yaml`.
- Each step writes its own artifact; re-running a downstream step doesn't
  re-trigger upstream work.
- The MD5 check on download makes re-runs free past the initial pull.
- Output schema is documented in `README.md` and stable across reruns
  (parquet column order, dtypes).

## Things I left out on purpose

- Batch correction across cohorts (ComBat, mutual nearest neighbors). Decision
  punted to the modeling stage — different consumers want different
  treatments.
- Variance-stabilizing transforms beyond `log1p`. Same reason.
- Single-cell / spatial RNA. Out of scope for the bulk-expression head.
- Methylation, CNV, mutations. Same pipeline pattern would work; not yet
  implemented.
