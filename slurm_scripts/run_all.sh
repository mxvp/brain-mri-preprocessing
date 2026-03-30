#!/bin/bash
# Submit all preprocessing jobs at once
# Run from repo root: bash slurm_scripts/run_all.sh

for ds in ixi oasis1 oasis2 adni stanford tcga schizo upenn; do
  N=$(wc -l < data/staging/$ds/files.txt)
  echo "Submitting $ds: $N volumes"
  sbatch --array=1-${N}%20 preprocess_slurm.sh data/staging/$ds/files.txt data/preprocessed/$ds/
done
