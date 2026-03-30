#!/bin/bash
# Submit one preprocessing job per dataset
# Run from repo root: bash slurm_scripts/run_all.sh

mkdir -p logs

for ds in ixi oasis1 oasis2 adni stanford tcga schizo upenn; do
  FILELIST="data/staging/$ds/files.txt"
  if [ ! -f "$FILELIST" ]; then
    echo "SKIP $ds: no files.txt"
    continue
  fi
  N=$(cat "$FILELIST" | wc -l | tr -d ' ')
  echo "Submitting $ds: $N volumes"
  sbatch --job-name="preproc-$ds" slurm_scripts/preprocess.sh "$FILELIST" "data/preprocessed/$ds/"
done
