#!/bin/bash
# Submit one preprocessing job per dataset
# Run from repo root: bash slurm_scripts/run_all.sh

mkdir -p logs

for ds in ixi oasis1 oasis2 adni stanford tcga schizo upenn ppmi; do
  MANIFEST="data/staging/$ds/manifest.json"
  if [ ! -f "$MANIFEST" ]; then
    echo "SKIP $ds: no manifest.json"
    continue
  fi
  N=$(python -c "import json; print(len(json.load(open('$MANIFEST'))))")
  echo "Submitting $ds: $N subjects"
  sbatch --job-name="preproc-$ds" slurm_scripts/preprocess.sh "$MANIFEST" "data/preprocessed/$ds/"
done
