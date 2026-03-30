#!/bin/bash
# Submit one preprocessing job per dataset
# Run from repo root: bash slurm_scripts/run_all.sh

for ds in ixi oasis1 oasis2 adni stanford tcga schizo upenn; do
  FILELIST="data/staging/$ds/files.txt"

  if [ ! -f "$FILELIST" ]; then
    echo "SKIP $ds: no files.txt"
    continue
  fi

  N=$(cat "$FILELIST" | wc -l | tr -d ' ')
  echo "Submitting $ds: $N volumes"

  sbatch --job-name="preproc-$ds" \
    <sanitized> --gres=gpu:1 --cpus-per-task=4 --mem=16G \
    --time=24:00:00 \
    --output=logs/${ds}_%j.out --error=logs/${ds}_%j.err \
    --wrap="cd <sanitized>/projects/brain-mri-preprocessing && source .venv/bin/activate && python preprocess.py --filelist $FILELIST --output data/preprocessed/$ds/"
done

mkdir -p logs
