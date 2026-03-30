#!/bin/bash
# Submit all preprocessing jobs
# Run from repo root: bash slurm_scripts/run_all.sh

cd <sanitized>/projects/brain-mri-preprocessing
source .venv/bin/activate

for ds in ixi oasis1 oasis2 adni stanford tcga schizo upenn; do
  FILELIST="data/staging/$ds/files.txt"
  OUTDIR="data/preprocessed/$ds"

  if [ ! -f "$FILELIST" ]; then
    echo "SKIP $ds: no files.txt"
    continue
  fi

  mkdir -p "$OUTDIR"

  echo "=== $ds ==="
  while read -r INPUT_FILE; do
    [ -z "$INPUT_FILE" ] && continue
    BASENAME=$(basename "$INPUT_FILE" .nii.gz)
    BASENAME=$(basename "$BASENAME" .nii)
    OUTPUT_FILE="${OUTDIR}/${BASENAME}_preprocessed.nii.gz"

    if [ -f "$OUTPUT_FILE" ]; then
      continue
    fi

    sbatch <sanitized> --gres=gpu:1 --cpus-per-task=4 --mem=16G --time=00:30:00 \
      --output=logs/${ds}_%j.out --error=logs/${ds}_%j.err \
      --wrap="cd <sanitized>/projects/brain-mri-preprocessing && source .venv/bin/activate && python preprocess.py '$INPUT_FILE' '$OUTPUT_FILE'"
  done < "$FILELIST"

  N=$(wc -l < "$FILELIST")
  echo "  Submitted $N jobs -> $OUTDIR"
done

mkdir -p logs
