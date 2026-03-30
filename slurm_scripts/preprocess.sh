#!/bin/bash
#SBATCH --job-name=mri-preprocess
#SBATCH <sanitized>
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# Usage:
#   sbatch --job-name=preproc-ixi slurm_scripts/preprocess.sh data/staging/ixi/manifest.json data/preprocessed/ixi/

MANIFEST=$1
OUTPUT_DIR=$2

cd <sanitized>/projects/brain-mri-preprocessing
source .venv/bin/activate
mkdir -p "$OUTPUT_DIR" logs

python preprocess.py --manifest "$MANIFEST" "$OUTPUT_DIR"
