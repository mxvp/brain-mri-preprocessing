#!/bin/bash
#SBATCH --job-name=mri-preprocess-cpu
#SBATCH <sanitized>
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

# CPU-only preprocessing (no GPU). For pre-skull-stripped datasets that only need registration.

MANIFEST=$1
OUTPUT_DIR=$2

cd <sanitized>/projects/brain-mri-preprocessing
source .venv/bin/activate
mkdir -p "$OUTPUT_DIR" logs

python preprocess.py --manifest "$MANIFEST" --output "$OUTPUT_DIR" --device cpu
