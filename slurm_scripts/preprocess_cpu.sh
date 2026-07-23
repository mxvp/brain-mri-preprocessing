#!/bin/bash
#SBATCH --job-name=mri-preprocess-cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# Add a partition for your cluster, e.g.:
##SBATCH --partition=your_partition_name

# CPU-only preprocessing (no GPU). For pre-skull-stripped datasets that only need registration.
# Submit from the repo root — SLURM_SUBMIT_DIR is set automatically to that dir.

MANIFEST=$1
OUTPUT_DIR=$2

cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
source .venv/bin/activate
mkdir -p "$OUTPUT_DIR" logs

python preprocess.py --manifest "$MANIFEST" --output "$OUTPUT_DIR" --device cpu
