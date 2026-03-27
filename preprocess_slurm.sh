#!/bin/bash
#SBATCH --job-name=mri-preprocess
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/preprocess_%A_%a.out
#SBATCH --error=logs/preprocess_%A_%a.err

# Usage:
#   python prepare.py ixi data/IXI staging/IXI        # generates staging/IXI/files.txt
#   N=$(wc -l < staging/IXI/files.txt)
#   sbatch --array=1-${N}%20 preprocess_slurm.sh staging/IXI/files.txt /path/to/output/

INPUT_LIST=$1
OUTPUT_DIR=$2

mkdir -p "$OUTPUT_DIR" logs

# Get the file for this array task (files.txt is 1-indexed)
INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")

if [ -z "$INPUT_FILE" ]; then
    echo "No file at line ${SLURM_ARRAY_TASK_ID}, exiting"
    exit 0
fi

BASENAME=$(basename "$INPUT_FILE" .nii.gz)
BASENAME=$(basename "$BASENAME" .nii)
OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_preprocessed.nii.gz"

if [ -f "$OUTPUT_FILE" ]; then
    echo "Already exists: $OUTPUT_FILE"
    exit 0
fi

echo "Processing: $INPUT_FILE -> $OUTPUT_FILE"
python preprocess.py "$INPUT_FILE" "$OUTPUT_FILE"
