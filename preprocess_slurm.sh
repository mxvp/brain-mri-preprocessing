#!/bin/bash
#SBATCH --job-name=mri-preprocess
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --array=0-99%20
#SBATCH --output=logs/preprocess_%A_%a.out
#SBATCH --error=logs/preprocess_%A_%a.err

# Usage:
#   Create a file list first:
#     find /path/to/raw/ -name "*.nii.gz" | sort > input_files.txt
#   Then submit:
#     sbatch preprocess_slurm.sh input_files.txt /path/to/output/
#   Adjust --array upper bound to match number of files (0-N where N = nfiles-1)

INPUT_LIST=$1
OUTPUT_DIR=$2

mkdir -p "$OUTPUT_DIR" logs

# Get the file for this array task
INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$INPUT_LIST")

if [ -z "$INPUT_FILE" ]; then
    echo "No file at line ${SLURM_ARRAY_TASK_ID}, exiting"
    exit 0
fi

BASENAME=$(basename "$INPUT_FILE" .nii.gz)
OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_preprocessed.nii.gz"

# Skip if already processed
if [ -f "$OUTPUT_FILE" ]; then
    echo "Already exists: $OUTPUT_FILE"
    exit 0
fi

echo "Processing: $INPUT_FILE -> $OUTPUT_FILE"

python preprocess.py "$INPUT_FILE" "$OUTPUT_FILE"
