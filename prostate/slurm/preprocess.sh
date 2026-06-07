#!/bin/bash
# Preprocess all complete prostate patients on a single GPU node.
#
# Submit:  sbatch prostate/slurm/preprocess.sh
# Monitor: tail -f prostate/slurm/logs/preprocess_<jobid>.out
#
# Layout assumptions:
#   - dcm2niix output:  $GROUP_OAK/maxvpuyv/projects/prostate_dt/data/.../nifti
#   - nnU-Net model:    $GROUP_SCRATCH/$USER/models/radboud_prostate_gland
#   - Outputs:          $GROUP_SCRATCH/$USER/prostate/preprocessed
#   - Conda env:        brain  (has SimpleITK, nibabel, nnunetv2, torch)

#SBATCH --job-name=prostate_pp
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=prostate/slurm/logs/preprocess_%j.out
#SBATCH --error=prostate/slurm/logs/preprocess_%j.err

set -euo pipefail

REPO=<sanitized>/projects/brain-mri-preprocessing
DATA=<sanitized>/projects/prostate_dt/data/Prostate-MRI-US-Biopsy-NBIA-manifest_v2_20231020/nifti
OUT=$GROUP_SCRATCH/$USER/prostate/preprocessed
MODEL=$GROUP_SCRATCH/$USER/models/radboud_prostate_gland
CFG=$SLURM_SUBMIT_DIR/_prostate_runtime_config.yaml

mkdir -p "$OUT" "$(dirname "$CFG")"

# Materialize the config with absolute paths.
cat > "$CFG" <<EOF
paths:
  nifti_root: $DATA
  output_root: $OUT
series:
  t2:    ['t2.*prostate', 't2_spc', 'ax_t2', 't2_tra']
  adc:   ['adc']
  dwi:   ['calc_bval', 'high_b', 'tracew']
n4: {shrink_factor: 4, n_iterations: [50,50,50,50], apply_to: [t2, adc]}
registration: {type: rigid, reference: t2, interpolator: linear}
segmentation:
  backend: nnunet
  nnunet_checkpoint: $MODEL
  nnunet_folds: all
  nnunet_checkpoint_name: checkpoint_final.pth
crop: {margin_mm: [20, 20, 6]}
resample: {voxel_spacing: [0.5, 0.5, 3.0]}
normalize:
  t2:  {method: zscore_nonzero}
  dwi: {method: zscore_nonzero}
  adc: {method: clip_then_zscore, clip_low: 0.5, clip_high: 99.5}
output: {channels: [t2, dwi, adc, mask], dtype: float32, format: nifti}
EOF

source <sanitized>/anaconda3/etc/profile.d/conda.sh
conda activate brain

cd "$REPO"
echo "[$(date)] starting preprocess"
echo "device: $(python -c 'import torch; print(torch.cuda.get_device_name() if torch.cuda.is_available() else \"cpu\")')"
python -m prostate inventory --config "$CFG"
python -m prostate preprocess --config "$CFG"
echo "[$(date)] done"
ls -lh "$OUT" | head
echo "total outputs: $(ls "$OUT" | wc -l)"
