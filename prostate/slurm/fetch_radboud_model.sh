#!/bin/bash
# Get the Radboud prostate-gland nnU-Net weights into a usable on-disk format.
#
# The official release is a Docker container (joeranbosma/picai_unet) that
# bundles the model. We extract the nnU-Net model directory from the
# container so nnUNetPredictor can use it directly — no Docker at inference.
#
# Output: $GROUP_SCRATCH/$USER/models/radboud_prostate_gland/  (model dir)
#
# Run interactively (needs sufficient disk in $GROUP_SCRATCH + apptainer or
# docker on the host):
#     bash prostate/slurm/fetch_radboud_model.sh
#
# If you don't have container tools available, see the alternative path at
# the bottom of this file.
set -euo pipefail

DEST=${RADBOUD_MODEL_DIR:-$GROUP_SCRATCH/$USER/models/radboud_prostate_gland}
mkdir -p "$DEST"

# --- Option A: apptainer (singularity) -- preferred on Sherlock ----------
if command -v apptainer >/dev/null; then
    SIF=$DEST/_radboud.sif
    if [ ! -f "$SIF" ]; then
        echo "Pulling container..."
        apptainer pull "$SIF" docker://joeranbosma/picai_unet:latest
    fi
    echo "Extracting model files from container into $DEST"
    apptainer exec "$SIF" bash -c \
        "cp -r /opt/algorithm/results/nnUNet/3d_fullres/Task2203_picai_baseline/* $DEST/"
    rm -f "$SIF"
    echo "Done. nnunet_checkpoint = $DEST"
    exit 0
fi

# --- Option B: docker ---------------------------------------------------
if command -v docker >/dev/null; then
    IMG=joeranbosma/picai_unet:latest
    docker pull "$IMG"
    CID=$(docker create "$IMG")
    docker cp "$CID:/opt/algorithm/results/nnUNet/3d_fullres/Task2203_picai_baseline/." "$DEST/"
    docker rm "$CID"
    echo "Done. nnunet_checkpoint = $DEST"
    exit 0
fi

cat <<MSG
Neither apptainer nor docker is available. Two manual alternatives:

  1. Run the official Docker container for inference (loses the local-CLI
     simplicity but works without manual model extraction). See
     https://github.com/DIAGNijmegen/AbdomenMRUS-prostate-segmentation#inference

  2. Train or borrow a separate prostate gland nnU-Net model and point
     segmentation.nnunet_checkpoint at that model folder. The Python API
     in segment.py is generic — any whole-gland binary nnU-Net works.

Exiting without weights.
MSG
exit 1
