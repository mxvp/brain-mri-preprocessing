"""Render warped tumor masks over their preprocessed volumes for visual QC.

Companion to warp_masks.py. Samples subjects from each cohort, finds the
slice carrying the most tumor, and overlays the warped mask on the
contrast-enhanced volume — where enhancing tumor is actually visible, so a
mask that landed in the wrong place is obvious rather than merely plausible.

Sampling is seeded, and specific subjects can be forced in with --include,
so the ones flagged by warp_report.csv can be inspected alongside a random
baseline in the same figure.

Usage:
    python qc_mask_overlay.py <masks_root> <preprocessed_root> out.png
    python qc_mask_overlay.py <masks_root> <preprocessed_root> out.png \
        --n-per-cohort 6 --include UPENN-GBM-00510_11 UCSF-PDGM-0094
"""
from __future__ import annotations

import argparse
import csv
import logging
import random
from pathlib import Path

log = logging.getLogger(__name__)

# Which preprocessed modality to draw underneath. Contrast-enhanced, so the
# enhancing rim the mask claims to outline is visible in the same image.
COHORT_UNDERLAY = {
    "upenn": ("UPenn_{sid}_t1gd_preprocessed.nii.gz", "UPenn_{sid}_segm_preprocessed.nii.gz"),
    "ucsf":  ("{sid}_t1c_preprocessed.nii.gz",        "{sid}_segm_preprocessed.nii.gz"),
}

# BraTS labels: 1 necrotic core, 2 edema, 4 enhancing. Index 3 is unused.
LABEL_COLORS = [(0, 0, 0, 0), (0, 1, 1, 0.55), (1, 1, 0, 0.40),
                (0, 0, 0, 0), (1, 0, 0, 0.65)]


def _read_report(masks_root: Path, cohort: str) -> dict[str, dict]:
    path = masks_root / cohort / "warp_report.csv"
    if not path.exists():
        return {}
    with open(path) as f:
        return {r["subject"]: r for r in csv.DictReader(f)}


def _pick(masks_root: Path, cohort: str, n: int, include: list[str], seed: int) -> list[str]:
    """Forced subjects first, then a seeded random sample of the rest."""
    # Recover subject ids by stripping the fixed prefix/suffix of the template.
    prefix, suffix = COHORT_UNDERLAY[cohort][1].split("{sid}")
    available = sorted(
        p.name[len(prefix):-len(suffix)]
        for p in (masks_root / cohort).glob("*_segm_preprocessed.nii.gz")
        if p.name.startswith(prefix) and p.name.endswith(suffix)
    )
    forced = [s for s in include if s in available]
    rest = [s for s in available if s not in forced]
    rng = random.Random(seed)
    sampled = rng.sample(rest, min(n, len(rest)))
    return forced + sampled


def render(masks_root: Path, prep_root: Path, out_png: Path,
           n_per_cohort: int = 6, include: list[str] | None = None, seed: int = 0):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import nibabel as nib
    import numpy as np
    from matplotlib.colors import ListedColormap

    include = include or []
    cmap = ListedColormap(LABEL_COLORS)

    panels = []
    for cohort in ("upenn", "ucsf"):
        if not (masks_root / cohort).is_dir():
            log.warning(f"no {cohort} dir under {masks_root}, skipping")
            continue
        report = _read_report(masks_root, cohort)
        img_tmpl, mask_tmpl = COHORT_UNDERLAY[cohort]
        for sid in _pick(masks_root, cohort, n_per_cohort, include, seed):
            img_p = prep_root / cohort / img_tmpl.format(sid=sid)
            msk_p = masks_root / cohort / mask_tmpl.format(sid=sid)
            if not img_p.exists():
                log.warning(f"{sid}: no underlay at {img_p.name}")
                continue
            panels.append((cohort, sid, img_p, msk_p, report.get(sid, {})))

    if not panels:
        raise SystemExit("nothing to render — check the masks/preprocessed paths")

    ncol = 4
    nrow = (len(panels) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.5 * nrow),
                             constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, (cohort, sid, img_p, msk_p, rep) in zip(axes, panels):
        img = np.asarray(nib.load(img_p).dataobj).astype(np.float32)
        msk = np.asarray(nib.load(msk_p).dataobj).astype(np.float32)
        # Slice with the most tumor — the most informative single view.
        per_slice = (msk > 0).sum(axis=(0, 1))
        z = int(per_slice.argmax()) if per_slice.max() else img.shape[2] // 2

        sl_i, sl_m = img[:, :, z], msk[:, :, z]
        nz = sl_i[sl_i != 0]
        vmin, vmax = (np.percentile(nz, [2, 98]) if nz.size else (0, 1))
        ax.imshow(np.rot90(sl_i), cmap="gray", vmin=vmin, vmax=vmax)
        ax.imshow(np.rot90(sl_m), cmap=cmap, vmin=0, vmax=4, interpolation="nearest")

        pct = rep.get("pct_inside_brain", "?")
        labs = sorted(int(v) for v in np.unique(msk) if v)
        flag = "" if rep.get("labels_match", "True") == "True" else "  [LABELS CHANGED]"
        ax.set_title(f"{sid}\nz={z}  inside={pct}%  labels={labs}{flag}", fontsize=6.5)
        ax.axis("off")

    for ax in axes[len(panels):]:
        ax.axis("off")

    fig.suptitle("Warped tumor masks over preprocessed volumes  "
                 "(cyan=necrotic  yellow=edema  red=enhancing)", fontsize=10)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=115)
    log.info(f"wrote {out_png}  ({len(panels)} panels)")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Visual QC for warped tumor masks")
    p.add_argument("masks_root", type=Path, help="dir containing upenn/ and ucsf/ mask subdirs")
    p.add_argument("preprocessed_root", type=Path, help="dir containing upenn/ and ucsf/ preprocessed subdirs")
    p.add_argument("out_png", type=Path)
    p.add_argument("--n-per-cohort", type=int, default=6, help="random subjects per cohort (on top of --include)")
    p.add_argument("--include", nargs="*", default=[], help="subject ids to force into the figure")
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    render(args.masks_root, args.preprocessed_root, args.out_png,
           args.n_per_cohort, args.include, args.seed)
