"""Dataset registry for brain MRI file discovery and format conversion.

Each dataset subclass knows how to find brain MRI files in its specific directory
structure, convert them to NIfTI if needed, and group them by subject with a
designated center modality (T1/T1Gd) for atlas registration.

Usage:
    from datasets import REGISTRY
    dataset = REGISTRY["ixi"]()
    subjects = dataset.prepare(input_dir, output_dir)
"""

import json
import logging
import subprocess
from abc import ABC, abstractmethod
from fnmatch import fnmatch
from pathlib import Path

import ants
import nibabel as nib
import numpy as np

log = logging.getLogger(__name__)


# --- Conversion helpers ---

def _oasis_analyze_to_nifti(hdr_path: Path, output_path: Path):
    """Convert OASIS Analyze (.hdr/.img) to NIfTI with correct orientation.

    OASIS-1/2 sagittal MPRAGE data is stored as (SI, AP, LR) but the Analyze
    header misreports the orientation. We permute axes, set correct metadata,
    and pre-register to SRI24 using direct ANTs (brainles registration doesn't
    converge well from the Analyze starting position).
    """
    # Load SRI24 atlas from brainles installation
    import brainles_preprocessing.registration as _reg
    sri24_paths = list(Path(_reg.__file__).parent.rglob("sri24.nii"))
    if not sri24_paths:
        raise RuntimeError("SRI24 atlas not found. Run brainles once to download it.")
    atlas = ants.image_read(str(sri24_paths[0]))

    img = ants.image_read(str(hdr_path))
    arr = img.numpy()  # (256, 256, 128) = (SI, AP, LR)

    # Permute to (LR, AP, SI)
    arr_ras = np.transpose(arr, (2, 1, 0)).astype(np.float32).copy()
    spacing = (img.spacing[2], img.spacing[1], img.spacing[0])

    new_img = ants.from_numpy(
        arr_ras,
        origin=(arr_ras.shape[0] * spacing[0], 0, 0),
        spacing=spacing,
        direction=atlas.direction,
    )

    # Pre-register to SRI24 — direct ANTs handles this reliably
    result = ants.registration(fixed=atlas, moving=new_img, type_of_transform="Affine")
    ants.image_write(result["warpedmovout"], str(output_path))


def _dcm2niix(dicom_dir: Path, output_dir: Path, filename: str) -> Path:
    result = subprocess.run(
        ["dcm2niix", "-z", "y", "-f", filename, "-o", str(output_dir), str(dicom_dir)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"dcm2niix failed: {result.stderr}")
    # Pick the first output, skip _Eq_* duplicates
    candidates = sorted(
        f for f in output_dir.glob(f"{filename}*.nii.gz")
        if "_Eq_" not in f.name
    )
    if not candidates:
        raise RuntimeError(f"dcm2niix produced no output for {dicom_dir}")
    return candidates[0]


def _is_3d_volume(path: Path, min_slices: int = 50, max_spacing: float = 3.0) -> bool:
    """Check if a NIfTI volume is a usable 3D acquisition (not 2D thick-slice)."""
    img = nib.load(path)
    shape = img.header.get_data_shape()[:3]
    zooms = img.header.get_zooms()[:3]
    return all(s >= min_slices for s in shape) and all(z <= max_spacing for z in zooms)


# --- Subject manifest type ---

# Each subject entry:
# {
#   "subject_id": str,
#   "center": {"modality": str, "path": str},
#   "moving": [{"modality": str, "path": str}, ...]
# }


class Dataset(ABC):
    """Base class. Subclass per dataset, implement prepare()."""

    name: str = ""

    @abstractmethod
    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        """Find brain MRI files, convert to NIfTI if needed, group by subject.

        Returns list of subject dicts with center + moving modalities.
        """
        ...


# --- Helper to identify center modality ---

T1_NAMES = {"t1", "t1gd", "t1_gd", "t1cd", "t1c", "t1w", "mprage", "fspgr",
            "spgr", "bravo", "ffe", "tfe"}
T2_NAMES = {"t2", "t2w", "flair"}


def _pick_center(modalities: dict[str, Path]) -> tuple[str, dict, list[dict]]:
    """Given {modality_name: path}, pick the best center modality.
    Prefers T1 > T1Gd > T1c. Returns (subject_id_unused, center_dict, moving_list)."""
    # Priority order for center
    for pref in ["t1", "t1w", "t1gd", "t1_gd", "t1cd", "t1c"]:
        if pref in modalities:
            center_mod = pref
            center_path = modalities[center_mod]
            moving = [{"modality": m, "path": str(p)} for m, p in modalities.items() if m != center_mod]
            return {"modality": center_mod, "path": str(center_path)}, moving

    # Fallback: first T1-like modality
    for mod, path in modalities.items():
        if any(t in mod.lower() for t in T1_NAMES):
            moving = [{"modality": m, "path": str(p)} for m, p in modalities.items() if m != mod]
            return {"modality": mod, "path": str(path)}, moving

    # Last resort: first modality
    first_mod = next(iter(modalities))
    moving = [{"modality": m, "path": str(p)} for m, p in modalities.items() if m != first_mod]
    return {"modality": first_mod, "path": str(modalities[first_mod])}, moving


# --- Dataset implementations ---

class IXI(Dataset):
    """IXI — raw .nii.gz, structural modalities (T1, T2) in filename."""

    name = "ixi"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        files = sorted(
            f for f in input_dir.glob("*.nii.gz")
            if "-T1." in f.name or "-T2." in f.name
        )
        # Group by subject: IXI002-Guys-0828-T1.nii.gz -> subject=IXI002-Guys-0828
        subjects = {}
        for f in files:
            parts = f.stem.replace(".nii", "").rsplit("-", 1)
            subject_id = parts[0]
            modality = parts[1].lower()
            subjects.setdefault(subject_id, {})[modality] = f

        results = []
        for subject_id, mods in sorted(subjects.items()):
            center, moving = _pick_center(mods)
            results.append({"subject_id": subject_id, "center": center, "moving": moving})

        log.info(f"IXI: {len(results)} subjects, {sum(1 + len(s['moving']) for s in results)} volumes")
        return results


class OASIS1(Dataset):
    """OASIS-1 — Analyze format (.hdr/.img), sagittal MPRAGE, axis permutation required."""

    name = "oasis1"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        output_dir.mkdir(parents=True, exist_ok=True)
        results = []
        for subject_dir in sorted(input_dir.rglob("OAS1_*_MR*")):
            raw_dir = subject_dir / "RAW"
            if not raw_dir.exists():
                continue
            hdr = raw_dir / f"{subject_dir.name}_mpr-1_anon.hdr"
            if not hdr.exists():
                hdrs = sorted(raw_dir.glob("*mpr-1*.hdr"))
                if not hdrs:
                    log.warning(f"No mpr-1 for {subject_dir.name}")
                    continue
                hdr = hdrs[0]
            out = output_dir / f"{subject_dir.name}.nii.gz"
            if not out.exists():
                _oasis_analyze_to_nifti(hdr, out)
                log.info(f"Converted {subject_dir.name}")
            results.append({
                "subject_id": subject_dir.name,
                "center": {"modality": "t1", "path": str(out)},
                "moving": [],
            })
        log.info(f"OASIS-1: {len(results)} subjects")
        return results


class OASIS2(Dataset):
    """OASIS-2 — Analyze format, same sagittal MPRAGE fix as OASIS-1."""

    name = "oasis2"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        output_dir.mkdir(parents=True, exist_ok=True)
        results = []
        for subject_dir in sorted(input_dir.rglob("OAS2_*_MR*")):
            raw_dir = subject_dir / "RAW"
            if not raw_dir.exists():
                continue
            hdrs = sorted(raw_dir.glob("mpr-1*.*hdr"))
            if not hdrs:
                log.warning(f"No mpr-1 for {subject_dir.name}")
                continue
            out = output_dir / f"{subject_dir.name}.nii.gz"
            if not out.exists():
                _oasis_analyze_to_nifti(hdrs[0], out)
                log.info(f"Converted {subject_dir.name}")
            results.append({
                "subject_id": subject_dir.name,
                "center": {"modality": "t1", "path": str(out)},
                "moving": [],
            })
        log.info(f"OASIS-2: {len(results)} subjects")
        return results


class OASIS3(Dataset):
    """OASIS-3 — BIDS NIfTI from oasis-scripts download (T1w)."""

    name = "oasis3"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        files = sorted(input_dir.rglob("*T1w.nii.gz"))
        results = [
            {"subject_id": f.stem.replace(".nii", ""), "center": {"modality": "t1", "path": str(f)}, "moving": []}
            for f in files
        ]
        log.info(f"OASIS-3: {len(results)} subjects")
        return results


class PPMI(Dataset):
    """PPMI — DICOM, multiple naming conventions. Only 3D structural sequences."""

    name = "ppmi"
    STRUCTURAL_KEYWORDS = [
        "T1", "MPRAGE", "FSPGR", "SPGR", "BRAVO", "FFE", "TFE",
        "T2", "FLAIR",
    ]
    EXCLUDE_KEYWORDS = ["FMRI", "BOLD", "REST", "DTI", "DWI", "ASL", "LOC",
                        "FIELD", "B0", "B1", "GRE-MT", "GRE_MT", "GRE-NM",
                        "GRE_-_MT", "GRE-NMMT", "CALIBRATION", "SCOUT",
                        "SETTER", "SCREEN", "PLANE"]

    def _is_structural_dir(self, name: str) -> bool:
        upper = name.upper()
        if any(kw in upper for kw in self.EXCLUDE_KEYWORDS):
            return False
        return any(kw in upper for kw in self.STRUCTURAL_KEYWORDS)

    def _infer_modality(self, dirname: str) -> str:
        upper = dirname.upper()
        if "FLAIR" in upper:
            return "flair"
        if any(kw in upper for kw in ["T2", "TSE", "FSE"]):
            return "t2"
        return "t1"

    def _process_modality_dirs(self, subject_id: str, parent_dir: Path,
                                output_dir: Path, subject_mods: dict):
        struct_dirs = [d for d in parent_dir.iterdir() if d.is_dir() and self._is_structural_dir(d.name)]
        for struct_dir in sorted(struct_dirs):
            modality = self._infer_modality(struct_dir.name)
            mod_label = struct_dir.name.replace(" ", "_").replace("-", "_").lower()
            for date_dir in sorted(struct_dir.iterdir()):
                if not date_dir.is_dir():
                    continue
                dcm_dirs = [d for d in date_dir.iterdir() if d.is_dir()]
                if not dcm_dirs:
                    continue
                dcm_dir = dcm_dirs[0]
                date_str = date_dir.name[:10].replace("-", "")
                filename = f"PPMI_{subject_id}_{mod_label}_{date_str}"
                out = output_dir / f"{filename}.nii.gz"
                if not out.exists():
                    try:
                        _dcm2niix(dcm_dir, output_dir, filename)
                    except Exception:
                        log.exception(f"PPMI {subject_id} {date_dir.name}")
                        continue
                if out.exists() and _is_3d_volume(out):
                    key = f"{modality}_{date_str}"
                    subject_mods[key] = (modality, out)
                elif out.exists():
                    log.info(f"Skipping 2D: {out.name}")

    def _scan_root(self, root: Path, output_dir: Path, results: list[dict]):
        for subject_dir in sorted(root.iterdir()):
            if not subject_dir.is_dir():
                continue
            subject_id = subject_dir.name
            subject_mods = {}
            subdirs = [d for d in subject_dir.iterdir() if d.is_dir()]
            if subdirs and all(d.name.isdigit() for d in subdirs[:5]):
                for sub_dir in sorted(subdirs):
                    self._process_modality_dirs(subject_id, sub_dir, output_dir, subject_mods)
            else:
                self._process_modality_dirs(subject_id, subject_dir, output_dir, subject_mods)

            if not subject_mods:
                continue

            # Group by date, pick center per session
            mods_by_name = {k: path for k, (mod, path) in subject_mods.items()}
            mods_typed = {k: (mod, path) for k, (mod, path) in subject_mods.items()}

            # Find T1 center
            t1_entries = {k: p for k, (m, p) in mods_typed.items() if m == "t1"}
            other_entries = {k: (m, p) for k, (m, p) in mods_typed.items() if m != "t1"}

            for t1_key, t1_path in t1_entries.items():
                date = t1_key.split("_")[-1]
                moving = [{"modality": m, "path": str(p)}
                          for k, (m, p) in other_entries.items() if date in k]
                results.append({
                    "subject_id": f"PPMI_{subject_id}_{date}",
                    "center": {"modality": "t1", "path": str(t1_path)},
                    "moving": moving,
                })

            # T1-less sessions (standalone T2/FLAIR) — register independently
            t1_dates = {k.split("_")[-1] for k in t1_entries}
            for k, (m, p) in other_entries.items():
                date = k.split("_")[-1]
                if date not in t1_dates:
                    results.append({
                        "subject_id": f"PPMI_{subject_id}_{date}",
                        "center": {"modality": m, "path": str(p)},
                        "moving": [],
                    })

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        output_dir.mkdir(parents=True, exist_ok=True)
        results = []
        subdirs = sorted(d for d in input_dir.iterdir() if d.is_dir())
        if subdirs and subdirs[0].name.startswith("PPMI"):
            for sub_root in subdirs:
                inner = sub_root / "PPMI"
                if inner.is_dir():
                    log.info(f"Scanning {sub_root.name}/PPMI/")
                    self._scan_root(inner, output_dir, results)
                else:
                    log.info(f"Scanning {sub_root.name}/")
                    self._scan_root(sub_root, output_dir, results)
        else:
            self._scan_root(input_dir, output_dir, results)
        log.info(f"PPMI: {len(results)} sessions")
        return results


class ADNI(Dataset):
    """ADNI — T1.nii per subject directory."""

    name = "adni"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            for name in ("T1.nii.gz", "T1.nii"):
                t1 = subject_dir / name
                if t1.exists():
                    results.append({
                        "subject_id": f"ADNI_{subject_dir.name}",
                        "center": {"modality": "t1", "path": str(t1)},
                        "moving": [],
                    })
                    break
        log.info(f"ADNI: {len(results)} subjects")
        return results


class Schizo(Dataset):
    """SCHIZO — COBRE dataset. T1/T2 per subject (already skull-stripped)."""

    name = "schizo"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            t1 = subject_dir / "T1.nii.gz"
            if not t1.exists():
                continue
            moving = []
            t2 = subject_dir / "T2.nii.gz"
            if t2.exists():
                moving.append({"modality": "t2", "path": str(t2)})
            results.append({
                "subject_id": f"SCHIZO_{subject_dir.name}",
                "center": {"modality": "t1", "path": str(t1)},
                "moving": moving,
            })
        log.info(f"SCHIZO: {len(results)} subjects")
        return results


class Stanford(Dataset):
    """Stanford — brain tumor dataset with T1Gd, FLAIR per subject."""

    name = "stanford"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            t1gd = subject_dir / "T1Gd.nii.gz"
            if not t1gd.exists():
                continue
            moving = []
            flair = subject_dir / "FLAIR.nii.gz"
            if flair.exists():
                moving.append({"modality": "flair", "path": str(flair)})
            results.append({
                "subject_id": f"Stanford_{subject_dir.name}",
                "center": {"modality": "t1gd", "path": str(t1gd)},
                "moving": moving,
            })
        log.info(f"Stanford: {len(results)} subjects")
        return results


class TCGA(Dataset):
    """TCGA — brain tumor dataset with t1, t1Gd, t2, flair per subject."""

    name = "tcga"
    MODALITY_SUFFIXES = {
        "_t1.nii.gz": "t1",
        "_t1Gd.nii.gz": "t1gd",
        "_t2.nii.gz": "t2",
        "_flair.nii.gz": "flair",
    }

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            mods = {}
            for f in sorted(subject_dir.glob("*.nii.gz")):
                for suffix, mod_name in self.MODALITY_SUFFIXES.items():
                    if f.name.endswith(suffix):
                        mods[mod_name] = f
            if not mods:
                continue
            center, moving = _pick_center(mods)
            results.append({
                "subject_id": f"TCGA_{subject_dir.name}",
                "center": center,
                "moving": moving,
            })
        log.info(f"TCGA: {len(results)} subjects")
        return results


class UCSF(Dataset):
    """UCSF-PDGM — already registered to SRI24. Structural modalities with bias correction."""

    name = "ucsf"
    MODALITY_MAP = {
        "_T1_bias.nii.gz": "t1",
        "_T1c_bias.nii.gz": "t1c",
        "_T2_bias.nii.gz": "t2",
        "_FLAIR_bias.nii.gz": "flair",
    }

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            mods = {}
            for suffix, mod_name in self.MODALITY_MAP.items():
                matches = sorted(subject_dir.glob(f"*{suffix}"))
                if matches:
                    mods[mod_name] = matches[0]
            if not mods:
                continue
            center, moving = _pick_center(mods)
            subject_id = subject_dir.name.replace("_nifti", "")
            results.append({
                "subject_id": subject_id,
                "center": center,
                "moving": moving,
            })
        log.info(f"UCSF: {len(results)} subjects (already in SRI24)")
        return results


class UPenn(Dataset):
    """UPenn-GBM — raw structural NIfTI: T1, T1GD, T2, FLAIR per subject."""

    name = "upenn"
    MODALITY_SUFFIXES = {
        "_T1.nii.gz": "t1",
        "_T1GD.nii.gz": "t1gd",
        "_T2.nii.gz": "t2",
        "_FLAIR.nii.gz": "flair",
    }

    def prepare(self, input_dir: Path, output_dir: Path) -> list[dict]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            mods = {}
            for f in sorted(subject_dir.glob("*.nii.gz")):
                for suffix, mod_name in self.MODALITY_SUFFIXES.items():
                    if f.name.endswith(suffix):
                        mods[mod_name] = f
            if not mods:
                continue
            center, moving = _pick_center(mods)
            results.append({
                "subject_id": f"UPenn_{subject_dir.name}",
                "center": center,
                "moving": moving,
            })
        log.info(f"UPenn: {len(results)} subjects")
        return results


REGISTRY: dict[str, type[Dataset]] = {
    "ixi": IXI,
    "oasis1": OASIS1,
    "oasis2": OASIS2,
    "oasis3": OASIS3,
    "ppmi": PPMI,
    "adni": ADNI,
    "schizo": Schizo,
    "stanford": Stanford,
    "tcga": TCGA,
    "ucsf": UCSF,
    "upenn": UPenn,
}
