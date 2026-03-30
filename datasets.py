"""Dataset registry for brain MRI file discovery and format conversion.

Each dataset subclass knows how to find brain MRI files in its specific directory
structure and convert them to NIfTI if needed. To add a new dataset, subclass
Dataset and register it in REGISTRY.

Usage:
    from datasets import REGISTRY
    dataset = REGISTRY["ixi"]()
    files = dataset.prepare(input_dir, output_dir)
"""

import logging
import subprocess
from abc import ABC, abstractmethod
from fnmatch import fnmatch
from pathlib import Path

import nibabel as nib
import numpy as np

log = logging.getLogger(__name__)


def _analyze_to_nifti(hdr_path: Path, output_path: Path):
    img = nib.load(hdr_path)
    data = img.get_fdata().squeeze().astype(np.float32)
    nib.save(nib.Nifti1Image(data, img.affine), output_path)


def _dcm2niix(dicom_dir: Path, output_dir: Path, filename: str) -> Path:
    result = subprocess.run(
        ["dcm2niix", "-z", "y", "-f", filename, "-o", str(output_dir), str(dicom_dir)],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"dcm2niix failed: {result.stderr}")
    candidates = sorted(output_dir.glob(f"{filename}*.nii.gz"))
    if not candidates:
        raise RuntimeError(f"dcm2niix produced no output for {dicom_dir}")
    return candidates[0]


class Dataset(ABC):
    """Base class. Subclass per dataset, implement prepare()."""

    name: str = ""

    @abstractmethod
    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        """Find T1w files, convert to NIfTI if needed.

        For datasets needing conversion, writes NIfTI to output_dir.
        For datasets already in NIfTI, returns paths to original files.

        Returns list of NIfTI paths ready for preprocessing.
        """
        ...


class IXI(Dataset):
    """IXI dataset — raw .nii.gz, structural modalities (T1, T2) in filename."""

    name = "ixi"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        files = sorted(
            f for f in input_dir.glob("*.nii.gz")
            if "-T1." in f.name or "-T2." in f.name
        )
        log.info(f"IXI: found {len(files)} structural volumes")
        return files


class OASIS1(Dataset):
    """OASIS-1 — Analyze format (.hdr/.img), multiple MPR runs per session."""

    name = "oasis1"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
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
                _analyze_to_nifti(hdr, out)
                log.info(f"Converted {subject_dir.name}")
            results.append(out)
        log.info(f"OASIS-1: {len(results)} T1w volumes")
        return results


class OASIS2(Dataset):
    """OASIS-2 — Analyze format, RAW/mpr-*.nifti.hdr."""

    name = "oasis2"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
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
                _analyze_to_nifti(hdrs[0], out)
                log.info(f"Converted {subject_dir.name}")
            results.append(out)
        log.info(f"OASIS-2: {len(results)} T1w volumes")
        return results


class OASIS3(Dataset):
    """OASIS-3 — BIDS NIfTI from oasis-scripts download (T1w)."""

    name = "oasis3"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        files = sorted(input_dir.rglob("*T1w.nii.gz"))
        log.info(f"OASIS-3: found {len(files)} T1w volumes")
        return files


class PPMI(Dataset):
    """PPMI — DICOM, multiple T1 naming conventions across sites.

    Handles two folder layouts:
      PPMI/{subject}/{modality}/...           (from zip 1)
      PPMI2/{subject}/{sub_id}/{modality}/... (from zip 2)
    """

    name = "ppmi"
    # Structural modality keywords (case-insensitive)
    STRUCTURAL_KEYWORDS = [
        "T1", "MPRAGE", "FSPGR", "SPGR", "BRAVO", "FFE", "TFE",  # T1 variants
        "T2", "FLAIR",                                              # T2/FLAIR
    ]
    # Keywords to exclude (functional, diffusion, fieldmaps, localizers, etc.)
    EXCLUDE_KEYWORDS = ["FMRI", "BOLD", "REST", "DTI", "DWI", "ASL", "LOC",
                        "FIELD", "B0", "B1", "GRE-MT", "GRE_MT", "GRE-NM",
                        "GRE_-_MT", "GRE-NMMT", "CALIBRATION", "SCOUT",
                        "SETTER", "SCREEN", "PLANE"]

    def _is_structural_dir(self, name: str) -> bool:
        upper = name.upper()
        if any(kw in upper for kw in self.EXCLUDE_KEYWORDS):
            return False
        return any(kw in upper for kw in self.STRUCTURAL_KEYWORDS)

    def _process_modality_dirs(self, subject_id: str, parent_dir: Path,
                                output_dir: Path, results: list[Path]):
        struct_dirs = [d for d in parent_dir.iterdir() if d.is_dir() and self._is_structural_dir(d.name)]
        if not struct_dirs:
            return

        for struct_dir in sorted(struct_dirs):
            # Sanitize modality name for filename
            modality = struct_dir.name.replace(" ", "_").replace("-", "_").lower()
            for date_dir in sorted(struct_dir.iterdir()):
                if not date_dir.is_dir():
                    continue
                dcm_dirs = [d for d in date_dir.iterdir() if d.is_dir()]
                if not dcm_dirs:
                    continue

                dcm_dir = dcm_dirs[0]  # contains .dcm files
                date_str = date_dir.name[:10].replace("-", "")
                filename = f"PPMI_{subject_id}_{modality}_{date_str}"
                out = output_dir / f"{filename}.nii.gz"

                if not out.exists():
                    try:
                        _dcm2niix(dcm_dir, output_dir, filename)
                    except Exception:
                        log.exception(f"PPMI {subject_id} {date_dir.name}")
                        continue
                if out.exists():
                    results.append(out)

    def _scan_root(self, root: Path, output_dir: Path, results: list[Path]):
        for subject_dir in sorted(root.iterdir()):
            if not subject_dir.is_dir():
                continue
            subject_id = subject_dir.name

            # Check if this has extra nesting (PPMI2 layout: subject/{sub_id}/{modality})
            subdirs = [d for d in subject_dir.iterdir() if d.is_dir()]
            if subdirs and all(d.name.isdigit() for d in subdirs[:5]):
                for sub_dir in sorted(subdirs):
                    self._process_modality_dirs(subject_id, sub_dir, output_dir, results)
            else:
                self._process_modality_dirs(subject_id, subject_dir, output_dir, results)

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        results = []

        # Handle multiple subdirs (PPMI/, PPMI2/) or direct subject dirs
        subdirs = sorted(d for d in input_dir.iterdir() if d.is_dir())
        if subdirs and subdirs[0].name.startswith("PPMI"):
            for sub_root in subdirs:
                # PPMI2/PPMI/{subjects} — unwrap extra PPMI dir if present
                inner = sub_root / "PPMI"
                if inner.is_dir():
                    log.info(f"Scanning {sub_root.name}/PPMI/")
                    self._scan_root(inner, output_dir, results)
                else:
                    log.info(f"Scanning {sub_root.name}/")
                    self._scan_root(sub_root, output_dir, results)
        else:
            self._scan_root(input_dir, output_dir, results)

        log.info(f"PPMI: {len(results)} T1w volumes")
        return results


class ADNI(Dataset):
    """ADNI (data_fusion) — T1.nii per subject directory."""

    name = "adni"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            for name in ("T1.nii.gz", "T1.nii"):
                t1 = subject_dir / name
                if t1.exists():
                    results.append(t1)
                    break
        log.info(f"ADNI: found {len(results)} T1w volumes")
        return results


class Schizo(Dataset):
    """SCHIZO — COBRE dataset. T1/T2 per subject (already skull-stripped, raw unavailable)."""

    name = "schizo"
    MODALITIES = ["T1.nii.gz", "T2.nii.gz"]

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            for mod in self.MODALITIES:
                f = subject_dir / mod
                if f.exists():
                    results.append(f)
        log.info(f"SCHIZO: found {len(results)} volumes (pre-skull-stripped)")
        return results


class Stanford(Dataset):
    """Stanford — brain tumor dataset with T1Gd, FLAIR per subject."""

    name = "stanford"
    # Raw modality files (not masks/intermediates)
    MODALITIES = ["T1Gd.nii.gz", "FLAIR.nii.gz"]

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            for mod in self.MODALITIES:
                f = subject_dir / mod
                if f.exists():
                    results.append(f)
        log.info(f"Stanford: found {len(results)} volumes")
        return results


class TCGA(Dataset):
    """TCGA — brain tumor dataset with t1, t1Gd, t2, flair per subject."""

    name = "tcga"
    # Pattern: TCGA-XX-XXXX_date_modality.nii.gz
    MODALITY_SUFFIXES = ["_t1.nii.gz", "_t1Gd.nii.gz", "_t2.nii.gz", "_flair.nii.gz"]

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            for f in sorted(subject_dir.glob("*.nii.gz")):
                if any(f.name.endswith(s) for s in self.MODALITY_SUFFIXES):
                    results.append(f)
        log.info(f"TCGA: found {len(results)} volumes")
        return results


class UCSF(Dataset):
    """UCSF-PDGM — already registered to SRI24. Structural modalities with bias correction."""

    name = "ucsf"
    # Prefer bias-corrected versions, all structural modalities
    MODALITIES = ["_T1_bias.nii.gz", "_T1c_bias.nii.gz", "_T2_bias.nii.gz", "_FLAIR_bias.nii.gz"]

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            for mod in self.MODALITIES:
                matches = sorted(subject_dir.glob(f"*{mod}"))
                results.extend(matches)
        log.info(f"UCSF: found {len(results)} volumes (already in SRI24)")
        return results


class UPenn(Dataset):
    """UPenn-GBM — raw structural NIfTI: T1, T1GD, T2, FLAIR per subject."""

    name = "upenn"
    MODALITY_SUFFIXES = ["_T1.nii.gz", "_T1GD.nii.gz", "_T2.nii.gz", "_FLAIR.nii.gz"]

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            for f in sorted(subject_dir.glob("*.nii.gz")):
                if any(f.name.endswith(s) for s in self.MODALITY_SUFFIXES):
                    results.append(f)
        log.info(f"UPenn: found {len(results)} volumes")
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
