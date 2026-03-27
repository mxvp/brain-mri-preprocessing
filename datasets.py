"""Dataset registry for T1w file discovery and format conversion.

Each dataset subclass knows how to find T1w files in its specific directory
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
    """IXI dataset — raw T1w .nii.gz, filename contains '-T1'."""

    name = "ixi"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        files = sorted(input_dir.glob("*-T1.nii.gz"))
        log.info(f"IXI: found {len(files)} T1w volumes")
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
    """PPMI — DICOM, T1w dirs match *T1*weighted* or *MPRAGE*."""

    name = "ppmi"
    T1_PATTERNS = ["*T1*weighted*", "*T1-weighted*", "*MPRAGE*"]

    def _is_t1_dir(self, name: str) -> bool:
        name_upper = name.upper()
        return any(fnmatch(name_upper, p.upper()) for p in self.T1_PATTERNS)

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        output_dir.mkdir(parents=True, exist_ok=True)
        results = []

        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            subject_id = subject_dir.name

            t1_dirs = [d for d in subject_dir.iterdir() if d.is_dir() and self._is_t1_dir(d.name)]
            if not t1_dirs:
                continue

            for t1_dir in sorted(t1_dirs):
                for date_dir in sorted(t1_dir.iterdir()):
                    if not date_dir.is_dir():
                        continue
                    dcm_dirs = [d for d in date_dir.iterdir() if d.is_dir() and list(d.glob("*.dcm"))[:1]]
                    if not dcm_dirs:
                        continue

                    dcm_dir = dcm_dirs[0]
                    date_str = date_dir.name[:10].replace("-", "")
                    filename = f"PPMI_{subject_id}_{date_str}"
                    out = output_dir / f"{filename}.nii.gz"

                    if not out.exists():
                        try:
                            _dcm2niix(dcm_dir, output_dir, filename)
                        except Exception:
                            log.exception(f"PPMI {subject_id} {date_dir.name}")
                            continue
                    if out.exists():
                        results.append(out)

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
    """SCHIZO (data_fusion) — T1.nii.gz per subject, already skull-stripped."""

    name = "schizo"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            t1 = subject_dir / "T1.nii.gz"
            if t1.exists():
                results.append(t1)
        log.info(f"SCHIZO: found {len(results)} T1w volumes")
        return results


class Stanford(Dataset):
    """Stanford (data_fusion) — T1.nii.gz per subject, already skull-stripped."""

    name = "stanford"

    def prepare(self, input_dir: Path, output_dir: Path) -> list[Path]:
        results = []
        for subject_dir in sorted(input_dir.iterdir()):
            if not subject_dir.is_dir():
                continue
            t1 = subject_dir / "T1.nii.gz"
            if t1.exists():
                results.append(t1)
        log.info(f"Stanford: found {len(results)} T1w volumes")
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
}
