"""Parametrized tests for parse_latent_path.

Covers every cohort-ID variant we know about, plus negative cases. New
cohorts/encodings should add a parametrize row here so we catch silent
join failures before they leak into modeling.
"""
import pytest

from multimodal.load import parse_latent_path


# (path, expected_project_id, expected_submitter_id)
ID_CASES = [
    # CPTAC: BIDS sub- prefix + nnUNet _0000 channel suffix
    ("temp/cptac/sub-C3L00016_0000.nii.gz",                  "CPTAC-3",   "C3L-00016"),
    ("temp/cptac/sub-C3N02256_0000.nii.gz",                  "CPTAC-3",   "C3N-02256"),
    # CPTAC: already dashed (e.g. GDC manifest stem)
    ("anywhere/C3L-00016",                                    "CPTAC-3",   "C3L-00016"),
    ("anywhere/C3N-12345",                                    "CPTAC-3",   "C3N-12345"),

    # TCGA: dashless (encoded form from colleague's nnUNet dump)
    ("imagesTr/sub-TCGA020003_0000.nii.gz",                  "TCGA",      "TCGA-02-0003"),
    ("imagesTr/sub-TCGAHTA61C_0000.nii.gz",                  "TCGA",      "TCGA-HT-A61C"),
    # TCGA: already dashed
    ("preprocessed/tcga/TCGA-06-1084_t1c_preprocessed.nii.gz", "TCGA",    "TCGA-06-1084"),
    ("anywhere/TCGA-AB-1234.something",                       "TCGA",     "TCGA-AB-1234"),

    # UPENN: both encoded variants
    ("preprocessed/upenn/UPenn_UPENN-GBM-00502_11_t1_preprocessed.nii.gz",
     "UPENN-GBM", "UPENN-GBM-00502"),
    ("data/UPENN/normalize/sub-UPENNGBM00001_ses-6_ce-GD_run-1_T1w_n4_register_ss_normalize_t1.nii.gz",
     "UPENN-GBM", "UPENN-GBM-00001"),

    # UCSF-PDGM
    ("preprocessed/ucsf/UCSF-PDGM-0078_flair_preprocessed.nii.gz",
     "UCSF-PDGM", "UCSF-PDGM-0078"),
]


@pytest.mark.parametrize("path,project,submitter", ID_CASES)
def test_id_parsing(path, project, submitter):
    got_project, got_submitter, _ = parse_latent_path(path)
    assert got_project == project, f"project mismatch on {path}"
    assert got_submitter == submitter, f"submitter mismatch on {path}"


# (path, expected_modality)
MODALITY_CASES = [
    # explicit token in filename
    ("anywhere/sub-X_t1c_preprocessed.nii.gz",   "t1c"),
    ("anywhere/sub-X_t1_preprocessed.nii.gz",    "t1"),
    ("anywhere/sub-X_t2_preprocessed.nii.gz",    "t2"),
    ("anywhere/sub-X_flair_preprocessed.nii.gz", "flair"),
    # T1Gd / T1ce should normalize to t1c
    ("anywhere/sub-X_T1Gd_preprocessed.nii.gz",  "t1c"),
    ("anywhere/sub-X_T1ce_preprocessed.nii.gz",  "t1c"),
    # No explicit token — must NOT infer from _0000 channel suffix
    ("temp/cptac/sub-C3L00016_0000.nii.gz",      None),
    ("imagesTr/sub-TCGA020003_0000.nii.gz",      None),
]


@pytest.mark.parametrize("path,modality", MODALITY_CASES)
def test_modality_parsing(path, modality):
    _, _, got = parse_latent_path(path)
    assert got == modality, f"modality mismatch on {path}: got {got}, want {modality}"


# Paths that should NOT match any cohort.
NEGATIVE_CASES = [
    "preprocessed/adni/ADNI_002_S_0295_t1_preprocessed.nii.gz",
    "preprocessed/ppmi/PPMI_171162_t1_preprocessed.nii.gz",
    "preprocessed/hbn/sub-NDARAU939WUK_t1_preprocessed.nii.gz",
    "preprocessed/ixi/IXI565-HH-2534_t1_preprocessed.nii.gz",
    "/random/path/no-cohort-id.nii.gz",
]


@pytest.mark.parametrize("path", NEGATIVE_CASES)
def test_negative_paths_return_none(path):
    project, submitter, _ = parse_latent_path(path)
    assert project is None, f"unexpected match {project} on {path}"
    assert submitter is None
