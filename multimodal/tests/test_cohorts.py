"""Smoke tests for cohort ID extraction.

These hit the real local metadata files — they verify shape and ID format,
not exhaustive content.
"""
import pytest

from multimodal import cohorts


def test_config_loads():
    cfg = cohorts.load_cohorts_config()
    assert "cohorts" in cfg
    assert {"tcga_gbm", "tcga_lgg", "cptac_gbm", "upenn_gbm", "ucsf_pdgm"} <= set(cfg["cohorts"])


def test_load_imaging_subjects_returns_long_format():
    df = cohorts.load_imaging_subjects()
    if df.empty:
        pytest.skip("no local metadata; nothing to test")
    assert {"cohort", "imaging_id", "submitter_id", "project_id"} <= set(df.columns)


def test_tcga_ids_have_correct_prefix():
    df = cohorts.load_imaging_subjects()
    if df.empty:
        pytest.skip("no local metadata")
    tcga = df[df["cohort"].str.startswith("tcga")]
    if len(tcga):
        assert tcga["submitter_id"].str.startswith("TCGA-").all()


def test_cptac_ids_have_correct_prefix():
    df = cohorts.load_imaging_subjects()
    if df.empty:
        pytest.skip("no local metadata")
    cptac = df[df["cohort"] == "cptac_gbm"]
    if len(cptac):
        assert cptac["submitter_id"].str.match(r"C3[LN]-\d+").all()
