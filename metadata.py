"""Build harmonized clinical/phenotypic metadata tables per dataset.

Each extractor returns a DataFrame with a canonical column set (see SCHEMA)
so all datasets can be concatenated into one master sheet later.

Usage:
    python metadata.py adni --loni data/metadata/ADNI/LONI --output data/metadata/ADNI/adni_master.csv
"""

import argparse
import logging
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


log = logging.getLogger(__name__)


# Canonical columns for the master sheet. Per-scan granularity.
SCHEMA = [
    "dataset",
    "subject_id",
    "scan_id",
    "visit_code",
    "scan_date",
    "age_at_scan",
    "sex",              # M / F
    "handedness",       # R / L / A
    "education_years",
    "ethnicity",
    "race",
    "site",
    "dx",               # harmonized: CN / MCI / AD / ... / HC / ASD / ADHD / SCZ / PD / GBM / ...
    "dx_detail",        # dataset-specific label (e.g. "EMCI", "LMCI", "Glioblastoma IDH-wt")
    "mmse",             # ADNI + OASIS
    "cdr_global",       # ADNI + OASIS
    "moca",             # ADNI (ADNIGO+) + PPMI
    "mds_updrs_iii",    # PPMI (motor severity, 0-132)
    "gds",              # ADNI + PPMI (depression, 0-15)
    "scopa_aut",        # PPMI (autonomic symptoms)
    "apoe",             # "E3/E3" / "E3/E4" / ...
    "tumor_grade",      # WHO CNS grade 2/3/4 (UCSF) — gliomas only
    "idh_status",       # "wildtype" / "mutant" (UCSF + UPENN)
    "mgmt_status",      # "methylated" / "unmethylated" (UCSF + UPENN)
    "scanner_manufacturer",
    "scanner_model",
    "field_strength",
    "modality",         # T1 / T1c / T2 / FLAIR / ...
]


# =============================================================================
# ADNI
# =============================================================================

# ADNI DXSUM: DIAGNOSIS column codes
# 1 = CN (Normal), 2 = MCI, 3 = AD
ADNI_DX_MAP = {1: "CN", 2: "MCI", 3: "AD"}

# PTGENDER: 1 = Male, 2 = Female
ADNI_SEX_MAP = {1: "M", 2: "F"}

# PTHAND: 1 = Right, 2 = Left
ADNI_HAND_MAP = {1: "R", 2: "L"}

# PTETHCAT: 1 = Hisp/Latino, 2 = Not Hisp/Latino, 3 = Unknown
ADNI_ETHNICITY_MAP = {1: "Hispanic/Latino", 2: "Not Hispanic/Latino", 3: "Unknown"}

# PTRACCAT: 1 = Am Indian/Alaska Native, 2 = Asian, 3 = Native Hawaiian/Pac Isl,
#           4 = Black, 5 = White, 6 = More than one, 7 = Unknown
ADNI_RACE_MAP = {
    "1": "American Indian/Alaska Native",
    "2": "Asian",
    "3": "Native Hawaiian/Pacific Islander",
    "4": "Black",
    "5": "White",
    "6": "More than one",
    "7": "Unknown",
}


def _read_csv(path: Path, **kwargs) -> pd.DataFrame:
    log.info(f"Reading {path.name}")
    return pd.read_csv(path, low_memory=False, **kwargs)


def _nearest_visit_join(left: pd.DataFrame, right: pd.DataFrame,
                         keys: list[str], date_col_left: str,
                         date_col_right: str, tolerance_days: int = 180) -> pd.DataFrame:
    """Join left to right on keys, matching the right row with the closest date
    (within tolerance). Used for joining scans to clinical visits."""
    left = left.copy()
    right = right.copy()
    left[date_col_left] = pd.to_datetime(left[date_col_left], errors="coerce").astype("datetime64[ns]")
    right[date_col_right] = pd.to_datetime(right[date_col_right], errors="coerce").astype("datetime64[ns]")
    # merge_asof needs both sides sorted and non-null on the date column
    right = right.dropna(subset=[date_col_right])
    left_has_date = left[date_col_left].notna()
    left_sorted = left[left_has_date].sort_values(date_col_left)
    right_sorted = right.sort_values(date_col_right)
    merged = pd.merge_asof(
        left_sorted,
        right_sorted,
        left_on=date_col_left,
        right_on=date_col_right,
        by=keys,
        tolerance=pd.Timedelta(days=tolerance_days),
        direction="nearest",
    )
    # Re-attach rows that had null dates on the left (with no right-side match)
    if (~left_has_date).any():
        merged = pd.concat([merged, left[~left_has_date]], ignore_index=True)
    return merged


def _extract_adni_xml(xml_dir: Path) -> pd.DataFrame:
    """Parse per-scan ADNI IDA XMLs (one per preprocessed scan).

    Each XML has: researchGroup, sex, age, APOE, MMSE/CDR/GDS/FAQ/NPI at scan
    visit, scanner protocol. Richer than LONI CSVs for our preprocessed subset.

    Note: the XMLs use `<project xmlns="">` which overrides the outer namespace,
    so <subject> and descendants are in the empty namespace (no prefix).
    """
    xml_files = sorted(xml_dir.glob("ADNI_*.xml"))
    log.info(f"Parsing {len(xml_files)} ADNI XMLs from {xml_dir}")

    rows = []
    for xml_path in xml_files:
        try:
            root = ET.parse(xml_path).getroot()
        except ET.ParseError as e:
            log.warning(f"Parse error {xml_path.name}: {e}")
            continue

        project = root.find(".//project")
        if project is None:
            continue
        subj = project.find("subject")
        if subj is None:
            continue

        # Filename looks like ADNI_{subject}_{pipeline}_S{series}_I{imageUID}.xml
        # so we always have a fallback for scan_id.
        fname_stem = xml_path.stem
        scan_id_fallback = fname_stem.rsplit("_I", 1)[-1] if "_I" in fname_stem else ""

        row = {
            "subject_id": (subj.findtext("subjectIdentifier") or "").strip(),
            "research_group": (subj.findtext("researchGroup") or "").strip(),
            "sex": (subj.findtext("subjectSex") or "").strip(),
            "site": (project.findtext("siteKey") or "").strip(),
            "scan_id": scan_id_fallback,
        }
        # APOE alleles
        for info in subj.findall("subjectInfo"):
            item = info.get("item", "")
            if item == "APOE A1":
                row["apoe_a1"] = info.text
            elif item == "APOE A2":
                row["apoe_a2"] = info.text
        # Visit + assessment scores
        visit = subj.find("visit")
        if visit is not None:
            row["visit_code"] = visit.findtext("visitIdentifier") or ""
            for assessment in visit.findall("assessment"):
                name = assessment.get("name", "")
                score_el = assessment.find(".//assessmentScore")
                if score_el is None or score_el.text is None:
                    continue
                try:
                    score = float(score_el.text)
                except ValueError:
                    continue
                if name == "MMSE":
                    row["mmse"] = score
                elif name == "CDR":
                    row["cdr_global"] = score
                elif name == "GDSCALE":
                    row["gds"] = score
                elif name == "Functional Assessment Questionnaire":
                    row["faq"] = score
                elif name.startswith("Neuropsychiatric"):
                    row["npiq"] = score

        # Scan / study info
        study = subj.find("study")
        if study is not None:
            row["age_at_scan"] = study.findtext("subjectAge")
            row["weight_kg"] = study.findtext("weightKg")
            series = study.find("series")
            if series is not None:
                row["scan_date"] = series.findtext("dateAcquired")
                row["series_id"] = series.findtext("seriesIdentifier")
                # Image UID from the derived product (overrides filename fallback)
                derived = series.find(".//derivedProduct")
                if derived is not None:
                    uid = derived.findtext("imageUID")
                    if uid:
                        row["scan_id"] = uid
                # Scanner protocol (from original related image)
                for proto in series.findall(".//protocolTerm/protocol"):
                    term = proto.get("term", "")
                    if term == "Manufacturer":
                        row["scanner_manufacturer"] = proto.text
                    elif term == "Mfg Model":
                        row["scanner_model"] = proto.text
                    elif term == "Field Strength":
                        row["field_strength"] = proto.text
                    elif term == "Weighting":
                        row["modality"] = proto.text  # T1 / T2 / ...

        rows.append(row)

    df = pd.DataFrame(rows)
    log.info(f"Parsed {len(df)} XML records")

    # Harmonize
    df["dataset"] = "ADNI"
    # Ensure scan_id has the I prefix (filename fallback already does; derivedProduct doesn't)
    df["scan_id"] = df["scan_id"].astype(str)
    df["scan_id"] = df["scan_id"].where(df["scan_id"].str.startswith("I"), "I" + df["scan_id"])
    df["apoe"] = df.get("apoe_a1", "").fillna("") + df.get("apoe_a2", "").fillna("")
    df["dx_detail"] = df["research_group"]  # CN/EMCI/LMCI/AD/SMC
    df["dx"] = df["research_group"].map({
        "CN": "CN", "SMC": "CN",
        "EMCI": "MCI", "LMCI": "MCI", "MCI": "MCI",
        "AD": "AD",
    })
    df["age_at_scan"] = pd.to_numeric(df["age_at_scan"], errors="coerce").round(1)
    df["field_strength"] = pd.to_numeric(df["field_strength"], errors="coerce")

    # Merge with LONI demographics for education/handedness/ethnicity/race
    # (not in XML). Caller can do this; for standalone XML output leave NA.
    for col in ["handedness", "education_years", "ethnicity", "race", "moca",
                "mds_updrs_iii", "scopa_aut"]:
        df[col] = pd.NA
    # GDS is parsed from XML assessment scores already
    if "gds" not in df.columns:
        df["gds"] = pd.NA

    return df[SCHEMA]


def extract_adni(loni_dir: Path, xml_dir: Path) -> pd.DataFrame:
    """Build ADNI master: start from per-scan XMLs (one per preprocessed scan,
    with research group, APOE, and visit-specific MMSE/CDR/GDS/FAQ/NPI),
    enrich with LONI demographics (education, handedness, ethnicity, race),
    then fill MMSE/CDR/MoCA gaps from the longitudinal LONI assessment CSVs."""
    xml_df = _extract_adni_xml(xml_dir)

    # --- Demographics (static per subject) ---
    ptdemog = _read_csv(next(loni_dir.glob("PTDEMOG_*.csv")))
    ptdemog_static = ptdemog.sort_values("VISDATE").groupby("PTID").first().reset_index()
    ptdemog_static["handedness"] = ptdemog_static["PTHAND"].map(ADNI_HAND_MAP)
    ptdemog_static["ethnicity_str"] = ptdemog_static["PTETHCAT"].map(ADNI_ETHNICITY_MAP)
    ptdemog_static["race_str"] = ptdemog_static["PTRACCAT"].map(ADNI_RACE_MAP)
    demog = ptdemog_static[["PTID", "handedness", "PTEDUCAT", "ethnicity_str", "race_str"]].rename(
        columns={"PTID": "subject_id", "PTEDUCAT": "education_years",
                 "ethnicity_str": "ethnicity", "race_str": "race"}
    )
    xml_df = xml_df.drop(columns=["handedness", "education_years", "ethnicity", "race"])
    xml_df = xml_df.merge(demog, on="subject_id", how="left")

    # --- MMSE from LONI longitudinal (prefer over XML value) ---
    xml_df["scan_date"] = pd.to_datetime(xml_df["scan_date"], errors="coerce")
    mmse = _read_csv(next(loni_dir.glob("MMSE_*.csv")))
    mmse_df = mmse[["PTID", "VISDATE", "MMSCORE"]].dropna(subset=["MMSCORE"])
    mmse_df = mmse_df[mmse_df["MMSCORE"] >= 0].rename(
        columns={"PTID": "subject_id", "VISDATE": "mmse_date", "MMSCORE": "mmse_loni"}
    )
    xml_df = _nearest_visit_join(xml_df, mmse_df, keys=["subject_id"],
                                  date_col_left="scan_date", date_col_right="mmse_date",
                                  tolerance_days=180)
    xml_df["mmse"] = xml_df["mmse_loni"].fillna(xml_df["mmse"])

    # --- CDR from LONI longitudinal ---
    cdr = _read_csv(next(loni_dir.glob("CDR_*.csv")))
    cdr_df = cdr[["PTID", "VISDATE", "CDGLOBAL"]].dropna(subset=["CDGLOBAL"])
    cdr_df = cdr_df[cdr_df["CDGLOBAL"] >= 0].rename(
        columns={"PTID": "subject_id", "VISDATE": "cdr_date", "CDGLOBAL": "cdr_loni"}
    )
    xml_df = _nearest_visit_join(xml_df, cdr_df, keys=["subject_id"],
                                  date_col_left="scan_date", date_col_right="cdr_date",
                                  tolerance_days=180)
    xml_df["cdr_global"] = xml_df["cdr_loni"].fillna(xml_df["cdr_global"])

    # --- MoCA from LONI longitudinal (XML doesn't have it) ---
    try:
        moca = _read_csv(next(loni_dir.glob("MOCA_*.csv")))
        if "MOCA" in moca.columns:
            moca_df = moca[["PTID", "VISDATE", "MOCA"]].dropna(subset=["MOCA"]).rename(
                columns={"PTID": "subject_id", "VISDATE": "moca_date", "MOCA": "moca_loni"}
            )
            xml_df = _nearest_visit_join(xml_df, moca_df, keys=["subject_id"],
                                          date_col_left="scan_date", date_col_right="moca_date",
                                          tolerance_days=180)
            xml_df["moca"] = xml_df["moca_loni"]
    except (StopIteration, KeyError):
        pass

    # Sanitize ADNI sentinel values (-1 = "not done / missing")
    for col in ("mmse", "cdr_global", "moca"):
        xml_df.loc[xml_df[col] < 0, col] = pd.NA

    # Re-stringify scan_date for output
    xml_df["scan_date"] = xml_df["scan_date"].dt.strftime("%Y-%m-%d")

    return xml_df[SCHEMA]


# =============================================================================
# PPMI
# =============================================================================

# PPMI Demographics: SEX 0 = Female, 1 = Male (per PPMI data dictionary)
PPMI_SEX_MAP = {0: "F", 1: "M"}

# PPMI Demographics: HANDED 1 = Right, 2 = Left, 3 = Mixed
PPMI_HAND_MAP = {1: "R", 2: "L", 3: "A"}

# PPMI COHORT code → harmonized dx
PPMI_COHORT_DX_MAP = {
    "Parkinson's Disease": "PD",
    "Healthy Control": "HC",
    "Prodromal": "Prodromal",
    "SWEDD": "SWEDD",
}

# GDS-15 positive items (reverse-coded: "No"=1 counts as depressed)
GDS_POSITIVE_ITEMS = ["GDSSATIS", "GDSGSPIR", "GDSHAPPY", "GDSALIVE", "GDSENRGY"]
GDS_NEGATIVE_ITEMS = ["GDSDROPD", "GDSEMPTY", "GDSBORED", "GDSAFRAD", "GDSHLPLS",
                      "GDSHOME", "GDSMEMRY", "GDSWRTLS", "GDSHOPLS", "GDSBETER"]


def _ppmi_race(row) -> str:
    """Collapse PPMI boolean race columns into a single label."""
    flags = {
        "White": row.get("RAWHITE"),
        "Black": row.get("RABLACK"),
        "Asian": row.get("RAASIAN"),
        "Native Hawaiian/Pacific Islander": row.get("RAHAWOPI"),
        "American Indian/Alaska Native": row.get("RAINDALS"),
        "Other": row.get("RANOS"),
    }
    positive = [k for k, v in flags.items() if v == 1]
    if len(positive) == 0:
        if row.get("RAUNKNOWN") == 1:
            return "Unknown"
        return pd.NA
    if len(positive) == 1:
        return positive[0]
    return "More than one"


def _compute_gds_total(gds: pd.DataFrame) -> pd.Series:
    """GDS-15 total: sum of negative items as-is + reverse-coded positive items."""
    neg = gds[GDS_NEGATIVE_ITEMS].sum(axis=1, min_count=10)
    pos_reversed = (1 - gds[GDS_POSITIVE_ITEMS]).sum(axis=1, min_count=5)
    return neg + pos_reversed


def extract_ppmi(loni_dir: Path) -> pd.DataFrame:
    """Build PPMI master from LONI CSVs (per-scan from Xing Core Lab MRI metadata).

    Expects in `loni_dir`:
        Xing_Core_Lab_-_MRI_acquisition_metadata_*.csv  — seed (per-scan)
        Participant_Status_*.csv                         — cohort
        Demographics_*.csv                               — sex, race, handedness, dob
        Age_at_visit_*.csv                               — exact age per visit
        Primary_Research_Diagnosis_*.csv                 — research dx per visit
        MDS-UPDRS_Part_III_*.csv                         — NP3TOT motor score
        Montreal_Cognitive_Assessment__MoCA__*.csv       — MCATOT
        Geriatric_Depression_Scale__Short_Version__*.csv — GDS items
        SCOPA-AUT_*.csv                                  — autonomic symptoms
        iu_genetic_consensus_*.csv                       — APOE + PD risk variants
    """
    # --- Seed: Participant_Status (ALL 8,417 PPMI subjects) ---
    ps = _read_csv(next(loni_dir.glob("Participant_Status_*.csv")))
    df = ps[["PATNO", "COHORT_DEFINITION", "ENROLL_AGE", "ENROLL_DATE"]].rename(
        columns={"COHORT_DEFINITION": "cohort"}
    )
    log.info(f"PPMI: seeding on Participant_Status, {len(df)} subjects")

    # --- Demographics (static per subject) ---
    demog = _read_csv(next(loni_dir.glob("Demographics_*.csv")))
    demog_static = demog.sort_values("INFODT").groupby("PATNO").first().reset_index()
    demog_static["sex"] = demog_static["SEX"].map(PPMI_SEX_MAP)
    demog_static["handedness"] = demog_static["HANDED"].map(PPMI_HAND_MAP)
    demog_static["race"] = demog_static.apply(_ppmi_race, axis=1)
    demog_static["ethnicity"] = demog_static["HISPLAT"].map(
        {1: "Hispanic/Latino", 0: "Not Hispanic/Latino"}
    )
    demog_static["BIRTHDT"] = pd.to_datetime(demog_static["BIRTHDT"], errors="coerce", format="%m/%Y")
    df = df.merge(
        demog_static[["PATNO", "sex", "handedness", "race", "ethnicity", "BIRTHDT"]],
        on="PATNO", how="left"
    )

    # --- Education (Socio-Economics, static per subject) ---
    socio_files = list(loni_dir.glob("Socio-Economics_*.csv"))
    if socio_files:
        socio = _read_csv(socio_files[0])
        socio_static = socio.dropna(subset=["EDUCYRS"]).sort_values("INFODT") \
            .groupby("PATNO").last().reset_index()
        df = df.merge(
            socio_static[["PATNO", "EDUCYRS"]].rename(columns={"EDUCYRS": "education_years"}),
            on="PATNO", how="left"
        )

    # --- Baseline clinical (take EVENT_ID=BL or earliest available visit) ---
    def _baseline_per_subject(path: Path, col: str, out_name: str) -> pd.DataFrame:
        """Take the BL (baseline) row per subject, falling back to earliest INFODT."""
        raw = _read_csv(path)
        raw = raw.dropna(subset=[col])
        if raw.empty:
            return pd.DataFrame(columns=["PATNO", out_name])
        # Prefer BL
        bl = raw[raw["EVENT_ID"] == "BL"]
        rest = raw[raw["EVENT_ID"] != "BL"].sort_values("INFODT").groupby("PATNO").first().reset_index()
        merged = pd.concat([bl, rest]).drop_duplicates("PATNO", keep="first")
        return merged[["PATNO", col]].rename(columns={col: out_name})

    # MDS-UPDRS III (motor) — baseline
    updrs_path = None
    for p in [loni_dir / "Motor_MDS_UPDRS" / "MDS-UPDRS_Part_III_11Apr2026.csv",
              *loni_dir.glob("MDS-UPDRS_Part_III_*.csv")]:
        if p.exists():
            updrs_path = p
            break
    if updrs_path:
        updrs_bl = _baseline_per_subject(updrs_path, "NP3TOT", "mds_updrs_iii")
        df = df.merge(updrs_bl, on="PATNO", how="left")

    # MoCA — baseline
    moca_path = next(loni_dir.glob("Montreal_Cognitive_Assessment*.csv"))
    moca_bl = _baseline_per_subject(moca_path, "MCATOT", "moca")
    df = df.merge(moca_bl, on="PATNO", how="left")

    # GDS-15 — baseline (compute total first)
    gds_raw = _read_csv(next(loni_dir.glob("Geriatric_Depression_Scale*.csv")))
    gds_raw["gds_total"] = _compute_gds_total(gds_raw)
    gds_clean = gds_raw[["PATNO", "EVENT_ID", "INFODT", "gds_total"]].dropna(subset=["gds_total"])
    bl = gds_clean[gds_clean["EVENT_ID"] == "BL"]
    rest = gds_clean[gds_clean["EVENT_ID"] != "BL"].sort_values("INFODT").groupby("PATNO").first().reset_index()
    gds_bl = pd.concat([bl, rest]).drop_duplicates("PATNO", keep="first")
    df = df.merge(gds_bl[["PATNO", "gds_total"]].rename(columns={"gds_total": "gds"}),
                   on="PATNO", how="left")

    # SCOPA-AUT — baseline
    scopa_raw = _read_csv(next(loni_dir.glob("SCOPA-AUT_*.csv")))
    scopa_items = [c for c in scopa_raw.columns if c.startswith("SCAU") and pd.api.types.is_numeric_dtype(scopa_raw[c])]
    scopa_raw["scopa_sum"] = scopa_raw[scopa_items].sum(axis=1, min_count=10)
    scopa_clean = scopa_raw[["PATNO", "EVENT_ID", "INFODT", "scopa_sum"]].dropna(subset=["scopa_sum"])
    bl = scopa_clean[scopa_clean["EVENT_ID"] == "BL"]
    rest = scopa_clean[scopa_clean["EVENT_ID"] != "BL"].sort_values("INFODT").groupby("PATNO").first().reset_index()
    scopa_bl = pd.concat([bl, rest]).drop_duplicates("PATNO", keep="first")
    df = df.merge(scopa_bl[["PATNO", "scopa_sum"]].rename(columns={"scopa_sum": "scopa_aut"}),
                   on="PATNO", how="left")

    # --- Genetic consensus (static) — APOE + PD risk variants ---
    gen_files = list(loni_dir.glob("iu_genetic_consensus_*.csv"))
    if gen_files:
        gen = _read_csv(gen_files[0])
        gen_df = gen[["PATNO", "APOE", "LRRK2", "GBA", "VPS35", "SNCA", "PRKN", "PARK7", "PINK1"]].rename(
            columns={"APOE": "apoe"}
        )
        df = df.merge(gen_df, on="PATNO", how="left")

    # Age: use ENROLL_AGE from Participant_Status (per-subject, baseline age)
    df["age_at_scan"] = pd.to_numeric(df["ENROLL_AGE"], errors="coerce").round(1)

    # --- Canonical output (one row per subject) ---
    out = pd.DataFrame({
        "dataset": "PPMI",
        "subject_id": df["PATNO"].astype(str),
        "scan_id": df["PATNO"].astype(str),  # subject-level: scan_id == subject_id
        "visit_code": "baseline",
        "scan_date": pd.NA,  # no per-scan date (use file's year-month for linking)
        "age_at_scan": df["age_at_scan"],
        "sex": df["sex"],
        "handedness": df["handedness"],
        "education_years": df.get("education_years"),
        "ethnicity": df["ethnicity"],
        "race": df["race"],
        "site": df["PATNO"].astype(str).str[:2],  # PPMI site from first 2 digits of PATNO
        "dx": df["cohort"].map(PPMI_COHORT_DX_MAP),
        "dx_detail": df["cohort"],
        "mmse": pd.NA,  # not in PPMI battery
        "cdr_global": pd.NA,  # not in PPMI battery
        "moca": df["moca"],
        "mds_updrs_iii": df.get("mds_updrs_iii"),
        "gds": df.get("gds"),
        "scopa_aut": df.get("scopa_aut"),
        "apoe": df.get("apoe"),
        "scanner_manufacturer": pd.NA,
        "scanner_model": pd.NA,
        "field_strength": pd.NA,
        "modality": "T1",
    })
    return out[SCHEMA]


# =============================================================================
# Tumor datasets (subject-level — one imaging session per subject)
# =============================================================================


def _na_row(**overrides) -> dict:
    """Template row with every SCHEMA column set to NA, then overridden."""
    row = {col: pd.NA for col in SCHEMA}
    row.update(overrides)
    return row


def _cbioportal_to_wide(json_path: Path) -> pd.DataFrame:
    """Pivot cBioPortal long-format clinical JSON into wide DataFrame (one row per patient)."""
    import json
    with open(json_path) as f:
        records = json.load(f)
    long = pd.DataFrame(records)[["patientId", "clinicalAttributeId", "value"]]
    wide = long.pivot_table(index="patientId", columns="clinicalAttributeId",
                              values="value", aggfunc="first").reset_index()
    return wide


def extract_tcga(meta_dir: Path) -> pd.DataFrame:
    """TCGA-GBM/LGG: one row per subject. Merges local tumor genetics with
    demographics pulled from cBioPortal (gbm_tcga + lgg_tcga clinical JSON)."""
    local = _read_csv(meta_dir / "metadata" / "TCGA_metadata.csv")
    grade_to_dx = {"TCGA-GBM": "GBM", "TCGA-LGG": "LGG"}

    # Merge cBioPortal demographics (AGE, SEX, RACE, ETHNICITY, KPS, OS)
    cbio_dir = meta_dir / "cbioportal"
    cbio_parts = []
    for name in ("gbm_clinical.json", "lgg_clinical.json"):
        p = cbio_dir / name
        if p.exists():
            cbio_parts.append(_cbioportal_to_wide(p))
    if cbio_parts:
        cbio = pd.concat(cbio_parts, ignore_index=True).rename(columns={"patientId": "eid"})
        df = local.merge(cbio, on="eid", how="left")
    else:
        df = local

    sex_map = {"Male": "M", "Female": "F", "MALE": "M", "FEMALE": "F"}
    out = pd.DataFrame([
        _na_row(
            dataset="TCGA",
            subject_id=row["eid"],
            scan_id=row["eid"],
            age_at_scan=pd.to_numeric(row.get("AGE"), errors="coerce"),
            sex=sex_map.get(str(row.get("SEX", "")).strip()),
            ethnicity=str(row.get("ETHNICITY", "")).strip() or pd.NA,
            race=str(row.get("RACE", "")).strip().title() or pd.NA,
            dx=grade_to_dx.get(row["grade"], "Glioma"),
            dx_detail=f"{row['grade']} | IDH1={row['IDH1']} | IDH2={row['IDH2']} | ATRX={'mut' if row['ATRX_bin'] else 'wt'} | CDKN2A={'mut' if row['CDKN2A_bin'] else 'wt'}",
            modality="T1,T1c,T2,FLAIR",
        )
        for _, row in df.iterrows()
    ])
    # Clean up empty-string / "Nan" → NA from the string casts above
    for col in ("ethnicity", "race"):
        out.loc[out[col].isin(["", "nan", "Nan", "NaN", "None"]), col] = pd.NA
    return out[SCHEMA]


def extract_upenn(meta_dir: Path) -> pd.DataFrame:
    """UPENN-GBM: one row per subject, all grade IV glioblastoma.
    Emits IDH + MGMT molecular markers (no grade column — all GBM)."""
    df = _read_csv(meta_dir / "UPENN-GBM_clinical_info_v2.1.csv")

    def _idh(s):
        if pd.isna(s):
            return pd.NA
        s = str(s).strip().lower()
        if s == "wildtype":
            return "wildtype"
        if s == "mutated":
            return "mutant"
        return pd.NA  # "NOS/NEC" etc.

    def _mgmt(s):
        if pd.isna(s):
            return pd.NA
        s = str(s).strip().lower()
        if s == "methylated":
            return "methylated"
        if s == "unmethylated":
            return "unmethylated"
        return pd.NA  # "Not Available", "Indeterminate"

    out = pd.DataFrame([
        _na_row(
            dataset="UPENN",
            subject_id=row["ID"],
            scan_id=row["ID"],
            age_at_scan=pd.to_numeric(row["Age_at_scan_years"], errors="coerce"),
            sex={"M": "M", "F": "F"}.get(str(row["Gender"]).strip()),
            dx="GBM",
            dx_detail=f"IDH1={row['IDH1']} | MGMT={row['MGMT']} | KPS={row['KPS']}",
            tumor_grade=4,  # UPENN-GBM is by definition grade IV
            idh_status=_idh(row["IDH1"]),
            mgmt_status=_mgmt(row["MGMT"]),
            modality="T1,T1c,T2,FLAIR",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


def extract_ucsf(meta_dir: Path) -> pd.DataFrame:
    """UCSF-PDGM: one row per subject, split by WHO CNS grade.
    Grade IV → GBM; Grade II/III → Glioma (LGG + grade III astrocytoma).
    Also emits molecular markers: tumor_grade, idh_status, mgmt_status."""
    df = _read_csv(meta_dir / "UCSF-PDGM-metadata_v2.csv")

    def _dx_from_grade(grade) -> str:
        try:
            g = int(grade)
        except (ValueError, TypeError):
            return "Glioma"
        return "GBM" if g == 4 else "Glioma"

    def _grade(g):
        try:
            return int(g)
        except (ValueError, TypeError):
            return pd.NA

    def _idh(s):
        if pd.isna(s):
            return pd.NA
        s = str(s).lower()
        if "wildtype" in s:
            return "wildtype"
        if "mut" in s or "r132" in s or "r172" in s:
            return "mutant"
        return pd.NA

    def _mgmt(s):
        if pd.isna(s):
            return pd.NA
        s = str(s).lower().strip()
        if s == "positive" or "methylated" in s and "unmeth" not in s:
            return "methylated"
        if s == "negative" or "unmethylated" in s:
            return "unmethylated"
        return pd.NA

    out = pd.DataFrame([
        _na_row(
            dataset="UCSF",
            subject_id=row["ID"],
            scan_id=row["ID"],
            age_at_scan=pd.to_numeric(row["Age at MRI"], errors="coerce"),
            sex=str(row["Sex"]).strip() if pd.notna(row["Sex"]) else pd.NA,
            dx=_dx_from_grade(row["WHO CNS Grade"]),
            dx_detail=f"WHO={row['WHO CNS Grade']} | {row['Final pathologic diagnosis (WHO 2021)']} | MGMT={row['MGMT status']} | 1p19q={row['1p/19q']} | IDH={row['IDH']}",
            tumor_grade=_grade(row["WHO CNS Grade"]),
            idh_status=_idh(row["IDH"]),
            mgmt_status=_mgmt(row["MGMT status"]),
            modality="T1,T1c,T2,FLAIR",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# SCHIZO (COBRE + MCIC merged from SchizConnect)
# =============================================================================


def extract_schizo(meta_dir: Path) -> pd.DataFrame:
    """SCHIZO: one row per subject, study + age + sex + dx."""
    df = pd.read_csv(meta_dir / "SCHIZO_metadata.tsv", sep="\t")
    out = pd.DataFrame([
        _na_row(
            dataset="SCHIZO",
            subject_id=row["eid"],
            scan_id=row["eid"],
            age_at_scan=pd.to_numeric(row["age"], errors="coerce"),
            sex={"male": "M", "female": "F"}.get(str(row["sex"]).lower().strip()),
            site=row["study"],  # MCICShare / COBRE
            dx={
                "Schizophrenia_Broad": "SCZ",
                "Schizophrenia_Strict": "SCZ",
                "No_Known_Disorder": "HC",
                "Bipolar_Disorder_Broad": "BP",
                "Bipolar_Disorder_Strict": "BP",
            }.get(str(row["dx"]).strip(), row["dx"]),
            dx_detail=row["dx"],
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# IXI (healthy lifespan — 3 London sites)
# =============================================================================


def extract_ixi(meta_dir: Path) -> pd.DataFrame:
    """IXI: 619 healthy subjects, age/sex/height/weight."""
    df = pd.read_excel(meta_dir / "IXI.xls")
    df["subject_id"] = "IXI" + df["IXI_ID"].astype(str).str.zfill(3)
    out = pd.DataFrame([
        _na_row(
            dataset="IXI",
            subject_id=row["subject_id"],
            scan_id=row["subject_id"],
            scan_date=row["STUDY_DATE"].strftime("%Y-%m-%d") if pd.notna(row["STUDY_DATE"]) else pd.NA,
            age_at_scan=round(row["AGE"], 1) if pd.notna(row["AGE"]) else pd.NA,
            sex={1: "M", 2: "F"}.get(row["SEX_ID (1=m, 2=f)"]),
            dx="HC",
            dx_detail="Healthy",
            modality="T1,T2",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# OASIS-1 (cross-sectional) and OASIS-2 (longitudinal)
# =============================================================================


def extract_oasis1(meta_dir: Path) -> pd.DataFrame:
    """OASIS-1: one row per MR session (OAS1_XXXX_MR1), CDR + MMSE."""
    xlsx_files = list(meta_dir.glob("oasis_cross-sectional-*.xlsx"))
    df = pd.read_excel(xlsx_files[0])
    cdr_to_dx = {0.0: "CN", 0.5: "MCI", 1.0: "AD", 2.0: "AD"}
    out = pd.DataFrame([
        _na_row(
            dataset="OASIS-1",
            subject_id=row["ID"].rsplit("_", 1)[0],  # OAS1_0001
            scan_id=row["ID"],                       # OAS1_0001_MR1
            visit_code="MR1",
            age_at_scan=row["Age"],
            sex={"M": "M", "F": "F"}.get(row["M/F"]),
            handedness=row["Hand"],
            education_years=row["Educ"],
            dx=cdr_to_dx.get(row["CDR"], "CN") if pd.notna(row["CDR"]) else "CN",
            dx_detail=f"CDR={row['CDR']}" if pd.notna(row["CDR"]) else "CDR=unk",
            mmse=row["MMSE"],
            cdr_global=row["CDR"],
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


def extract_oasis2(meta_dir: Path) -> pd.DataFrame:
    """OASIS-2: longitudinal, one row per MR session (OAS2_XXXX_MRn)."""
    xlsx_files = list(meta_dir.glob("oasis_longitudinal*.xlsx"))
    df = pd.read_excel(xlsx_files[0])
    group_to_dx = {
        "Nondemented": "CN",
        "Demented": "AD",
        "Converted": "MCI",  # CN → AD converters
    }
    out = pd.DataFrame([
        _na_row(
            dataset="OASIS-2",
            subject_id=row["Subject ID"],
            scan_id=row["MRI ID"],
            visit_code=f"MR{row['Visit']}",
            age_at_scan=row["Age"],
            sex={"M": "M", "F": "F"}.get(row["M/F"]),
            handedness=row["Hand"],
            education_years=row["EDUC"],
            dx=group_to_dx.get(row["Group"], "CN"),
            dx_detail=row["Group"],
            mmse=row["MMSE"],
            cdr_global=row["CDR"],
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# ABIDE I, ABIDE II (autism consortium)
# =============================================================================


def extract_abide(meta_dir: Path) -> pd.DataFrame:
    """ABIDE I: one row per subject, DX_GROUP 1=ASD 2=Control."""
    df = _read_csv(meta_dir / "Phenotypic_V1_0b.csv")
    out = pd.DataFrame([
        _na_row(
            dataset="ABIDE-I",
            subject_id=str(row["SUB_ID"]),
            scan_id=row["FILE_ID"] if pd.notna(row["FILE_ID"]) else str(row["SUB_ID"]),
            age_at_scan=row["AGE_AT_SCAN"],
            sex={1: "M", 2: "F"}.get(row["SEX"]),
            handedness={"R": "R", "L": "L", "Ambi": "A", "Mixed": "A"}.get(str(row["HANDEDNESS_CATEGORY"]).strip()),
            site=row["SITE_ID"],
            dx={1: "ASD", 2: "HC"}.get(row["DX_GROUP"]),
            dx_detail=str(row["DSM_IV_TR"]) if pd.notna(row["DSM_IV_TR"]) else pd.NA,
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


def extract_abide2(meta_dir: Path) -> pd.DataFrame:
    """ABIDE II: one row per subject (composite)."""
    df = _read_csv(meta_dir / "ABIDEII_Composite_Phenotypic.csv", encoding="latin-1")
    age_col = "AGE_AT_SCAN " if "AGE_AT_SCAN " in df.columns else "AGE_AT_SCAN"
    out = pd.DataFrame([
        _na_row(
            dataset="ABIDE-II",
            subject_id=str(row["SUB_ID"]),
            scan_id=str(row["SUB_ID"]),
            age_at_scan=row[age_col],
            sex={1: "M", 2: "F"}.get(row["SEX"]),
            handedness={"R": "R", "L": "L", "Ambi": "A", "Mixed": "A"}.get(str(row["HANDEDNESS_CATEGORY"]).strip()),
            site=row["SITE_ID"],
            dx={1: "ASD", 2: "HC"}.get(row["DX_GROUP"]),
            dx_detail=str(row.get("PDD_DSM_IV_TR", "")) or pd.NA,
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# ADHD-200
# =============================================================================


def extract_adhd200(meta_dir: Path) -> pd.DataFrame:
    """ADHD200: one row per subject, DX 0=Control 1=ADHD-C 2=ADHD-H 3=ADHD-I."""
    df = pd.read_csv(meta_dir / "adhd200_preprocessed_phenotypics.tsv", sep="\t")
    df["DX"] = pd.to_numeric(df["DX"], errors="coerce")  # "pending" → NaN
    dx_map = {0: "HC", 1: "ADHD", 2: "ADHD", 3: "ADHD"}
    dx_detail_map = {0: "Typically Developing", 1: "ADHD-Combined",
                     2: "ADHD-Hyperactive/Impulsive", 3: "ADHD-Inattentive"}
    out = pd.DataFrame([
        _na_row(
            dataset="ADHD200",
            subject_id=str(row["ScanDir ID"]).zfill(7),
            scan_id=str(row["ScanDir ID"]).zfill(7),
            age_at_scan=pd.to_numeric(row["Age"], errors="coerce"),
            sex={0: "F", 1: "M"}.get(row["Gender"]) if pd.notna(row["Gender"]) else pd.NA,
            handedness={1: "R", 0: "L"}.get(row["Handedness"]) if pd.notna(row["Handedness"]) else pd.NA,
            site=str(row["Site"]),
            dx=dx_map.get(int(row["DX"])) if pd.notna(row["DX"]) else pd.NA,
            dx_detail=dx_detail_map.get(int(row["DX"])) if pd.notna(row["DX"]) else pd.NA,
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# CoRR (consortium for reliability/reproducibility)
# =============================================================================


def extract_corr(meta_dir: Path) -> pd.DataFrame:
    """CoRR: subject + session, mixed healthy. SEX/HANDEDNESS are strings
    with '#' as missing sentinel."""
    df = _read_csv(meta_dir / "CoRR_AggregatedPhenotypicData.csv")
    sex_map = {"1": "M", "2": "F", 1: "M", 2: "F"}
    hand_map = {"R": "R", "L": "L", "A": "A",
                "Right": "R", "Left": "L", "Ambidextrous": "A"}
    out = pd.DataFrame([
        _na_row(
            dataset="CoRR",
            subject_id=str(row["SUBID"]),
            scan_id=f"{row['SUBID']}_{row['SESSION']}",
            visit_code=row["SESSION"],
            age_at_scan=pd.to_numeric(row["AGE_AT_SCAN_1"], errors="coerce"),
            sex=sex_map.get(str(row["SEX"]).strip()) if str(row["SEX"]).strip() != "#" else pd.NA,
            handedness=hand_map.get(str(row["HANDEDNESS"]).strip())
                if str(row["HANDEDNESS"]).strip() != "#" else pd.NA,
            site=row["SITE"],
            dx="HC",
            dx_detail="Healthy",
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# FCON1000 (33 sites, basic participants.tsv)
# =============================================================================


def extract_fcon1000(meta_dir: Path) -> pd.DataFrame:
    """FCON1000: concatenate 33 per-site BIDS participants.tsv files."""
    rows = []
    for tsv in sorted((meta_dir / "per_site").glob("*.tsv")):
        site = tsv.stem.replace("_participants", "")
        try:
            df = pd.read_csv(tsv, sep="\t")
        except Exception as e:
            log.warning(f"FCON1000 skip {tsv.name}: {e}")
            continue
        for _, row in df.iterrows():
            pid = row.get("participant_id") or row.get("Subject")
            if pd.isna(pid):
                continue
            age = pd.to_numeric(row.get("age"), errors="coerce")
            if pd.notna(age) and (age < 0 or age > 120):
                age = pd.NA  # 9999 / -1 / etc sentinel values
            sex_raw = str(row.get("sex", "")).strip().upper()
            sex = {"M": "M", "F": "F", "MALE": "M", "FEMALE": "F"}.get(sex_raw)
            hand = str(row.get("handed", row.get("handedness", ""))).strip().upper()
            hand = {"R": "R", "L": "L", "RIGHT": "R", "LEFT": "L", "A": "A"}.get(hand)
            rows.append(_na_row(
                dataset="FCON1000",
                subject_id=str(pid),
                scan_id=str(pid),
                age_at_scan=age,
                sex=sex,
                handedness=hand,
                site=site,
                dx="HC",
                dx_detail="Healthy",
                modality="T1",
            ))
    return pd.DataFrame(rows)[SCHEMA]


# =============================================================================
# HBN (Healthy Brain Network — pediatric psychiatric)
# =============================================================================


def extract_hbn(meta_dir: Path) -> pd.DataFrame:
    """HBN: 4 site-level BIDS participants.tsv (basic: age/sex/EHQ)."""
    rows = []
    for tsv in sorted((meta_dir / "per_site").glob("*.tsv")):
        site = tsv.stem.replace("_participants", "")  # Site-CBIC, Site-CUNY, etc
        # Some HBN files are actually CSV despite .tsv extension — auto-detect
        try:
            df = pd.read_csv(tsv, sep="\t")
            if "participant_id" not in df.columns:
                df = pd.read_csv(tsv, sep=",")
        except Exception as e:
            log.warning(f"HBN skip {tsv.name}: {e}")
            continue
        for _, row in df.iterrows():
            pid = row["participant_id"]
            age = pd.to_numeric(row.get("Age"), errors="coerce")
            sex_raw = row.get("Sex")
            sex = {0: "M", 1: "F", "0": "M", "1": "F", "M": "M", "F": "F"}.get(
                sex_raw if pd.notna(sex_raw) else None
            )
            ehq = pd.to_numeric(row.get("EHQ_Total"), errors="coerce")
            hand = pd.NA
            if pd.notna(ehq):
                hand = "R" if ehq > 40 else "L" if ehq < -40 else "A"
            rows.append(_na_row(
                dataset="HBN",
                subject_id=str(pid),
                scan_id=str(pid),
                age_at_scan=age,
                sex=sex,
                handedness=hand,
                site=site,
                dx="Pediatric-Mixed",
                dx_detail="HBN cohort (ADHD/ASD/anxiety/LD — full phenotype needs DUA)",
                modality="T1",
            ))
    return pd.DataFrame(rows)[SCHEMA]


# =============================================================================
# NKI, NKI2 (Rockland Sample)
# =============================================================================


def extract_nki(meta_dir: Path) -> pd.DataFrame:
    """NKI Rockland: basic participants.tsv."""
    df = pd.read_csv(meta_dir / "participants.tsv", sep="\t")
    out = pd.DataFrame([
        _na_row(
            dataset="NKI",
            subject_id=str(row["participant_id"]),
            scan_id=str(row["participant_id"]),
            age_at_scan=pd.to_numeric(row.get("age"), errors="coerce"),
            sex={"MALE": "M", "FEMALE": "F", "M": "M", "F": "F"}.get(
                str(row.get("sex", "")).strip().upper()
            ),
            handedness={"RIGHT": "R", "LEFT": "L", "AMBIDEXTROUS": "A"}.get(
                str(row.get("handedness", "")).strip().upper()
            ),
            site="NKI",
            dx="HC",
            dx_detail="NKI Rockland (full battery needs DUA)",
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


def extract_nki2(meta_dir: Path) -> pd.DataFrame:
    """NKI2: BIDS participants.tsv (ID-only). Backfills age/sex/handedness
    from NKI1 for the 114 overlapping subjects; the other ~142 stay blank
    until the Rockland DUA goes through."""
    df = pd.read_csv(meta_dir / "participants.tsv", sep="\t")
    df["subject_id"] = df["participant_id"].astype(str).str.replace("sub-", "", regex=False)

    # Backfill from NKI1 on subject_id
    nki1_path = meta_dir.parent / "NKI" / "participants.tsv"
    nki1 = pd.read_csv(nki1_path, sep="\t") if nki1_path.exists() else pd.DataFrame()
    if not nki1.empty:
        nki1 = nki1.rename(columns={"participant_id": "subject_id"})
        df = df.merge(nki1[["subject_id", "age", "sex", "handedness"]],
                       on="subject_id", how="left")

    sex_map = {"MALE": "M", "FEMALE": "F"}
    hand_map = {"RIGHT": "R", "LEFT": "L", "AMBIDEXTROUS": "A"}
    out = pd.DataFrame([
        _na_row(
            dataset="NKI2",
            subject_id=row["subject_id"],
            scan_id=row["subject_id"],
            age_at_scan=pd.to_numeric(row.get("age"), errors="coerce"),
            sex=sex_map.get(str(row.get("sex", "")).strip().upper()),
            handedness=hand_map.get(str(row.get("handedness", "")).strip().upper()),
            site="NKI",
            dx="HC",
            dx_detail="NKI Rockland 2 (demographics backfilled from NKI1 where available)",
            modality="T1",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# HCP-YA (Human Connectome Project, Young Adult)
# =============================================================================

# HCP age bins (exact age is Restricted — we get 5-year bins)
HCP_AGE_BIN_MID = {"22-25": 23.5, "26-30": 28, "31-35": 33, "36+": 37}


def extract_hcp_ya(meta_dir: Path) -> pd.DataFrame:
    """HCP-YA: 1,206 healthy subjects from BALSA Column Selector export.
    Age is a 5-year bin (Restricted otherwise). Harmonize cognition (MMSE,
    NIH Toolbox composites), NEO-FFI personality."""
    # Find the behavioral CSV (exported filename has timestamp)
    behav_files = list(meta_dir.glob("HCP_YA_subjects_*.csv"))
    if not behav_files:
        raise FileNotFoundError(f"No HCP_YA_subjects_*.csv in {meta_dir}")
    df = _read_csv(behav_files[0])

    out = pd.DataFrame([
        _na_row(
            dataset="HCP-YA",
            subject_id=str(row["Subject"]),
            scan_id=str(row["Subject"]),  # per-subject (1-2 T1s each)
            age_at_scan=HCP_AGE_BIN_MID.get(row["Age"], pd.NA),
            sex=row["Gender"],
            site="WashU",
            dx="HC",
            dx_detail=f"Age bin {row['Age']} | T1 count {row['T1_Count']}",
            mmse=row.get("MMSE_Score"),
            modality="T1,T2",
        )
        for _, row in df.iterrows()
    ])
    return out[SCHEMA]


# =============================================================================
# Concat: combine all per-dataset masters into a single big master
# =============================================================================


def concat_all(metadata_root: Path) -> pd.DataFrame:
    """Find every *_master.csv under metadata_root and concatenate them.
    All files must share the same SCHEMA."""
    csvs = sorted(metadata_root.rglob("*_master.csv"))
    if not csvs:
        raise FileNotFoundError(f"No *_master.csv found under {metadata_root}")

    parts = []
    for p in csvs:
        # Skip the combined output if it's inside the search tree
        if p.name == "all_datasets_master.csv":
            continue
        df = pd.read_csv(p, low_memory=False)
        log.info(f"  {p.parent.name:10s}: {len(df):>6d} rows")
        parts.append(df)

    combined = pd.concat(parts, ignore_index=True)
    # Enforce schema column order (missing cols get added as NA)
    for col in SCHEMA:
        if col not in combined.columns:
            combined[col] = pd.NA
    return combined[SCHEMA]


# =============================================================================
# CLI
# =============================================================================


# Map: dataset name → extractor function that takes a single Path
EXTRACTORS = {
    "tcga":     extract_tcga,
    "upenn":    extract_upenn,
    "ucsf":     extract_ucsf,
    "schizo":   extract_schizo,
    "ixi":      extract_ixi,
    "oasis1":   extract_oasis1,
    "oasis2":   extract_oasis2,
    "abide":    extract_abide,
    "abide2":   extract_abide2,
    "adhd200":  extract_adhd200,
    "corr":     extract_corr,
    "fcon1000": extract_fcon1000,
    "hbn":      extract_hbn,
    "nki":      extract_nki,
    "nki2":     extract_nki2,
    "hcp-ya":   extract_hcp_ya,
    "ppmi":     extract_ppmi,
}
# ADNI and 'concat' need special argument handling → not in EXTRACTORS dict
ALL_DATASETS = list(EXTRACTORS.keys()) + ["adni", "concat"]


def main():
    parser = argparse.ArgumentParser(description="Build harmonized metadata tables")
    parser.add_argument("dataset", choices=ALL_DATASETS)
    parser.add_argument("--input", "-i", type=Path,
                        help="Dataset metadata directory (or metadata root for 'concat')")
    parser.add_argument("--loni", type=Path, help="ADNI LONI CSV dir")
    parser.add_argument("--xml", type=Path, help="ADNI per-scan XML dir")
    parser.add_argument("--output", "-o", type=Path, required=True, help="Output CSV path")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")

    if args.dataset == "concat":
        if not args.input:
            parser.error("--input required for concat (metadata root dir)")
        df = concat_all(args.input)
    elif args.dataset == "adni":
        if not args.loni or not args.xml:
            parser.error("--loni and --xml both required for adni")
        df = extract_adni(args.loni, args.xml)
    else:
        if not args.input:
            parser.error(f"--input required for {args.dataset}")
        df = EXTRACTORS[args.dataset](args.input)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    log.info(f"Wrote {len(df)} rows to {args.output}")

    # Quick summary
    log.info(f"Subjects: {df['subject_id'].nunique()}")
    log.info(f"Dx distribution:\n{df['dx'].value_counts(dropna=False).to_string()}")
    log.info(f"Sex distribution:\n{df['sex'].value_counts(dropna=False).to_string()}")
    if df["age_at_scan"].notna().any():
        log.info(f"Age: mean={df['age_at_scan'].mean():.1f}, "
                 f"min={df['age_at_scan'].min():.1f}, max={df['age_at_scan'].max():.1f}")


if __name__ == "__main__":
    main()
