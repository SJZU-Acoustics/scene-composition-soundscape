"""Build the analysis-ready survey table for the paper "The limits of visual
scene composition in explaining and predicting urban soundscape evaluation"
from the public data deposit.

Inputs:
  - data/SSID_China_urban_soundscape_survey_2018_2025.xlsx -- the survey
    workbook, downloaded from Mendeley Data
    (https://doi.org/10.17632/wr5hrbpjbs); sheet Survey_2018_2025 holds the
    2,821 respondent records from the 13 China sites
  - _shared/legacy_0519_ids.csv (shipped) -- the 2,373 respondent IDs of the
    earlier curated sample used by the group's previous papers, with the age
    values digitised for that sample
  - _shared/site_typology.csv (shipped) -- site functional typology

Output:
  analysis_01_clean_core/outcome/clean_data/clean_core_analysis.csv
  (the single table every analysis reads) plus audit tables in
  analysis_01_clean_core/outcome/audit_tables/.

The transformation code is carried verbatim from the working pipeline so the
public build reproduces the paper's table exactly.

Key rules:
- Sample flags: 问卷问题 (C), 音频问题 (D), 音频是否分析 (E), 照片是否分析 (F).
  include_core = not questionnaire-flagged, EXCEPT rows in the legacy 0519
  sample (IDs 2149/2432/2683 were manually vetted and kept by the 0519
  curation, so legacy membership overrides the flag).
- legacy_0519 marks the 2,373 rows used by projects 10/12/14 (snapshot in
  _shared/legacy_0519_ids.csv, taken from P14's clean core).
- ISO coordinates recomputed from the eight PAQ items with the SATP cmn
  adjusted angles (Aletta et al. 2024, Table 4), soundscapy-style centred and
  normalised to [-1, 1]. Since 0716 the workbook's 修正PISO/EISO are themselves
  recomputed on the 116 corrected Zhongshan-2018 rows, so recompute and
  workbook should agree everywhere to rounding (asserted in the audit).
- RAQ: the six-batch answer-scale reversal is CONFIRMED and applied inside the
  workbook since 0716 (565 rows, by corrected batch assignment) - the builder
  reads RAQ as-is and must NOT reverse again. raq_batch_reversed keeps the
  batch flag for sensitivity analyses; a coupling guard raises if any large
  batch still looks flipped.
- Visitor status (round-3 settled rules, 2026-07-17): 身份=0 -> missing
  (global rule); waves 2-6 code 1 city resident / 2 tourist / 3 other; waves
  7-10 code 1 nearby resident / 2 city non-nearby / 3 tourist / 4 other.
  Harmonised three-tier visitor_group (local/tourist/other): old 1 <-> new
  1+2, old 2 <-> new 3, old 3 <-> new 4; residency_fine preserves the
  nearby / city-non-nearby split available only in the new-coding waves.
- Activity (round-3): waves 9/10 asked an open activity question; their text
  answers stay in the workbook. activity_state is the numeric 状态
  (1 staying / 2 arriving / 3 leaving) with text answers set missing;
  activity_text carries the open answers as a separate variable.
- Acoustic battery: available-channel rule (Prof Zhang, 2026-08-07, replacing
  the former binaural gate): the analysis only ever uses the two-ear average,
  so binaural recording is not required and mono rows enter with their single
  valid channel. Value = the workbook 均值 where both channels exist (since the
  0724 master it is recomputed as (L+R)/2 and blanked on monaural rows —
  asserted below), otherwise the one valid channel (single-ear sessions at
  中山公园/中山广场/大东广场/奥林匹克公园/劝工园/方塔园, incl. the R11
  left-only row 1358). audio_ok still gates; acoustic_channels (2/1/0) carries
  provenance and binaural_ok is retained for sensitivity splits. On binaural
  rows the left channel alone correlates 0.988 with the two-ear mean (LAeq;
  mean offset 0.000 dB, SD 1.0 dB), so single-channel values are comparable.
- Photo variables: the 150 ADE20K area shares carried through as-is, plus
  theory-driven scene composites and a row-sum audit column.
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = REPO_ROOT / "data" / "SSID_China_urban_soundscape_survey_2018_2025.xlsx"
SHARED_DIR = REPO_ROOT / "_shared"
LEGACY_IDS_CSV = SHARED_DIR / "legacy_0519_ids.csv"
SITE_TYPOLOGY_CSV = SHARED_DIR / "site_typology.csv"
OUTCOME_DIR = REPO_ROOT / "analysis_01_clean_core" / "outcome"
CLEAN_DIR = OUTCOME_DIR / "clean_data"
AUDIT_DIR = OUTCOME_DIR / "audit_tables"

SHEET = "Survey_2018_2025"

WAVE_LABELS = {
    1: "2018 autumn", 2: "2019 spring", 3: "2019 summer", 4: "2019 autumn",
    5: "2020 autumn", 6: "2021 autumn", 7: "2023 spring", 8: "2023 summer",
    9: "2024 spring", 10: "2025 spring",
}

# Column labels in the deposit sheet (row 1 header) -> clean names. The
# deposit omits three internal bookkeeping columns of the working master
# (各自编号, 问卷编号, 音频文件格式); no analysis reads them.
RENAME = {
    "ID": "id",
    "问卷问题\n1-声景维度相加=10": "questionnaire_problem",
    "音频问题": "audio_problem",
    "音频是否分析": "audio_ok",
    "照片是否分析": "photo_ok",
    "时间编码": "wave_code",
    "问卷日期": "survey_date",
    "环境状况": "scene_condition",
    "问卷版本": "survey_version",
    "问卷地点": "site_code_raw",
    "地点": "site_label",
    "性别": "gender_raw",
    "实际年龄": "age_raw",
    "标准年龄划分": "age_group_std",
    "细分年龄段": "age_group_fine",
    "粗分年龄段": "age_group_coarse",
    "受教育程度（修改3-大专在读     4-大学在读）": "education_raw",
    "身份": "visitor_status_raw",
    "职业情况": "occupation_raw",
    "状态": "activity_state_raw",
    "结伴状况": "companionship_raw",
    "交通噪声 ": "traffic_noise",
    "人的声音": "human_sounds",
    "自然声 ": "natural_sounds",
    "其他噪声 ": "other_noise",
    "愉悦感": "isop_equal_raw",
    "事件感": "isoe_equal_raw",
    "标准化ISO-P": "isop_equal_norm",
    "标准化ISO-E": "isoe_equal_norm",
    "修正\nPISO": "isop_src",
    "修正\nEISO": "isoe_src",
    "PAQ01.愉快的\np": "pleasant",
    "PAQ02.混乱的\nch": "chaotic",
    "PAQ03.有活力的\nv": "vibrant",
    "PAQ04.平淡无事的\nu": "uneventful",
    "PAQ05.平静的\nca": "calm",
    "PAQ06.恼人的\na": "annoying",
    "PAQ07.多事件的\ne": "eventful",
    "PAQ08.单调的\nm": "monotonous",
    "声环境综合评价": "overall_soundscape",
    "声环境协调度": "harmony",
    "吵闹程度": "perceived_loudness",
    "到访频率": "visit_frequency_current",
    "希望频率": "visit_frequency_desired",
    "WHO情绪自评价总分": "who5_total",
    "幸福感分组": "who5_group_code",
    "RAQ01.放松感": "raq_relaxation",
    "RAQ02.吸引力": "raq_attraction",
    "RAQ03.趣味性": "raq_interest",
    "RAQ04.舒适性": "raq_comfort",
}

WHO_ITEM_COLS = {
    "WHO01.我感到很高兴愉快，精神饱满。": "who01",
    "WHO0.我感到平静放松。": "who02",
    "WHO03.我感到精力充沛。": "who03",
    "WHO04.我早上醒来感到神清气爽，焕然一新。\n（早上醒来神清气爽）": "who04",
    "WHO05.日常生活中充满了令我感兴趣的事情。\n（充满兴趣）": "who05",
}

# Acoustic 均值 columns -> clean names (average channel, binaural-gated).
ACOUSTIC_MEAN = {
    "LAeq-均值": "laeq_a",
    "L5-均值": "l5_a",
    "L10-均值": "l10_a",
    "L50-均值": "l50_a",
    "L90-均值": "l90_a",
    "L95-均值": "l95_a",
    "L5-L95-均值": "l5_l95_a",
    "N-均值": "loudness_n",
    "N5-均值": "loudness_n5",
    "N10-均值": "loudness_n10",
    "N50-均值": "loudness_n50",
    "N90-均值": "loudness_n90",
    "N95-均值": "loudness_n95",
    "N5/N95-均值": "loudness_n5_n95",
    "S-均值": "sharpness_s",
    "S5-均值": "sharpness_s5",
    "S95-均值": "sharpness_s95",
    "R-均值": "roughness_r",
    "R10-均值": "roughness_r10",
    "R50-均值": "roughness_r50",
    "F-均值": "fluctuation_f",
    "F10-均值": "fluctuation_f10",
    "F50-均值": "fluctuation_f50",
    "TC-均值": "tonality_tc",
    "TA-均值": "tonality_ta",
}

# RAQ batch reversal - CONFIRMED by the data-curation student (0716 reply) and
# applied INSIDE the workbook since 0716 (565 rows, by corrected batch
# assignment). The builder therefore reads RAQ as-is: applying 6-x here again
# would double-reverse. The batch set is kept only as the raq_batch_reversed
# provenance flag (A11 drop-batch sensitivity); wave codes below are the
# corrected (0716+) assignments.
RAQ_REVERSED_BATCHES = {
    ("中山公园", 7), ("大东广场", 2), ("大东广场", 4), ("大东广场", 5),
    ("劳模园", 9), ("方塔园", 8),
}

# The 150 ADE20K photo columns run from 'wall' to 'flag' in sheet order.
PHOTO_FIRST, PHOTO_LAST = "wall", "flag"

# Theory-driven scene composites (area-share sums). Categories not listed in
# any group fall into comp_other (largely indoor ADE20K classes).
PHOTO_COMPOSITES = {
    "comp_greenery": ["tree", "grass", "plant", "palm", "flower", "field"],
    "comp_sky": ["sky"],
    "comp_water": ["water", "sea", "river", "lake", "waterfall", "fountain",
                    "swimmingpool"],
    "comp_paving": ["floor", "road", "sidewalk", "path", "earth", "sand",
                     "dirttrack", "runway", "land", "stairs", "stairway",
                     "step"],
    "comp_built": ["wall", "building", "house", "skyscraper", "hovel", "tower",
                    "bridge", "fence", "railing", "bannister", "column",
                    "ceiling", "door", "windowpane", "awning", "canopy",
                    "booth", "tent", "grandstand", "stage", "pier",
                    "escalator"],
    "comp_people": ["person"],
    "comp_vehicles": ["car", "bus", "truck", "van", "bicycle", "minibike",
                       "boat", "ship", "airplane"],
    "comp_furniture": ["bench", "chair", "seat", "armchair", "swivelchair",
                        "stool", "sofa", "table", "desk", "coffeetable",
                        "pole", "streetlight", "signboard", "trafficlight",
                        "sculpture", "ashcan", "basket", "bulletinboard",
                        "poster", "tradename", "flag", "lamp", "light",
                        "chandelier", "plaything", "ball", "barrel", "bag",
                        "bottle", "box", "food", "pot", "vase", "tray",
                        "plate", "glass", "towel", "blanket", "cushion",
                        "pillow", "apparel", "book", "clock", "fan",
                        "screen", "crtscreen", "monitor",
                        "televisionreceiver", "computer", "arcademachine"],
    "comp_terrain": ["rock", "mountain", "hill"],
    "comp_animal": ["animal"],
}

# --- SATP cmn adjusted-angle ISO coordinates (soundscapy add_iso_coords) ----
PAQ_PROJECTION_ORDER = ["pleasant", "vibrant", "eventful", "chaotic",
                        "annoying", "monotonous", "uneventful", "calm"]
SATP_CMN_ANGLES = (0, 18, 38, 154, 171, 196, 217, 318)
SATP_VAL_RANGE = (1, 5)


def satp_iso_coords(values: np.ndarray) -> tuple[float, float]:
    v = np.asarray(values, dtype=float)
    if np.any(np.isnan(v)):
        return np.nan, np.nan
    v = v - (v.max() + v.min()) / 2.0
    rad = np.deg2rad(np.asarray(SATP_CMN_ANGLES, dtype=float))
    scale = max(SATP_VAL_RANGE) - min(SATP_VAL_RANGE)
    iso_pl = np.sum(np.cos(rad) * v) / (scale / 2.0 * np.sum(np.abs(np.cos(rad))))
    iso_ev = np.sum(np.sin(rad) * v) / (scale / 2.0 * np.sum(np.abs(np.sin(rad))))
    return float(iso_pl), float(iso_ev)


def parse_age(raw) -> float:
    """Parse the handful of non-numeric 实际年龄 entries ('60+', '20-30', ...)."""
    if raw is None or raw == "":
        return np.nan
    if isinstance(raw, (int, float)):
        return float(raw)
    s = str(raw).strip()
    m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", s)
    if m:
        return (float(m.group(1)) + float(m.group(2))) / 2.0
    m = re.match(r"(\d+)", s)
    if m and any(tok in s for tok in ("+", "多", "左右", "岁")):
        return float(m.group(1))
    if "将近" in s:
        m2 = re.search(r"(\d+)", s)
        return float(m2.group(1)) - 2 if m2 else np.nan
    return np.nan  # 不愿透露 / 退休了 / 70年生 stay missing


def to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def read_sheet() -> pd.DataFrame:
    from python_calamine import CalamineWorkbook

    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"{DATA_PATH.name} not found in data/ -- download it from "
            "https://doi.org/10.17632/wr5hrbpjbs (see data/README.md)")
    wb = CalamineWorkbook.from_path(DATA_PATH)
    rows = wb.get_sheet_by_name(SHEET).to_python()
    header = [str(h) for h in rows[0]]
    id_col = header.index("ID")
    data = []
    for r in rows[1:]:
        r = list(r) + [""] * (len(header) - len(r))
        if r[id_col] in ("", None):
            continue
        data.append(r[:len(header)])
    df = pd.DataFrame(data, columns=header)
    if len(df) != 2821:
        raise RuntimeError(f"Expected 2,821 data records, got {len(df)}")
    return df


def build_core(raw: pd.DataFrame) -> pd.DataFrame:
    photo_cols = list(raw.columns[raw.columns.get_loc(PHOTO_FIRST):
                                  raw.columns.get_loc(PHOTO_LAST) + 1])
    if len(photo_cols) != 150:
        raise RuntimeError(f"Expected 150 photo columns, got {len(photo_cols)}")

    d = pd.DataFrame(index=raw.index)
    for src, dst in RENAME.items():
        if src not in raw.columns:
            raise RuntimeError(f"Missing expected column: {src!r}")
        d[dst] = raw[src]
    for src, dst in WHO_ITEM_COLS.items():
        d[dst] = to_num(raw[src])

    numeric_cols = [
        "id", "questionnaire_problem", "audio_ok", "photo_ok", "wave_code",
        "scene_condition", "survey_version", "site_code_raw", "gender_raw",
        "age_group_std", "age_group_fine", "age_group_coarse", "education_raw",
        "visitor_status_raw", "occupation_raw", "companionship_raw",
        "traffic_noise", "human_sounds", "natural_sounds", "other_noise",
        "isop_equal_raw", "isoe_equal_raw", "isop_equal_norm", "isoe_equal_norm",
        "isop_src", "isoe_src",
        "pleasant", "chaotic", "vibrant", "uneventful", "calm", "annoying",
        "eventful", "monotonous",
        "overall_soundscape", "harmony", "perceived_loudness",
        "visit_frequency_current", "visit_frequency_desired",
        "who5_total", "who5_group_code",
        "raq_relaxation", "raq_attraction", "raq_interest", "raq_comfort",
    ]
    for c in numeric_cols:
        d[c] = to_num(d[c])

    d["id"] = d["id"].astype(int)
    d["site_label"] = d["site_label"].astype(str).str.strip()
    d["wave_label"] = d["wave_code"].map(WAVE_LABELS)
    d["gender"] = d["gender_raw"].map({1: "Male", 2: "Female"})
    d["education_ord"] = d["education_raw"].where(d["education_raw"] != 8)
    d["age"] = d["age_raw"].map(parse_age)
    d["who5_group"] = d["who5_group_code"].map({1: "LW", 2: "HW"})

    # RAQ: workbook values are already batch-harmonised (0716+); keep the
    # batch flag for sensitivity work and mirror the values into the *_src
    # columns so older column references stay valid (src == primary now).
    rev = pd.Series(
        [(s, int(w)) in RAQ_REVERSED_BATCHES if pd.notna(w) else False
         for s, w in zip(d["site_label"], d["wave_code"])], index=d.index)
    d["raq_batch_reversed"] = rev.astype(int)
    for c in ["raq_relaxation", "raq_attraction", "raq_interest",
              "raq_comfort"]:
        d[f"{c}_src"] = d[c]
    d["raq_mean"] = d[["raq_relaxation", "raq_attraction", "raq_interest",
                       "raq_comfort"]].mean(axis=1)
    d["raq_mean_src"] = d["raq_mean"]

    # Guard against a missed or double reversal: every sizeable batch must now
    # couple positively with ISO pleasantness.
    for (s, w), g in d.groupby(["site_label", "wave_code"]):
        if len(g) < 30:
            continue
        r = g["raq_mean"].corr(g["isop_src"])
        if pd.notna(r) and r < -0.2:
            raise RuntimeError(
                f"Batch ({s}, {w}) still couples negatively with ISO-P "
                f"(r={r:.2f}) - RAQ reversal state is inconsistent")

    # Visitor status (round-3 settled rules, 2026-07-17). 0 = question absent
    # or unrecorded -> missing (global rule). Two coding eras:
    #   waves 2-6 (2019-2021): 1 city resident / 2 tourist / 3 other
    #   waves 7-10 (2023+):    1 nearby res. / 2 city non-nearby / 3 tourist / 4 other
    # The two hospital-visitor batches were resolved at session level in the
    # 0717 workbook (奥林匹克 2020-09-20 -> 1; 中山广场 2020-10-03 -> missing).
    vs = d["visitor_status_raw"].where(d["visitor_status_raw"] != 0)
    new_era = d["wave_code"] >= 7
    d["visitor_group"] = np.select(
        [(~new_era & vs.eq(1)) | (new_era & vs.isin([1, 2])),
         (~new_era & vs.eq(2)) | (new_era & vs.eq(3)),
         (~new_era & vs.eq(3)) | (new_era & vs.eq(4))],
        ["local", "tourist", "other"], default=None)
    d["residency_fine"] = np.select(
        [new_era & vs.eq(1), new_era & vs.eq(2)],
        ["nearby", "city_non_nearby"], default=None)
    bad_vs = vs.notna() & d["visitor_group"].isna()
    if bad_vs.any():
        raise RuntimeError(
            f"Unmapped 身份 values: {sorted(vs[bad_vs].unique())} "
            f"(ids {sorted(d.loc[bad_vs, 'id'].head(10))})")

    # Activity (round-3): numeric 状态 = 1 staying / 2 arriving / 3 leaving;
    # the waves-9/10 open text answers become a separate activity variable and
    # leave activity_state missing.
    act_raw = d["activity_state_raw"].astype(str).str.strip()
    act_num = to_num(d["activity_state_raw"])
    d["activity_state"] = act_num.where(act_num.isin([1, 2, 3]))
    d["activity_text"] = act_raw.where(
        act_num.isna() & ~act_raw.isin(["", "nan", "None"]))

    # Legacy 0519 membership + age fallback from the 0519 cleaned values
    legacy = pd.read_csv(LEGACY_IDS_CSV)
    d["legacy_0519"] = d["id"].isin(set(legacy["ID"])).astype(int)
    age_0519 = legacy.set_index("ID")["age_0519"]
    need_age = d["age"].isna() & d["legacy_0519"].eq(1)
    d.loc[need_age, "age"] = d.loc[need_age, "id"].map(age_0519)

    # ID 1204: the recorded 实际年龄 is an impossible value verified as a
    # faithful transcription of the paper questionnaire; the age is set
    # missing and the row retained. Applied after the 0519 fallback so no
    # path can refill it.
    d.loc[d["id"].eq(1204), "age"] = np.nan

    # Inclusion: questionnaire flag excludes, unless vetted legacy row
    d["questionnaire_flagged"] = d["questionnaire_problem"].eq(1).astype(int)
    d["include_core"] = ((d["questionnaire_flagged"] == 0)
                          | (d["legacy_0519"] == 1)).astype(int)

    # ISO recompute (SATP cmn adjusted angles, soundscapy normalisation)
    paq = d[PAQ_PROJECTION_ORDER].to_numpy()
    coords = np.array([satp_iso_coords(row) for row in paq])
    d["isop"] = coords[:, 0]
    d["isoe"] = coords[:, 1]
    complete = d[PAQ_PROJECTION_ORDER].notna().all(axis=1)
    if d.loc[complete, "isop"].isna().any():
        raise RuntimeError("ISO recompute failed on complete-PAQ rows")
    dev = (d.loc[complete, "isop"] - d.loc[complete, "isop_src"]).abs()
    iso_src_max_dev = float(dev.max())

    # Acoustics: available-channel battery (Prof Zhang, 2026-08-07 — binaural
    # not required; mono rows enter with their single valid channel)
    d["laeq_l"] = to_num(raw["LAeq-L"])
    d["laeq_r"] = to_num(raw["LAeq-R"])
    d["binaural_ok"] = (d["laeq_l"].notna() & d["laeq_r"].notna()).astype(int)
    d["acoustic_channels"] = (d["laeq_l"].notna().astype(int)
                              + d["laeq_r"].notna().astype(int))
    for src, dst in ACOUSTIC_MEAN.items():
        mean_vals = to_num(raw[src])
        lch = to_num(raw[src.replace("-均值", "-L")])
        rch = to_num(raw[src.replace("-均值", "-R")])
        if (mean_vals.notna() & (lch.isna() | rch.isna())).any():
            raise RuntimeError(f"{src} present without both channels")
        if (mean_vals.isna() & lch.notna() & rch.notna()).any():
            raise RuntimeError(f"{src} missing despite both channels")
        d[dst] = mean_vals.fillna(lch).fillna(rch).where(d["audio_ok"].eq(1))
    d["l10_l90_a"] = d["l10_a"] - d["l90_a"]
    d["loudness_n10_n90"] = d["loudness_n10"] - d["loudness_n90"]

    # Source-composition indices (P14 conventions)
    d["source_total"] = (d["traffic_noise"] + d["human_sounds"]
                          + d["natural_sounds"] + d["other_noise"])
    d["natural_share"] = d["natural_sounds"] / d["source_total"]
    d["noise_load"] = (d["traffic_noise"] + d["other_noise"]) / 2
    d["mech_vs_natural"] = ((d["traffic_noise"] + d["other_noise"]) / 2
                             - d["natural_sounds"])

    # Photo variables: shares as-is + composites + row-sum audit
    photo = raw[photo_cols].apply(to_num)
    photo.columns = [f"ph_{c}" for c in photo_cols]
    d = pd.concat([d, photo], axis=1)
    d["photo_row_sum"] = photo.sum(axis=1, min_count=1)
    d["photo_present"] = photo.notna().any(axis=1).astype(int)
    grouped = set()
    for comp, members in PHOTO_COMPOSITES.items():
        cols = [f"ph_{m}" for m in members]
        missing = [c for c in cols if c not in d.columns]
        if missing:
            raise RuntimeError(f"Composite {comp} missing columns: {missing}")
        d[comp] = d[cols].sum(axis=1, min_count=1)
        grouped.update(cols)
    rest = [c for c in photo.columns if c not in grouped]
    d["comp_other"] = d[rest].sum(axis=1, min_count=1)
    d["comp_nature"] = d["comp_greenery"] + d["comp_water"]

    # Site typology join (on 地点 name)
    typ = pd.read_csv(SITE_TYPOLOGY_CSV)
    d = d.merge(typ[["site_label", "site_code", "space_type", "space_class"]],
                on="site_label", how="left", validate="many_to_one")
    if d["space_class"].isna().any():
        bad = sorted(d.loc[d["space_class"].isna(), "site_label"].unique())
        raise RuntimeError(f"Sites missing from typology: {bad}")
    d["site_key"] = ("CN_" + d["site_code"].astype(int).astype(str).str.zfill(2)
                     + "_" + d["site_label"])
    d = d.copy()
    d.attrs["iso_src_max_dev"] = iso_src_max_dev
    return d


def export_audits(d: pd.DataFrame) -> None:
    inc = d[d["include_core"] == 1]
    # ISO provenance: since 0716 the workbook 修正PISO/EISO are recomputed on
    # the 116 corrected Zhongshan-2018 rows too, so workbook and recompute
    # must agree everywhere to 2-dp rounding.
    dev = np.maximum((d["isop"] - d["isop_src"]).abs(),
                     (d["isoe"] - d["isoe_src"]).abs())
    stale = dev > 0.02
    if stale.sum() > 0:
        raise RuntimeError(
            f"{int(stale.sum())} rows deviate >0.02 from workbook 修正PISO/EISO "
            f"- unexpected under the settled 0717 workbook")
    pd.DataFrame([
        {"field": "data_workbook", "value": DATA_PATH.name},
        {"field": "rows_raw", "value": len(d)},
        {"field": "iso_src_stale_rows_gt_0.02", "value": int(stale.sum())},
        {"field": "iso_src_max_dev",
         "value": float(dev.max())},
        {"field": "rows_include_core", "value": len(inc)},
        {"field": "rows_legacy_0519", "value": int(d["legacy_0519"].sum())},
        {"field": "rows_extra_included",
         "value": int(((d["legacy_0519"] == 0) & (d["include_core"] == 1)).sum())},
        {"field": "rows_extra_excluded_flag",
         "value": int(((d["legacy_0519"] == 0) & (d["include_core"] == 0)).sum())},
        {"field": "rows_acoustic",
         "value": int(inc["laeq_a"].notna().sum())},
        {"field": "rows_acoustic_binaural",
         "value": int((inc["laeq_a"].notna()
                       & inc["binaural_ok"].eq(1)).sum())},
        {"field": "rows_acoustic_single_channel",
         "value": int((inc["laeq_a"].notna()
                       & inc["binaural_ok"].eq(0)).sum())},
        {"field": "rows_photo_valid",
         "value": int((inc["photo_ok"].eq(1) & inc["photo_present"].eq(1)).sum())},
        {"field": "iso_src_max_abs_dev_recompute_vs_workbook",
         "value": d.attrs["iso_src_max_dev"]},
        {"field": "sites", "value": int(d["site_key"].nunique())},
    ] + [
        {"field": "note",
         "value": "isop/isoe (recomputed from corrected PAQ, SATP cmn angles,"
                  " soundscapy normalisation) are authoritative; isop_src/"
                  "isoe_src kept for provenance only"},
    ]).to_csv(AUDIT_DIR / "run_metadata.csv", index=False, encoding="utf-8-sig")

    by_site = (inc.groupby(["site_key", "site_label", "space_class",
                            "space_type"], dropna=False)
               .agg(rows=("id", "size"),
                    legacy_n=("legacy_0519", "sum"),
                    extra_n=("legacy_0519", lambda s: int((s == 0).sum())),
                    raq_mean=("raq_mean", "mean"),
                    isop_mean=("isop", "mean"),
                    isoe_mean=("isoe", "mean"),
                    laeq_n=("laeq_a", lambda s: int(s.notna().sum())),
                    laeq_mean=("laeq_a", "mean"),
                    photo_n=("photo_present", "sum"))
               .reset_index().sort_values("site_key"))
    by_site.to_csv(AUDIT_DIR / "sample_counts_by_site.csv", index=False,
                   encoding="utf-8-sig")

    by_wave = (inc.groupby(["wave_code", "wave_label"], dropna=False)
               .agg(rows=("id", "size"), legacy_n=("legacy_0519", "sum"),
                    sites=("site_key", "nunique"))
               .reset_index())
    by_wave.to_csv(AUDIT_DIR / "sample_counts_by_wave.csv", index=False,
                   encoding="utf-8-sig")

    key = ["raq_relaxation", "raq_attraction", "raq_interest", "raq_comfort",
           "traffic_noise", "human_sounds", "natural_sounds", "other_noise",
           "isop", "isoe", "overall_soundscape", "harmony",
           "perceived_loudness", "visit_frequency_current",
           "visit_frequency_desired", "who5_total", "age", "gender",
           "education_ord", "laeq_a", "loudness_n5", "sharpness_s",
           "roughness_r", "fluctuation_f", "tonality_tc", "comp_greenery",
           "comp_sky", "comp_people", "photo_row_sum"]
    rows = []
    for name, frame in [("include_core", inc), ("all_raw", d)]:
        for c in key:
            miss = int(frame[c].isna().sum())
            rows.append({"sample": name, "variable": c, "n": len(frame),
                         "n_missing": miss,
                         "pct_missing": round(miss / len(frame) * 100, 2)})
    pd.DataFrame(rows).to_csv(AUDIT_DIR / "missingness_key_variables.csv",
                              index=False, encoding="utf-8-sig")

    extras = d[d["legacy_0519"] == 0]
    (extras.groupby(["site_label", "wave_label", "include_core"])
     .agg(rows=("id", "size"))
     .reset_index()
     .to_csv(AUDIT_DIR / "extra_records_by_site_wave.csv", index=False,
             encoding="utf-8-sig"))

    flagged = d[d["questionnaire_flagged"] == 1]
    flagged[["id", "site_label", "wave_label", "legacy_0519", "include_core"]]\
        .to_csv(AUDIT_DIR / "questionnaire_flagged_rows.csv", index=False,
                encoding="utf-8-sig")

    over = inc[inc["photo_row_sum"] > 1.05]
    over[["id", "site_label", "photo_row_sum"]].sort_values("id")\
        .to_csv(AUDIT_DIR / "photo_rowsum_over_105.csv", index=False,
                encoding="utf-8-sig")

    # Visitor / activity derivation audit (round-3 rules)
    va = (d.assign(era=np.where(d["wave_code"] >= 7, "new_2023plus",
                                "old_2019_2021"))
          .groupby(["era", "visitor_group"], dropna=False)
          .agg(rows=("id", "size")).reset_index())
    va.to_csv(AUDIT_DIR / "visitor_group_by_era.csv", index=False,
              encoding="utf-8-sig")
    act = (d.groupby(["wave_code", "wave_label"], dropna=False)
           .agg(state_numeric=("activity_state", lambda s: int(s.notna().sum())),
                activity_text=("activity_text", lambda s: int(s.notna().sum())),
                fine_residency=("residency_fine", lambda s: int(s.notna().sum())))
           .reset_index())
    act.to_csv(AUDIT_DIR / "activity_residency_by_wave.csv", index=False,
               encoding="utf-8-sig")

    # RAQ batch-reversal evidence: per site-wave coupling with the
    # pleasantness family (src == primary since the workbook applies the
    # confirmed reversal; the table verifies every batch couples positively)
    ev = []
    for (s, w), g in d.groupby(["site_label", "wave_code"]):
        if len(g) < 10:
            continue
        ev.append({
            "site_label": s, "wave_code": w, "n": len(g),
            "flagged_reversed": int((s, int(w)) in RAQ_REVERSED_BATCHES),
            "raq_mean_src": round(g["raq_mean_src"].mean(), 3),
            "raq_mean_harmonised": round(g["raq_mean"].mean(), 3),
            "r_src_isop": round(g["raq_mean_src"].corr(g["isop"]), 3),
            "r_harmonised_isop": round(g["raq_mean"].corr(g["isop"]), 3),
            "r_src_oa": round(g["raq_mean_src"].corr(g["overall_soundscape"]), 3),
        })
    pd.DataFrame(ev).sort_values(["site_label", "wave_code"]).to_csv(
        AUDIT_DIR / "raq_batch_reversal_evidence.csv", index=False,
        encoding="utf-8-sig")


def main() -> None:
    CLEAN_DIR.mkdir(parents=True, exist_ok=True)
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore", FutureWarning)
    raw = read_sheet()
    d = build_core(raw)
    d.to_csv(CLEAN_DIR / "clean_core_analysis.csv", index=False,
             encoding="utf-8-sig")
    export_audits(d)
    print(f"rows: {len(d)}; include_core: {int(d['include_core'].sum())}; "
          f"legacy: {int(d['legacy_0519'].sum())}; "
          f"iso max dev vs workbook: {d.attrs['iso_src_max_dev']:.2e}")


if __name__ == "__main__":
    main()
