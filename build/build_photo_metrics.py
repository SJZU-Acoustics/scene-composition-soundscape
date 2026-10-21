"""Build the respondent-level photograph-metric table for the paper "The limits
of visual scene composition in explaining and predicting urban soundscape
evaluation" from the public data deposit.

Input:
  data/SSID_China_urban_soundscape_survey_2018_2025.xlsx -- sheet Photo_metrics
  (version 4 or later of the deposit; one row per record with an analysed
  photograph: panorama group, number of photographs, capture-system class and
  the thirteen photograph metrics reported in the paper)

Output:
  analysis_23_photo_quality_metrics/outcome/photo_metrics_per_id.csv
  (the table the photograph-metric analyses read, under the column names the
  analysis code uses)

The deposit carries the per-record values; the photographs themselves are not
public. analysis_23_photo_quality_metrics/code/extract_photo_metrics.py
documents how each metric is computed from a photograph.
"""

from __future__ import annotations

import csv
from pathlib import Path

from python_calamine import CalamineWorkbook

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = REPO_ROOT / "data" / "SSID_China_urban_soundscape_survey_2018_2025.xlsx"
OUT_PATH = (REPO_ROOT / "analysis_23_photo_quality_metrics" / "outcome"
            / "photo_metrics_per_id.csv")
SHEET = "Photo_metrics"

# deposit column -> analysis column
RENAME = {
    "ID": "id",
    "panorama_group": "panorama_group",
    "n_photos": "n_shots",
    "green_pixel_fraction": "veg_frac",
    "vegetation_saturation": "veg_sat",
    "vegetation_hue": "veg_hue",
    "excess_green": "exg_mean",
    "colourfulness": "colorfulness",
    "sky_pixel_fraction_top": "sky_frac_top",
    "sky_blueness": "sky_blue_ratio",
    "brightness": "lum_mean",
    "contrast": "lum_sd",
    "edge_density": "edge_density",
    "luminance_entropy": "entropy",
    "blur_laplacian_variance": "blur_lapvar",
    "highlight_clipping_fraction": "clip_hi",
}
INT_COLS = {"id", "n_shots"}


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(
            f"{DATA_PATH.name} not found in data/ -- download it from "
            "https://doi.org/10.17632/wr5hrbpjbs (see data/README.md)")
    wb = CalamineWorkbook.from_path(DATA_PATH)
    if SHEET not in wb.sheet_names:
        raise RuntimeError(
            f"sheet {SHEET} not found: the deposit file must be version 4 or later")
    rows = wb.get_sheet_by_name(SHEET).to_python()
    header = [str(h) for h in rows[0]]
    missing = [c for c in RENAME if c not in header]
    if missing:
        raise RuntimeError(f"columns missing from {SHEET}: {missing}")
    idx = {c: header.index(c) for c in RENAME}
    sys_idx = header.index("capture_system")
    out_cols = ["id", "photo_flag", "n_shots", "panorama_group", "dual_generation"] + [
        RENAME[c] for c in RENAME if RENAME[c] not in ("id", "n_shots", "panorama_group")]
    records = []
    for r in rows[1:]:
        if r[idx["ID"]] in ("", None):
            continue
        rec = {"photo_flag": 1,
               "dual_generation": int(r[sys_idx] == "dual-generation")}
        for c, name in RENAME.items():
            v = r[idx[c]]
            if name in INT_COLS:
                v = int(v)
            elif v in ("", None):
                v = "NA"
            rec[name] = v
        records.append(rec)
    if len(records) != 2269:
        raise RuntimeError(f"Expected 2,269 photograph records, got {len(records)}")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols)
        w.writeheader()
        for rec in records:
            w.writerow({k: (repr(v) if isinstance(v, float) else v)
                        for k, v in rec.items()})
    print(f"photograph records: {len(records)}; "
          f"panorama groups: {len({r['panorama_group'] for r in records})}")


if __name__ == "__main__":
    main()
