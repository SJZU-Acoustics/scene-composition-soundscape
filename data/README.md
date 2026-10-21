# Data folder

One file belongs here; it is not shipped with the repository.

| File | Status | Contents |
|---|---|---|
| `SSID_China_urban_soundscape_survey_2018_2025.xlsx` | **download** from Mendeley Data, DOI [10.17632/wr5hrbpjbs](https://doi.org/10.17632/wr5hrbpjbs) (CC BY 4.0), **version 4 or later** | The survey workbook. Sheet `Survey_2018_2025` holds the 2,821 respondent records from the 13 urban public spaces: questionnaire, demographics, the per-channel acoustic/psychoacoustic indicator battery and the semantic-segmentation area shares of the 150 scene categories. Sheet `Photo_metrics` (from version 4) holds, for the 2,269 records with an analysed photograph, the panorama group, the number of photographs, the capture-system class and the thirteen photograph metrics the paper reports |

Two small scaffold files ship in `_shared/`: `legacy_0519_ids.csv` (the 2,373
respondent IDs of the earlier curated sample, with the age values digitised
for that sample) and `site_typology.csv` (the site functional typology).

`build/build_core.py` and `build/build_photo_metrics.py` (both called
automatically by `run_all.R`) assemble these into
`analysis_01_clean_core/outcome/clean_data/clean_core_analysis.csv`, the survey
table every analysis reads, and
`analysis_23_photo_quality_metrics/outcome/photo_metrics_per_id.csv`, the
photograph-metric table the condition analyses join to it.
