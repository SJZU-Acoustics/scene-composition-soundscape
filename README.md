Code for reproducing the statistical analyses, figures and tables for the manuscript "The limits of visual scene composition in explaining and predicting urban soundscape evaluation".

## Requirements

- R 4.5+ (developed and verified on R 4.5)

- CRAN packages: `tidyverse`, `lmerTest`, `broom`, `broom.mixed`, `performance`, `ranger`, `patchwork`, `scales`, `ggtext`, `ragg`

  ```r
  install.packages(c("tidyverse", "lmerTest", "broom", "broom.mixed", "performance",
                     "ranger", "patchwork", "scales", "ggtext", "ragg"))
  ```

- Python 3.10+ with `pandas`, `numpy` and `python-calamine` (used once, to assemble the analysis tables), plus `requests` and `meteostat` for the optional date-covariate table:

  ```bash
  pip install pandas numpy python-calamine requests meteostat
  ```

- No non-standard hardware is required.

## Data

The analysis tables are assembled from the public survey deposit and two small scaffold files:

1. **Download** the survey workbook `SSID_China_urban_soundscape_survey_2018_2025.xlsx` from Mendeley Data, **DOI [10.17632/wr5hrbpjbs](https://doi.org/10.17632/wr5hrbpjbs)** (CC BY 4.0), **version 4 or later**, and place it in `data/`. Its sheet `Survey_2018_2025` holds the 2,821 respondent records from the 13 urban public spaces (questionnaire, demographics, the per-channel acoustic/psychoacoustic indicator battery and the semantic-segmentation area shares of the 150 scene categories); its sheet `Photo_metrics` holds the panorama group, the capture-system class and the thirteen photograph metrics of the 2,269 records with an analysed photograph.
2. `_shared/legacy_0519_ids.csv` (shipped) — the 2,373 respondent IDs of the earlier curated sample used by the group's previous papers, with the age values digitised for that sample.
3. `_shared/site_typology.csv` (shipped) — the site functional typology (park / civic square / residential; green / built).

`run_all.R` builds both analysis tables automatically on first run (or run `python3 build/build_core.py` and `python3 build/build_photo_metrics.py` yourself).

Two inputs are not in the repository:

- **Date covariates** (daily weather from Meteostat and daily air quality from the quotsoft.net/air archive, one row per survey date): `analysis_22_weather_aqi/code/fetch_date_covariates.py` builds them from their online sources, and `run_all.R` calls it when the file is absent. Without it, analyses 22 and 27 and Supplementary Tables S16 and S17 are skipped (see `analysis_22_weather_aqi/data/README.md`); nothing else depends on it.
- **The photographs.** The respondent-level photograph metrics are in the deposit; `analysis_23_photo_quality_metrics/code/extract_photo_metrics.py` documents how each metric is computed from a photograph but is not run by the pipeline. The capture-system and processing-generation statistics quoted in the Methods were computed on the photograph collection and are not regenerated here.

## File structure

- `run_all.R` — master script: assembles the analysis tables if needed, runs every analysis module that feeds a display item or a reported number, and regenerates all figures and tables.
- `build/build_core.py` — builds the survey table `analysis_01_clean_core/outcome/clean_data/clean_core_analysis.csv` from the deposit workbook; the transformation code (sample flags, SATP adjusted-angle ISO coordinates, visitor-status harmonisation, available-channel acoustic battery, scene composites) is carried verbatim from the working pipeline.
- `build/build_photo_metrics.py` — builds the photograph-metric table `analysis_23_photo_quality_metrics/outcome/photo_metrics_per_id.csv` from the deposit's `Photo_metrics` sheet.
- `_shared/helpers.R` — single data entry point (`load_rest_data()` for the survey table, `load_vis_data()` for the 2,262 photo-valid records this paper analyses), outcome and indicator metadata, and shared model helpers; `analysis_23_photo_quality_metrics/code/photo_helpers.R` joins the photograph metrics (`load_quality_data()`).
- `analysis_02` … `analysis_34` — one folder per analysis, run in numeric order. Numbers follow the working pipeline; **01 is the table build above, 23 is the photograph-metric layer (loader and extraction reference), and 07 and 09 are intentionally absent** (exploratory analyses not used in the manuscript). `02` photograph validity flags, row sums and wave structure · `03` reduction of the 150 scene categories to the eight composites · `04` within- and between-site variance of the composites · `05` composition family (eight composites on the evaluation outcomes) · `06` cross-modal congruence family · `08` variance partition, auditory measures versus composites · `10` space class versus composition · `11` element-presence family · `12` robustness (leave-one-site-out, site-cluster bootstrap, within/between decomposition, renormalised shares) · `13` constructed scene descriptors · `14` opportunity-conditioned element refits · `15` audio-visual match and conflict indices · `16` eventfulness pathway (indirect effects) · `17` season contrasts of the visual effects · `18` internal structure of greenery · `19` perceived-loudness and sound-level moderation · `20` season contrasts of the sound-source effects · `21` visitor-familiarity and visit-stage moderation · `22` date-level weather and air quality · `24` convergent validation of the photograph metrics · `25` chromatic season contrasts · `26` condition-metric family · `27` sky blueness by sky share, with weather adjustment · `28` permutation-calibrated grouped cross-validation (linear and random forest; composition and condition blocks) · `29` blueness integration (sound-adjusted, by season, by perceived loudness) · `30` curvature check · `31` panorama-group sensitivity and the green-per-greenery means · `32` within- and between-site components of the reported effects · `33` water effects with the auditory measures in the model · `34` prediction at held-out sites.
- `display/` — the manuscript display items: `style.R` (shared plot style), `build_display_data.R` (the per-panel data behind each figure, Tables 1–2), `fig1_boundary.R`, `fig2_quality.R`, `fig3_eventfulness.R`, `fig4_condition.R` (the main figures) and `si_tables.R` (the supplementary tables) as LaTeX fragments.

## Usage

From the repository root:

```bash
Rscript run_all.R
```

The full pipeline takes about five hours on a normal desktop, almost all of it in the permutation-calibrated cross-validation of analysis 28 (1,000 linear and 200 random-forest permutations, each scored on 100 folds); every other module runs in minutes. The deposit stores the photograph metrics to 15 significant digits, so the modules that read them reproduce the paper's values to display precision rather than bit for bit. Outputs are written to:

- `figures/` — the main figures (PNG, 600 dpi): `fig1_boundary.png` is Figure 1, `fig3_eventfulness.png` Figure 2, `fig2_quality.png` Figure 3 and `fig4_condition.png` Figure 4
- `tables/` — Tables 1–2; `tables/si/` — Supplementary Tables S1–S28 (LaTeX `tabular` fragments)
- `output/` — the per-panel data behind each figure (CSV)

## License

The code is released under the MIT License (see `LICENSE`). The scaffold files shipped in `_shared/` are shared under CC BY 4.0, matching the Mendeley Data deposit's licence.
