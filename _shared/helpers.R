# Shared helpers: the single data entry point for every analysis.
# load_rest_data() reads the analysis-ready survey table built by
# build/build_core.py from the public deposit (13 sites, 2,792 core records;
# ISO coordinates recomputed with the SATP adjusted angles; visitor status
# harmonised across questionnaire eras; acoustic indicators as the average of
# the valid channels; the 150 scene-category shares and their composites).
# load_vis_data() restricts it to the photo-valid records (n = 2,262) that
# this paper analyses. China arm only, one survey language, so plain
# z-scores are the working scale.
#
# Samples:
#   "core"   include_core == 1 (default)
#   "legacy" legacy_0519 == 1 (the earlier curated row set)
#   "all"    every raw row (audits only)

suppressPackageStartupMessages({
  library(tidyverse)
  library(lmerTest)
  library(broom)
  library(broom.mixed)
})

PROJECT_DIR <- "."  # run everything from the repository root
CLEAN_CSV <- file.path(PROJECT_DIR,
  "analysis_01_clean_core/outcome/clean_data/clean_core_analysis.csv")

RAQ_ITEMS <- c("raq_relaxation", "raq_attraction", "raq_interest", "raq_comfort")

OUTCOMES <- tribble(
  ~var,                      ~label,                    ~role,
  "raq_mean",                "Restorativeness (RAQ)",   "focal",
  "raq_relaxation",          "RAQ relaxation",          "focal_item",
  "raq_attraction",          "RAQ attraction",          "focal_item",
  "raq_interest",            "RAQ interest",            "focal_item",
  "raq_comfort",             "RAQ comfort",             "focal_item",
  "overall_soundscape",      "Overall evaluation",      "comparison",
  "isop",                    "ISO pleasantness",        "comparison",
  "isoe",                    "ISO eventfulness",        "comparison",
  "harmony",                 "Perceived harmony",       "comparison",
  "perceived_loudness",      "Perceived loudness",      "comparison",
  "visit_frequency_current", "Current visit frequency", "behaviour",
  "visit_frequency_desired", "Desired visit frequency", "behaviour",
)

ADJECTIVES <- c("pleasant", "vibrant", "eventful", "chaotic",
                "annoying", "monotonous", "uneventful", "calm")

SOURCES <- c("traffic_noise", "human_sounds", "natural_sounds", "other_noise")

ACOUSTIC <- tribble(
  ~var,               ~label,                      ~family,
  "laeq_a",           "LAeq",                      "level",
  "l10_l90_a",        "L10-L90 (level dynamics)",  "dynamics",
  "loudness_n",       "Loudness N",                "loudness",
  "loudness_n5",      "Loudness N5",               "loudness",
  "loudness_n10_n90", "N10-N90 (loudness range)",  "dynamics",
  "sharpness_s",      "Sharpness S",               "spectral",
  "roughness_r",      "Roughness R",               "temporal",
  "fluctuation_f",    "Fluctuation strength F",    "temporal",
  "tonality_tc",      "Tonality TC",               "spectral",
  "tonality_ta",      "Tonality TA",               "spectral",
)

VISUAL_COMPOSITES <- c("comp_greenery", "comp_sky", "comp_water", "comp_paving",
                       "comp_built", "comp_people", "comp_vehicles",
                       "comp_furniture", "comp_terrain", "comp_nature")

load_rest_data <- function(sample = c("core", "legacy", "all")) {
  sample <- match.arg(sample)
  d <- readr::read_csv(CLEAN_CSV, show_col_types = FALSE, guess_max = 3000)
  d <- d %>%
    mutate(
      site = factor(site_key),
      space_class = factor(space_class,
                           levels = c("park", "civic_square", "residential")),
      space_type = factor(space_type, levels = c("built", "green")),
      wave = factor(wave_code),
      gender = factor(gender, levels = c("Female", "Male")),
      visitor_group = factor(visitor_group,
                             levels = c("local", "tourist", "other")),
      residency_fine = factor(residency_fine,
                              levels = c("nearby", "city_non_nearby")),
      who5_z = as.numeric(scale(who5_total)),
      age_dec = (age - 40) / 10,
      photo_valid = as.integer(photo_ok == 1 & photo_present == 1)
    )
  if (sample == "core") d <- filter(d, include_core == 1)
  if (sample == "legacy") d <- filter(d, legacy_0519 == 1)
  zvars <- c(OUTCOMES$var, ADJECTIVES, SOURCES, ACOUSTIC$var,
             "natural_share", "noise_load", "mech_vs_natural",
             VISUAL_COMPOSITES)
  for (v in zvars) {
    d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
  }
  d
}

fdr <- function(p) p.adjust(p, method = "BH")

# survey_date appears in four formats across waves: "2019-04-27",
# "20230520.0", "2024.05.05", "2025.4.28"
parse_survey_date <- function(x) {
  x <- str_remove(str_trim(as.character(x)), "\\.0$")
  x <- str_replace_all(x, "\\.", "-")
  suppressWarnings(coalesce(as.Date(x, format = "%Y-%m-%d"),
                            as.Date(x, format = "%Y%m%d")))
}

tidy_term <- function(fit, term_name) {
  broom.mixed::tidy(fit, effects = "fixed") %>%
    filter(term == term_name) %>%
    select(term, estimate, std.error, statistic, df, p.value)
}

P22_DIR <- "."  # repository root (run everything from there)

# Photo-valid analytic view: core sample, valid photo rows only, plus
# row-sum-renormalised composite variants (comp_*_rn; the A02 decision) and
# their z-scores.
COMPS8 <- c("comp_greenery", "comp_sky", "comp_water", "comp_paving",
            "comp_built", "comp_people", "comp_vehicles", "comp_furniture")

load_vis_data <- function(sample = "core") {
  d <- load_rest_data(sample) %>%
    filter(photo_valid == 1) %>%
    mutate(across(all_of(c(COMPS8, "comp_terrain", "comp_animal",
                           "comp_other", "comp_nature")),
                  ~ .x / photo_row_sum, .names = "{.col}_rn"))
  for (v in paste0(c(COMPS8, "comp_nature"), "_rn")) {
    d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
  }
  d
}
