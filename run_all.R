# Master script -- run from the repository root:  Rscript run_all.R
# 1. Builds the analysis-ready survey table and the photograph-metric table from
#    the public data deposit (see data/README.md) unless they already exist.
# 2. Builds the survey-date covariate table from its online sources unless it
#    already exists (see analysis_22_weather_aqi/data/README.md).
# 3. Runs every analysis module that feeds a display item or a reported number.
# 4. Regenerates all main and supplementary figures and tables.

core <- "analysis_01_clean_core/outcome/clean_data/clean_core_analysis.csv"
if (!file.exists(core)) {
  message("== build/build_core.py (assembling the survey table)")
  status <- system2("python3", "build/build_core.py")
  if (status != 0 || !file.exists(core)) stop("build_core.py failed -- see data/README.md")
}
pm <- "analysis_23_photo_quality_metrics/outcome/photo_metrics_per_id.csv"
if (!file.exists(pm)) {
  message("== build/build_photo_metrics.py (assembling the photograph-metric table)")
  status <- system2("python3", "build/build_photo_metrics.py")
  if (status != 0 || !file.exists(pm)) stop("build_photo_metrics.py failed -- see data/README.md")
}
wx <- "analysis_22_weather_aqi/data/date_covariates.csv"
if (!file.exists(wx)) {
  message("== analysis_22_weather_aqi/code/fetch_date_covariates.py (date covariates)")
  status <- try(system2("python3", "analysis_22_weather_aqi/code/fetch_date_covariates.py"))
  if (inherits(status, "try-error") || status != 0 || !file.exists(wx))
    message("Date covariates unavailable: analyses 22 and 27 and Supplementary ",
            "Tables S16-S17 are skipped (see analysis_22_weather_aqi/data/README.md).")
}

modules <- c(
  "analysis_02_photo_qa/code/photo_qa.R",
  "analysis_03_dimension_reduction/code/dimension_reduction.R",
  "analysis_04_variance_anatomy/code/variance_anatomy.R",
  "analysis_05_visual_to_evaluation/code/visual_to_evaluation.R",
  "analysis_06_congruence/code/congruence.R",
  "analysis_08_variance_partitioning/code/variance_partitioning.R",
  "analysis_10_class_vs_composition/code/class_vs_composition.R",
  "analysis_11_element_candidates/code/element_candidates.R",
  "analysis_12_robustness/code/robustness.R",
  "analysis_13_context_composites/code/context_composites.R",
  "analysis_14_opportunity_conditioned/code/opportunity_conditioned.R",
  "analysis_15_congruence_conflict/code/congruence_conflict.R",
  "analysis_16_eventfulness_pathway/code/eventfulness_pathway.R",
  "analysis_17_season_phenology/code/season_phenology.R",
  "analysis_18_greenery_composition/code/greenery_composition.R",
  "analysis_19_loudness_conditioned_vision/code/loudness_conditioned_vision.R",
  "analysis_20_auditory_phenology/code/auditory_phenology.R",
  "analysis_21_familiarity_subgroups/code/familiarity_subgroups.R",
  "analysis_22_weather_aqi/code/weather_aqi.R",
  "analysis_24_metric_validation/code/validation.R",
  "analysis_25_chromatic_phenology/code/chromatic_phenology.R",
  "analysis_26_quality_outcome_scan/code/quality_scan.R",
  "analysis_27_sky_quality_interaction/code/sky_quality.R",
  "analysis_28_cv_readout/code/cv_readout.R",
  "analysis_29_blueness_integration/code/blueness_integration.R",
  "analysis_30_curvature_check/code/curvature_check.R",
  "analysis_31_panorama_clustering/code/panorama_clustering.R",
  "analysis_32_within_between/code/within_between.R",
  "analysis_33_water_sound_adjusted/code/water_sound_adjusted.R",
  "analysis_34_site_heldout_mapping/code/site_heldout_mapping.R"
)
needs_wx <- c("analysis_22_weather_aqi/code/weather_aqi.R",
              "analysis_27_sky_quality_interaction/code/sky_quality.R")
if (!file.exists(wx)) modules <- setdiff(modules, needs_wx)

display <- c(
  "display/build_display_data.R",
  "display/fig1_boundary.R",
  "display/fig2_quality.R",
  "display/fig3_eventfulness.R",
  "display/fig4_condition.R",
  "display/si_tables.R"
)

for (d in file.path(dirname(dirname(modules)), "outcome"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

for (s in c(modules, display)) {
  message("== ", s)
  status <- system2("Rscript", s)
  if (status != 0) stop("FAILED: ", s)
}
message("All display items regenerated (figures/, tables/ and output/).")
