# Shared loader for the photograph-metric layer (analyses 24-29, 31, 32).
# Joins the respondent-level photograph metrics (built from the deposit's
# Photo_metrics sheet by build/build_photo_metrics.py) onto the photo-valid
# analysis view. The thirteen metrics below are the ones the paper reports;
# extract_photo_metrics.py in this folder documents their computation.

source("_shared/helpers.R")

A23_DIR <- "analysis_23_photo_quality_metrics"

PHOTO_METRICS <- c("lum_mean", "lum_sd", "veg_frac", "veg_hue", "veg_sat",
                   "exg_mean", "sky_frac_top", "sky_blue_ratio",
                   "colorfulness", "edge_density", "blur_lapvar",
                   "clip_hi", "entropy")

# quality metrics carried into the A26 scan (chromatic / luminance / structure)
QUALITY_VARS <- tribble(
  ~var,             ~label,                          ~block,
  "veg_sat",        "Vegetation saturation",         "chromatic",
  "veg_hue_dev",    "Vegetation hue (yellow shift)", "chromatic",
  "exg_mean",       "Excess green",                  "chromatic",
  "colorfulness",   "Colourfulness (HS)",            "chromatic",
  "sky_blue_ratio", "Sky blueness",                  "sky",
  "lum_mean",       "Brightness",                    "luminance",
  "lum_sd",         "RMS contrast",                  "luminance",
  "edge_density",   "Edge density",                  "structure",
  "entropy",        "Luminance entropy",             "structure",
)

load_photo_metrics <- function() {
  pm <- readr::read_csv(file.path(A23_DIR, "outcome/photo_metrics_per_id.csv"),
                        show_col_types = FALSE) %>%
    mutate(id = as.integer(id)) %>%
    # one row per id: photo-flagged clusters only, and if a duplicate remains
    # (defensive), keep the largest cluster
    filter(photo_flag == 1) %>%
    arrange(id, desc(n_shots)) %>%
    distinct(id, .keep_all = TRUE) %>%
    mutate(
      # vegetation hue deviation from green reference (120 deg): positive =
      # yellow-shifted (autumn browning); circular-safe within the 60-170 band
      veg_hue_dev = 120 - veg_hue,
      equip_proxy = case_when(n_shots >= 7 ~ "dual",
                              n_shots == 4 ~ "shots4",
                              n_shots <= 3 ~ "shots3")
    )
  pm
}

# photo-valid core joined with metrics; adds z-scores for the quality vars
load_quality_data <- function() {
  d <- load_vis_data() %>%
    inner_join(load_photo_metrics() %>%
                 select(id, panorama_group, n_shots, dual_generation,
                        equip_proxy, all_of(PHOTO_METRICS), veg_hue_dev),
               by = "id")
  for (v in c(QUALITY_VARS$var, "veg_frac", "sky_frac_top")) {
    d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
  }
  d %>%
    mutate(season = factor(sub(".* ", "", wave_label),
                           levels = c("autumn", "spring", "summer")))
}
