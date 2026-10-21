# Analysis 31 - Panorama-group clustering sensitivity (post-review robustness).
# Photographs are one-per-measurement-point, so respondents sharing a panorama
# share their visual predictors exactly; the condition-metric models (A26/A29)
# already carry a (1|panorama_group) term but the composition/element/congruence
# headline models (A05/A06/A11/A15/A19) carry only wave + (1|site) (or
# site + site:wave for A19). Here every Table-2 headline effect outside the
# condition block is refitted with the panorama-group random effect added, and
# the original specification is refitted alongside as the reproduction check.
# No new tests: these are robustness refits of already-declared FDR survivors.
# Panorama groups follow A28's assignment (A23 metrics table; unmatched rows
# become singleton groups).
# Also writes the multi-season-subset green-per-greenery means by season
# (A24 ratio definition on A25's identification frame) for Figure 3c.

source("_shared/helpers.R")
set.seed(20260819)

ANALYSIS_DIR <- "analysis_31_panorama_clustering"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- data with panorama groups (A28's assignment) ---------------------------
pg <- readr::read_csv("analysis_23_photo_quality_metrics/outcome/photo_metrics_per_id.csv",
      show_col_types = FALSE) %>%
  filter(photo_flag == 1) %>% arrange(id, desc(n_shots)) %>%
  distinct(id, .keep_all = TRUE) %>% select(id, panorama_group)

d <- load_vis_data() %>%
  left_join(pg, by = "id") %>%
  mutate(panorama_group = ifelse(is.na(panorama_group), paste0("solo_", id),
                                 panorama_group),
         # A11/A19 element definitions
         water_any = as.integer(comp_water > 0.005),
         flower_any = as.integer(replace_na(ph_flower, 0) > 0.001),
         # A06/A15 constructed variables (identical construction and scaling)
         comp_traffic_vis = comp_vehicles + ph_road,
         comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
         audio_nature_balance = natural_sounds - (traffic_noise + other_noise) / 2,
         visual_nature_balance = comp_nature - comp_traffic_vis,
         audio_nature_balance_z = as.numeric(scale(audio_nature_balance)),
         visual_nature_balance_z = as.numeric(scale(visual_nature_balance)),
         idx_balance_match = audio_nature_balance_z * visual_nature_balance_z,
         idx_balance_match_z = as.numeric(scale(idx_balance_match)))
cat("n =", nrow(d), "| panorama groups:", n_distinct(d$panorama_group), "\n")

da <- d %>% filter(!is.na(laeq_a))
cat("acoustic subsample n =", nrow(da), "\n")

# ---- the headline refits ----------------------------------------------------
# Each spec: reported model (reproduction check) and the same + (1|panorama_group).
SPECS <- tribble(
  ~effect, ~data, ~formula_base, ~focal,
  "vehicles_isoe (A05)", "d",
    "isoe_z ~ comp_vehicles_z + wave + (1 | site)", "comp_vehicles_z",
  "water_oa (A11)", "d",
    "overall_soundscape_z ~ water_any + wave + (1 | site)", "water_any",
  "water_isop (A11)", "d",
    "isop_z ~ water_any + wave + (1 | site)", "water_any",
  "flowers_oa (A11)", "d",
    "overall_soundscape_z ~ flower_any + wave + (1 | site)", "flower_any",
  "traffic_congruence_oa (A06)", "d",
    "overall_soundscape_z ~ traffic_noise_z * comp_traffic_vis_z + wave + (1 | site)",
    "traffic_noise_z:comp_traffic_vis_z",
  "balance_match_isoe (A15)", "d",
    paste("isoe_z ~ idx_balance_match_z + audio_nature_balance_z +",
          "visual_nature_balance_z + wave + (1 | site)"), "idx_balance_match_z",
  "water_x_pl_oa (A19)", "da",
    "overall_soundscape_z ~ water_any * perceived_loudness_z + (1 | site) + (1 | site:wave)",
    "water_any:perceived_loudness_z")

fit_one <- function(f_str, dat, focal) {
  fit <- suppressMessages(lmer(as.formula(f_str), data = dat))
  tt <- broom.mixed::tidy(fit, effects = "fixed") %>% filter(term == focal)
  tibble(estimate = tt$estimate, se = tt$std.error, p.value = tt$p.value,
         n = nobs(fit), singular = isSingular(fit))
}

res <- pmap_dfr(SPECS, function(effect, data, formula_base, focal) {
  dat <- if (data == "da") da else d
  orig <- fit_one(formula_base, dat, focal) %>% mutate(spec = "as reported")
  pano <- fit_one(paste(formula_base, "+ (1 | panorama_group)"), dat, focal) %>%
    mutate(spec = "plus panorama RE")
  bind_rows(orig, pano) %>% mutate(effect = effect, .before = 1)
})
res <- res %>% select(effect, spec, estimate, se, p.value, n, singular)
write_csv(res, file.path(OUT, "panorama_sensitivity.csv"))
print(res, n = 20, width = Inf)

# ---- Figure 3c support: subset green-per-greenery means by season -----------
# A24's ratio (green pixels per unit segmented greenery, rows with greenery
# share > 0.05) on A25's identification frame (multi-season sites, autumn/
# spring), with A25's site-cluster bootstrap contrast alongside.
source("analysis_23_photo_quality_metrics/code/photo_helpers.R")
q <- load_quality_data()
season_tab <- table(as.character(q$site_label), q$season)
multi_sites <- rownames(season_tab)[rowSums(season_tab > 0) >= 2]
gap <- q %>%
  filter(site_label %in% multi_sites, season %in% c("autumn", "spring"),
         comp_greenery > 0.05) %>%
  mutate(ratio = veg_frac / comp_greenery)
gap_season <- gap %>% group_by(season) %>%
  summarise(n = n(), mean_ratio = mean(ratio), sd = sd(ratio), .groups = "drop")
write_csv(gap_season, file.path(OUT, "quality_gap_by_season_subset.csv"))
cat("\nSubset green-per-greenery by season:\n"); print(gap_season)
