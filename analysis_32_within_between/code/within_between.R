# Analysis 32 - Within-site / between-site decomposition of the Table-2 effects
# (post-review robustness, 2026-09-18). A12 ran the decomposition for visible
# vehicles -> ISO-E only; here every reported effect is refitted with each
# model term split into its site mean (between-site component) and the
# respondent's deviation from that mean (within-site component). For
# interaction and product terms the constituents and the product are all split,
# so the within-site coefficient of the product equals its site-fixed-effects
# estimate. No new tests: robustness refits of effects already significant in
# their own families.

source("analysis_23_photo_quality_metrics/code/photo_helpers.R")

ANALYSIS_DIR <- "analysis_32_within_between"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- data: identical constructions to A05/A06/A11/A15/A19/A26 ----------------
build <- function(dat) dat %>%
  mutate(water_any = as.integer(comp_water > 0.005),
         flower_any = as.integer(replace_na(ph_flower, 0) > 0.001),
         comp_traffic_vis = comp_vehicles + ph_road,
         comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
         audio_nature_balance = natural_sounds - (traffic_noise + other_noise) / 2,
         visual_nature_balance = comp_nature - comp_traffic_vis,
         audio_nature_balance_z = as.numeric(scale(audio_nature_balance)),
         visual_nature_balance_z = as.numeric(scale(visual_nature_balance)),
         idx_balance_match_z = as.numeric(scale(audio_nature_balance_z * visual_nature_balance_z)),
         traffic_congruence = traffic_noise_z * comp_traffic_vis_z,
         water_x_pl = water_any * perceived_loudness_z)
d <- build(load_vis_data())
da <- d %>% filter(!is.na(laeq_a))
dq <- build(load_quality_data()) %>% mutate(equip = factor(equip_proxy))
cat("n =", nrow(d), "| acoustic n =", nrow(da), "| condition n =", nrow(dq), "\n")

# ---- specifications: terms to split, focal term, remaining right-hand side ---
SPECS <- tribble(
  ~effect, ~data, ~outcome, ~split, ~focal, ~rest,
  "vehicles_isoe", "d", "isoe_z", "comp_vehicles_z", "comp_vehicles_z",
    "wave + (1 | site)",
  "water_oa", "d", "overall_soundscape_z", "water_any", "water_any",
    "wave + (1 | site)",
  "water_isop", "d", "isop_z", "water_any", "water_any",
    "wave + (1 | site)",
  "flowers_oa", "d", "overall_soundscape_z", "flower_any", "flower_any",
    "wave + (1 | site)",
  "traffic_congruence_oa", "d", "overall_soundscape_z",
    "traffic_noise_z comp_traffic_vis_z traffic_congruence", "traffic_congruence",
    "wave + (1 | site)",
  "balance_match_isoe", "d", "isoe_z",
    "idx_balance_match_z audio_nature_balance_z visual_nature_balance_z", "idx_balance_match_z",
    "wave + (1 | site)",
  "blueness_oa", "dq", "overall_soundscape_z", "sky_blue_ratio_z comp_sky_z", "sky_blue_ratio_z",
    "wave + equip + (1 | site) + (1 | panorama_group)",
  "water_x_pl_oa", "da", "overall_soundscape_z",
    "water_any perceived_loudness_z water_x_pl", "water_x_pl",
    "(1 | site) + (1 | site:wave)")

res <- pmap_dfr(SPECS, function(effect, data, outcome, split, focal, rest) {
  effect_name <- effect
  dat <- get(data); terms <- str_split(split, " ")[[1]]
  for (v in terms) {
    dat <- dat %>% group_by(site) %>%
      mutate("{v}_b" := mean(.data[[v]]), "{v}_w" := .data[[v]] - mean(.data[[v]])) %>%
      ungroup()
  }
  rhs <- paste(c(paste0(terms, "_w"), paste0(terms, "_b"), rest), collapse = " + ")
  fit <- suppressMessages(lmer(as.formula(paste(outcome, "~", rhs)), data = dat))
  broom.mixed::tidy(fit, effects = "fixed") %>%
    filter(term %in% paste0(focal, c("_w", "_b"))) %>%
    transmute(effect = effect_name,
              component = ifelse(str_ends(term, "_w"), "within-site", "between-site"),
              estimate, se = std.error, p.value, n = nobs(fit), singular = isSingular(fit))
})
write_csv(res, file.path(OUT, "within_between.csv"))
print(res, n = 20, width = Inf)
