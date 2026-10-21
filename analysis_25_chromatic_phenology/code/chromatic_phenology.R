# A25 - chromatic phenology: does the photo colour layer carry the visual
# quality that the segmentation misses (A17/A18's inferred chromatic axis)?
# Identification follows A17: multi-season sites only (>= 2 seasons), models
# with (1|site) + (1|site:wave), equipment proxy as fixed technical covariate.
# Pre-scoping exploratory layer (raw-photo delivery 2026-08-18).
suppressPackageStartupMessages({library(tidyverse); library(lmerTest)})
source("analysis_23_photo_quality_metrics/code/photo_helpers.R")

ADIR <- "analysis_25_chromatic_phenology"
OUT <- file.path(ADIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20260818)

d <- load_quality_data() %>% mutate(equip = factor(equip_proxy))
season_tab <- table(as.character(d$site_label), d$season)
multi_sites <- rownames(season_tab)[rowSums(season_tab > 0) >= 2]
dm <- d %>% filter(site_label %in% multi_sites, season %in% c("autumn", "spring")) %>%
  droplevels()
cat("Multi-season sites:", paste(multi_sites, collapse = ", "), "\n")
cat("Rows (autumn/spring):", nrow(dm), "-"); print(table(dm$season))

# --- 1. quantity constancy vs chromatic movement (declared family: 6 metrics) ---
metrics <- tribble(
  ~var,          ~label,                          ~type,
  "comp_greenery", "Greenery share (quantity)",   "quantity",
  "veg_frac",    "Green-pixel fraction",          "chromatic",
  "veg_sat",     "Vegetation saturation",         "chromatic",
  "veg_hue_dev", "Vegetation hue yellow-shift",   "chromatic",
  "exg_mean",    "Excess green",                  "chromatic",
  "colorfulness", "Colourfulness",                "chromatic",
)
season_effects <- pmap_dfr(metrics, function(var, label, type) {
  fit <- lmer(reformulate(c("season", "equip", "(1 | site)", "(1 | site:wave)"), var), data = dm)
  co <- summary(fit)$coefficients
  tibble(var = var, label = label, type = type,
         est_spring_minus_autumn = co["seasonspring", "Estimate"],
         se = co["seasonspring", "Std. Error"], p = co["seasonspring", "Pr(>|t|)"],
         sd_var = sd(dm[[var]], na.rm = TRUE)) %>%
    mutate(std_est = est_spring_minus_autumn / sd_var)
})
season_effects <- season_effects %>% mutate(q = p.adjust(p, "BH"))
write_csv(season_effects, file.path(OUT, "season_effects_on_metrics.csv"))
cat("\nSeason effects on visual metrics (spring - autumn; BH-FDR over 6):\n")
print(season_effects %>% select(label, type, est_spring_minus_autumn, std_est, p, q), n = 10)

# --- 2. the matched-quantity chromatic contrast: colour season effect holding
#        segmented greenery constant (the direct quality-blindness statement) ---
fit_m <- lmer(veg_sat ~ season + comp_greenery + equip + (1 | site) + (1 | site:wave), data = dm)
co_m <- summary(fit_m)$coefficients
cat("\nveg_sat ~ season at matched greenery share:\n"); print(co_m["seasonspring", , drop = FALSE])
fit_m2 <- lmer(veg_frac ~ season + comp_greenery + equip + (1 | site) + (1 | site:wave), data = dm)
co_m2 <- summary(fit_m2)$coefficients
cat("veg_frac ~ season at matched greenery share:\n"); print(co_m2["seasonspring", , drop = FALSE])
write_csv(bind_rows(
  as_tibble(co_m["seasonspring", , drop = FALSE]) %>% mutate(outcome = "veg_sat"),
  as_tibble(co_m2["seasonspring", , drop = FALSE]) %>% mutate(outcome = "veg_frac")),
  file.path(OUT, "matched_quantity_contrasts.csv"))

# --- 3. does chromatic quality explain the spring-only greenery slope (A17)? ---
# A17: greenery -> OA/ISO-P positive in spring only. Test: (a) replicate the
# greenery x season interaction; (b) add veg quality (veg_sat, veg_frac) and its
# main effect - does the interaction attenuate?; (c) quality x quantity model.
dm2 <- dm %>% mutate(greenery_z = as.numeric(scale(comp_greenery)),
                     veg_sat_z2 = as.numeric(scale(veg_sat)),
                     veg_frac_z2 = as.numeric(scale(veg_frac)))
res3 <- map_dfr(c("overall_soundscape_z", "isop_z"), function(oc) {
  f0 <- lmer(reformulate(c("greenery_z * season", "equip", "(1 | site)", "(1 | site:wave)"), oc), data = dm2)
  c0 <- summary(f0)$coefficients
  f1 <- lmer(reformulate(c("greenery_z * season", "veg_sat_z2", "veg_frac_z2", "equip",
                           "(1 | site)", "(1 | site:wave)"), oc), data = dm2)
  c1 <- summary(f1)$coefficients
  f2 <- lmer(reformulate(c("greenery_z", "season", "veg_frac_z2", "veg_sat_z2", "greenery_z:veg_sat_z2",
                           "equip", "(1 | site)", "(1 | site:wave)"), oc), data = dm2)
  c2 <- summary(f2)$coefficients
  tibble(outcome = oc,
         interact_base = c0["greenery_z:seasonspring", "Estimate"],
         p_base = c0["greenery_z:seasonspring", "Pr(>|t|)"],
         interact_quality_adj = c1["greenery_z:seasonspring", "Estimate"],
         p_quality_adj = c1["greenery_z:seasonspring", "Pr(>|t|)"],
         veg_sat_main = c1["veg_sat_z2", "Estimate"], p_veg_sat = c1["veg_sat_z2", "Pr(>|t|)"],
         veg_frac_main = c1["veg_frac_z2", "Estimate"], p_veg_frac = c1["veg_frac_z2", "Pr(>|t|)"],
         qual_x_quant = c2["greenery_z:veg_sat_z2", "Estimate"],
         p_qual_x_quant = c2["greenery_z:veg_sat_z2", "Pr(>|t|)"])
})
write_csv(res3, file.path(OUT, "greenery_season_interaction_vs_quality.csv"))
cat("\nGreenery x season vs chromatic quality (OA, ISO-P):\n"); print(res3, width = Inf)

# --- 4. bootstrap CI (1,000, site-cluster) for the headline ratio contrast ---
gap <- dm %>% filter(comp_greenery > 0.05) %>%
  mutate(ratio = veg_frac / comp_greenery)
boot_est <- function(dat) {
  sites <- unique(dat$site_label)
  s <- sample(sites, length(sites), replace = TRUE)
  bd <- map_dfr(seq_along(s), function(i) dat %>% filter(site_label == s[i]) %>%
                  mutate(site_label = paste0(s[i], "_", i)))
  m <- tryCatch(lmer(ratio ~ season + equip + (1 | site_label), data = bd),
                error = function(e) NULL)
  if (is.null(m)) return(NA_real_)
  co <- summary(m)$coefficients
  if (!"seasonspring" %in% rownames(co)) return(NA_real_)
  co["seasonspring", "Estimate"]
}
boots <- replicate(1000, boot_est(gap))
ci <- quantile(boots, c(.025, .975), na.rm = TRUE)
cat("\nSite-cluster bootstrap (1,000) for spring-autumn ratio contrast: ",
    sprintf("est %.3f, 95%% CI [%.3f, %.3f], n_ok %d\n",
            mean(boots, na.rm = TRUE), ci[1], ci[2], sum(!is.na(boots))))
# point estimate on the full subset (the bootstrap gives the interval)
fit_gap <- lmer(ratio ~ season + equip + (1 | site_label), data = gap)
est_gap <- summary(fit_gap)$coefficients["seasonspring", "Estimate"]
cat(sprintf("Point estimate on the full subset: %.3f\n", est_gap))
write_csv(tibble(est = est_gap, mean = mean(boots, na.rm = TRUE), lo = ci[1], hi = ci[2],
                 n_ok = sum(!is.na(boots))),
          file.path(OUT, "bootstrap_ratio_contrast.csv"))
