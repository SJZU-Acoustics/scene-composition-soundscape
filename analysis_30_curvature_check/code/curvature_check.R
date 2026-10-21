# A30 - curvature check on the quantity boundary. The boundary claim (A05/A08:
# within-site scene composition barely moves evaluation) rests on LINEAR slopes.
# A reviewer can object that a non-monotone dose (e.g. greenery good only past a
# threshold, or an inverted-U) is invisible to a linear scan. Declared family:
# 8 composites x 4 outcomes = 32 quadratic terms, BH-FDR. Models:
# outcome_z ~ comp_z + comp_z^2 + wave + (1|site). Survivors (if any) get a
# sound-battery-adjusted refit and a site-cluster bootstrap (1,000).
source("_shared/helpers.R")
set.seed(20260819)

ANALYSIS_DIR <- "analysis_30_curvature_check"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- load_vis_data()
CURV_OUTS <- c(overall_soundscape_z = "OA", isop_z = "ISO-P",
               isoe_z = "ISO-E", harmony_z = "Harmony")
SOUND <- c("traffic_noise_z", "human_sounds_z", "natural_sounds_z",
           "other_noise_z", "perceived_loudness_z")

scan <- crossing(comp = COMPS8, outcome = names(CURV_OUTS)) %>%
  pmap_dfr(function(comp, outcome) {
    zcol <- paste0(comp, "_z")
    dd <- d %>% mutate(cx = .data[[zcol]], cx2 = cx^2) %>%
      drop_na(all_of(c(outcome, zcol)))
    fit <- lmer(reformulate(c("cx", "cx2", "wave", "(1 | site)"), outcome), data = dd)
    co <- summary(fit)$coefficients
    tibble(composite = comp, outcome = outcome, n = nrow(dd),
           b_linear = co["cx", "Estimate"], p_linear = co["cx", "Pr(>|t|)"],
           b_quad = co["cx2", "Estimate"], se_quad = co["cx2", "Std. Error"],
           p_quad = co["cx2", "Pr(>|t|)"])
  }) %>%
  mutate(outcome_label = CURV_OUTS[outcome], q_quad = fdr(p_quad)) %>%
  arrange(q_quad)
write_csv(scan, file.path(OUT, "curvature_scan.csv"))
cat("A30 curvature scan (family =", nrow(scan), "); quadratic survivors q < 0.05:\n")
print(scan %>% select(composite, outcome_label, b_linear, p_linear, b_quad, p_quad, q_quad), n = 40)

# --- survivor follow-up: sound-adjusted refit + site-cluster bootstrap --------
surv <- scan %>% filter(q_quad < 0.05)
if (nrow(surv) > 0) {
  follow <- pmap_dfr(surv %>% select(composite, outcome), function(composite, outcome) {
    zcol <- paste0(composite, "_z")
    dd <- d %>% mutate(cx = .data[[zcol]], cx2 = cx^2) %>%
      drop_na(all_of(c(outcome, zcol)))
    fit <- lmer(reformulate(c("cx", "cx2", SOUND, "wave", "(1 | site)"), outcome),
                data = dd)
    co <- summary(fit)$coefficients
    sites <- unique(as.character(dd$site_label))
    boots <- replicate(1000, {
      s <- sample(sites, length(sites), replace = TRUE)
      bd <- map_dfr(seq_along(s), function(i) dd %>% filter(site_label == s[i]) %>%
                      mutate(site = paste0(s[i], "_", i)))
      f <- tryCatch(suppressMessages(lmer(
        reformulate(c("cx", "cx2", "wave", "(1 | site)"), outcome), data = bd)),
        error = function(e) NULL)
      if (is.null(f)) return(NA_real_)
      cc <- summary(f)$coefficients
      if ("cx2" %in% rownames(cc)) cc["cx2", "Estimate"] else NA_real_
    })
    ci <- quantile(boots, c(.025, .975), na.rm = TRUE)
    tibble(composite = composite, outcome = outcome,
           b_quad_sound_adj = co["cx2", "Estimate"],
           p_quad_sound_adj = co["cx2", "Pr(>|t|)"],
           boot_lo = ci[1], boot_hi = ci[2], boot_n_ok = sum(!is.na(boots)))
  })
  write_csv(follow, file.path(OUT, "curvature_survivor_followup.csv"))
  cat("\nSurvivor follow-up:\n"); print(follow, width = Inf)
} else {
  cat("\nNo curvature survivors: the linear-scan boundary is not a functional-form artefact.\n")
}
