# A26 - do photo-measured visual-quality metrics predict soundscape evaluation
# beyond the segmentation quantity layer?
# Declared family: 9 quality metrics x 4 outcomes (OA, ISO-P, ISO-E, harmony)
# = 36 tests, BH-FDR across the family. Models: outcome_z ~ metric_z +
# quantity control (per metric block) + wave + equip + (1|site) + (1|panorama_group).
# Survivors: site-cluster bootstrap (1,000) + LOSO. Sensitivity: exclude
# blur bottom-decile and clip_hi top-decile rows.
# Pre-scoping exploratory layer (raw-photo delivery 2026-08-18).
suppressPackageStartupMessages({library(tidyverse); library(lmerTest)})
source("analysis_23_photo_quality_metrics/code/photo_helpers.R")

ADIR <- "analysis_26_quality_outcome_scan"
OUT <- file.path(ADIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20260818)

d <- load_quality_data() %>%
  mutate(equip = factor(equip_proxy),
         veg_sat_z = as.numeric(scale(veg_sat)),
         veg_hue_dev_z = as.numeric(scale(veg_hue_dev)))

OUTCOMES4 <- c("overall_soundscape_z", "isop_z", "isoe_z", "harmony_z")
CONTROLS <- list(chromatic = "comp_greenery_z", sky = "comp_sky_z",
                 luminance = NULL, structure = NULL)

scan <- crossing(QUALITY_VARS, outcome = OUTCOMES4) %>%
  pmap_dfr(function(var, label, block, outcome) {
    ctrl <- CONTROLS[[block]]
    rhs <- c(paste0(var, "_z"), ctrl, "wave", "equip", "(1 | site)", "(1 | panorama_group)")
    dd <- d %>% drop_na(all_of(c(var, outcome)))
    fit <- tryCatch(lmer(reformulate(rhs, outcome), data = dd), error = function(e) NULL)
    if (is.null(fit)) return(tibble(var = var, outcome = outcome, est = NA_real_))
    co <- summary(fit)$coefficients
    tibble(var = var, label = label, block = block, outcome = outcome, n = nrow(dd),
           est = co[paste0(var, "_z"), "Estimate"], se = co[paste0(var, "_z"), "Std. Error"],
           p = co[paste0(var, "_z"), "Pr(>|t|)"])
  }) %>%
  mutate(q = p.adjust(p, "BH"))
write_csv(scan, file.path(OUT, "quality_outcome_scan.csv"))
cat("A26 scan (family n =", nrow(scan), "), FDR survivors q < 0.05:\n")
surv <- scan %>% filter(q < 0.05) %>% arrange(q)
print(surv %>% select(label, outcome, est, p, q, n), n = 40)
cat("\nNominal only (p < 0.05, q >= 0.05):\n")
print(scan %>% filter(p < 0.05, q >= 0.05) %>% arrange(p) %>%
        select(label, outcome, est, p, q), n = 40)

# --- joint quality-block test per outcome (LRT vs quantity-only baseline) ---
qvars_z <- paste0(QUALITY_VARS$var, "_z")
joint <- map_dfr(OUTCOMES4, function(oc) {
  dd <- d %>% drop_na(all_of(c(QUALITY_VARS$var, "comp_greenery", "comp_sky"))) %>%
    drop_na(all_of(sub("_z$", "", oc)))
  base <- lmer(reformulate(c("comp_greenery_z", "comp_sky_z", "wave", "equip",
                             "(1 | site)", "(1 | panorama_group)"), oc),
               data = dd, REML = FALSE)
  full <- lmer(reformulate(c("comp_greenery_z", "comp_sky_z", qvars_z, "wave", "equip",
                             "(1 | site)", "(1 | panorama_group)"), oc),
               data = dd, REML = FALSE)
  lrt <- anova(base, full)
  tibble(outcome = oc, n = nrow(dd), chisq = lrt$Chisq[2], df = lrt$Df[2],
         p_joint = lrt$`Pr(>Chisq)`[2],
         r2_marg_gain = performance::r2(full)$R2_marginal - performance::r2(base)$R2_marginal)
})
write_csv(joint, file.path(OUT, "joint_block_tests.csv"))
cat("\nJoint quality block (LRT vs quantity baseline):\n"); print(joint)

# --- survivor robustness: site-cluster bootstrap (1,000) + LOSO + blur/clip sensitivity ---
if (nrow(surv) > 0) {
  rob <- pmap_dfr(surv %>% select(var, block, outcome), function(var, block, outcome) {
    ctrl <- CONTROLS[[block]]
    rhs <- c(paste0(var, "_z"), ctrl, "wave", "equip", "(1 | site)", "(1 | panorama_group)")
    dd <- d %>% drop_na(all_of(c(var, sub("_z$", "", outcome))))
    # LOSO
    loso <- map_dbl(unique(as.character(dd$site_label)), function(s) {
      f <- tryCatch(lmer(reformulate(rhs, outcome), data = dd %>% filter(site_label != s)),
                    error = function(e) NULL)
      if (is.null(f)) return(NA_real_)
      summary(f)$coefficients[paste0(var, "_z"), "Estimate"]
    })
    # blur/clip sensitivity
    dd2 <- dd %>% filter(blur_lapvar > quantile(blur_lapvar, .1),
                         clip_hi < quantile(clip_hi, .9))
    f2 <- lmer(reformulate(rhs, outcome), data = dd2)
    co2 <- summary(f2)$coefficients
    # site-cluster bootstrap
    boots <- replicate(1000, {
      sites <- unique(as.character(dd$site_label))
      s <- sample(sites, length(sites), replace = TRUE)
      bd <- map_dfr(seq_along(s), function(i) dd %>% filter(site_label == s[i]) %>%
                      mutate(site = paste0(s[i], "_", i)))
      f <- tryCatch(suppressMessages(lmer(reformulate(rhs, outcome), data = bd)),
                    error = function(e) NULL)
      if (is.null(f)) NA_real_ else {
        co <- summary(f)$coefficients
        if (paste0(var, "_z") %in% rownames(co)) co[paste0(var, "_z"), "Estimate"] else NA_real_
      }
    })
    ci <- quantile(boots, c(.025, .975), na.rm = TRUE)
    tibble(var = var, outcome = outcome,
           loso_min = min(loso, na.rm = TRUE), loso_max = max(loso, na.rm = TRUE),
           loso_sign_stable = all(sign(loso[!is.na(loso)]) == sign(loso[which(!is.na(loso))[1]])),
           est_qualityfiltered = co2[paste0(var, "_z"), "Estimate"],
           p_qualityfiltered = co2[paste0(var, "_z"), "Pr(>|t|)"],
           boot_lo = ci[1], boot_hi = ci[2], boot_n_ok = sum(!is.na(boots)))
  })
  write_csv(rob, file.path(OUT, "survivor_robustness.csv"))
  cat("\nSurvivor robustness:\n"); print(rob, width = Inf)
}
