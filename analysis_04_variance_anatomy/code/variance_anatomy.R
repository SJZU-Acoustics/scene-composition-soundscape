# Analysis 04 — Variance anatomy gate: where does visual-composite variance
# live? ICC(site), wave share, and the usable within-site remainder.
# With 13 (here 12) sites, composites that vary only between sites are
# near-unidentifiable against site effects; this gate decides the estimand.

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_04_variance_anatomy"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data()

anatomy <- map_dfr(COMPS8, function(v) {
  m0 <- lmer(as.formula(paste0(v, "_z ~ (1 | site)")), data = d)
  vc0 <- as.data.frame(VarCorr(m0))
  icc_site <- vc0$vcov[1] / sum(vc0$vcov)
  m1 <- lmer(as.formula(paste0(v, "_z ~ (1 | site) + (1 | site:wave)")), data = d)
  vc1 <- as.data.frame(VarCorr(m1))
  v_sw <- vc1$vcov[match("site:wave", vc1$grp)]
  v_s <- vc1$vcov[match("site", vc1$grp)]
  v_r <- vc1$vcov[vc1$grp == "Residual"]
  tibble(composite = v, icc_site = icc_site,
         var_site = v_s / (v_s + v_sw + v_r),
         var_site_wave = v_sw / (v_s + v_sw + v_r),
         var_within = v_r / (v_s + v_sw + v_r),
         sd_within_raw = sd(resid(lm(as.formula(paste0(v, " ~ site + wave")),
                                     data = d))))
})
write_csv(anatomy, file.path(OUT, "variance_anatomy.csv"))

# Within-site distributional reach: per-site SD and range of each composite
reach <- map_dfr(COMPS8, function(v) {
  d %>% group_by(site_label) %>%
    summarise(sd = sd(.data[[v]]), iqr = IQR(.data[[v]]), .groups = "drop") %>%
    summarise(composite = v, median_within_site_sd = median(sd),
              min_sd = min(sd), max_sd = max(sd),
              median_iqr = median(iqr))
})
write_csv(reach, file.path(OUT, "within_site_reach.csv"))

print(anatomy, n = 8)
print(reach, n = 8)
