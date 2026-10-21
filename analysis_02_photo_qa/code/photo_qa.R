# Analysis 02 — Photo QA: validity flags, row sums (the 54 > 1.05 decision),
# zero rows, Fangta handling, batch (wave) structure.

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_02_photo_qa"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d_all <- load_rest_data("core")   # incl. photo-invalid rows, for the flag audit
d <- load_vis_data()

# --- Flag audit -----------------------------------------------------------------
flags <- d_all %>%
  mutate(photo_state = case_when(
    photo_ok == 1 & photo_present == 1 ~ "valid",
    photo_ok == 1 & photo_present == 0 ~ "flag_ok_but_no_data",
    photo_ok == 0 & photo_present == 1 ~ "data_but_flag_0",
    TRUE ~ "no_photo")) %>%
  count(site_label, photo_state) %>%
  pivot_wider(names_from = photo_state, values_from = n, values_fill = 0)
write_csv(flags, file.path(OUT, "photo_flag_audit_by_site.csv"))

# --- Row-sum distribution ---------------------------------------------------------
rs <- d %>%
  summarise(n = n(), min = min(photo_row_sum), q05 = quantile(photo_row_sum, .05),
            median = median(photo_row_sum), q95 = quantile(photo_row_sum, .95),
            max = max(photo_row_sum),
            n_below_095 = sum(photo_row_sum < 0.95),
            n_over_105 = sum(photo_row_sum > 1.05))
write_csv(rs, file.path(OUT, "row_sum_distribution.csv"))
write_csv(d %>% filter(photo_row_sum > 1.05) %>%
            select(id, site_label, wave_label, photo_row_sum) %>% arrange(id),
          file.path(OUT, "rows_over_105.csv"))

# --- Renormalisation decision: does it matter? ------------------------------------
# Correlate raw vs renormalised composites; refit the visual->isop screen both ways.
cor_tab <- map_dfr(COMPS8, function(v) {
  tibble(composite = v,
         r_raw_vs_rn = cor(d[[v]], d[[paste0(v, "_rn")]], use = "pairwise"))
})
write_csv(cor_tab, file.path(OUT, "raw_vs_renormalised_correlations.csv"))

both <- map_dfr(COMPS8, function(v) {
  f1 <- lmer(as.formula(paste0("isop_z ~ ", v, "_z + (1 | site)")), data = d)
  f2 <- lmer(as.formula(paste0("isop_z ~ ", v, "_rn_z + (1 | site)")), data = d)
  bind_rows(tidy_term(f1, paste0(v, "_z")) %>% mutate(variant = "raw"),
            tidy_term(f2, paste0(v, "_rn_z")) %>% mutate(variant = "renormalised")) %>%
    mutate(composite = v, .before = 1)
})
write_csv(both, file.path(OUT, "isop_screen_raw_vs_renormalised.csv"))

# --- Batch (wave) structure ---------------------------------------------------------
batch <- d %>%
  group_by(site_label, wave_label) %>%
  summarise(n = n(), row_sum_mean = mean(photo_row_sum),
            greenery = mean(comp_greenery), sky = mean(comp_sky),
            people = mean(comp_people), .groups = "drop")
write_csv(batch, file.path(OUT, "batch_structure.csv"))

# Wave effect on composites within site (is there a batch artefact?)
bt <- map_dfr(c("comp_greenery", "comp_sky", "comp_people", "comp_built"),
              function(v) {
  dd <- d %>% group_by(site_label) %>% filter(n_distinct(wave_code) > 1) %>% ungroup()
  fit0 <- lm(as.formula(paste0(v, " ~ site_label")), data = dd)
  fit1 <- lm(as.formula(paste0(v, " ~ site_label + factor(wave_code)")), data = dd)
  tibble(composite = v, r2_site = summary(fit0)$r.squared,
         r2_site_wave = summary(fit1)$r.squared,
         wave_partial_r2 = (summary(fit1)$r.squared - summary(fit0)$r.squared) /
           (1 - summary(fit0)$r.squared),
         anova_p = anova(fit0, fit1)$`Pr(>F)`[2])
})
write_csv(bt, file.path(OUT, "wave_batch_effects.csv"))

print(rs); print(cor_tab, n = 8); print(bt)
