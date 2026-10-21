# Analysis 12 — Robustness battery for the project-22 findings:
#   V1 vehicles -> ISO-E (main)          V2 traffic x visible traffic -> OA
#   V3 water_any -> OA                   V4 flower_any -> OA
# Checks: LOSO (12 sites), cluster bootstrap (1,000 site resamples), wave
# dropped vs kept, renormalised-composite variant (V1), Mundlak within-only (V1).

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_12_robustness"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(comp_traffic_vis = comp_vehicles + ph_road,
         comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
         water_any = as.integer(comp_water > 0.005),
         flower_any = as.integer(replace_na(ph_flower, 0) > 0.001))

EFFECTS <- list(
  V1_vehicles_isoe = list(
    f = "isoe_z ~ comp_vehicles_z + wave + (1 | site)", term = "comp_vehicles_z"),
  V2_traffic_congruence_oa = list(
    f = "overall_soundscape_z ~ traffic_noise_z * comp_traffic_vis_z + wave + (1 | site)",
    term = "traffic_noise_z:comp_traffic_vis_z"),
  V3_water_oa = list(
    f = "overall_soundscape_z ~ water_any + wave + (1 | site)", term = "water_any"),
  V4_flower_oa = list(
    f = "overall_soundscape_z ~ flower_any + wave + (1 | site)", term = "flower_any")
)

fit_term <- function(f, dd, term) {
  fit <- lmer(as.formula(f), data = dd)
  tidy_term(fit, term)
}

main <- imap_dfr(EFFECTS, function(spec, name) {
  bind_rows(
    fit_term(spec$f, d, spec$term) %>% mutate(check = "baseline"),
    fit_term(str_replace(spec$f, "\\+ wave ", ""), d, spec$term) %>%
      mutate(check = "no_wave")) %>%
    mutate(effect = name, .before = 1)
})
write_csv(main, file.path(OUT, "battery_main.csv"))

# V1 variants: renormalised composite; Mundlak within-only
v1_rn <- fit_term("isoe_z ~ comp_vehicles_rn_z + wave + (1 | site)",
                  d, "comp_vehicles_rn_z") %>% mutate(check = "renormalised")
dm <- d %>% group_by(site) %>%
  mutate(veh_b = mean(comp_vehicles_z), veh_w = comp_vehicles_z - veh_b) %>%
  ungroup()
v1_mund <- lmer(isoe_z ~ veh_w + veh_b + wave + (1 | site), data = dm)
v1_extra <- bind_rows(
  v1_rn,
  broom.mixed::tidy(v1_mund, effects = "fixed") %>%
    filter(term %in% c("veh_w", "veh_b")) %>%
    mutate(check = paste0("mundlak_", term)) %>%
    select(term, estimate, std.error, statistic, df, p.value, check))
write_csv(v1_extra %>% mutate(effect = "V1_vehicles_isoe", .before = 1),
          file.path(OUT, "v1_variants.csv"))

# LOSO
loso <- imap_dfr(EFFECTS, function(spec, name) {
  map_dfr(unique(as.character(d$site)), function(s) {
    est <- tryCatch(fit_term(spec$f, filter(d, site != s), spec$term),
                    error = function(e) tibble(estimate = NA, p.value = NA))
    tibble(effect = name, dropped_site = s, estimate = est$estimate[1],
           p.value = est$p.value[1])
  })
})
write_csv(loso, file.path(OUT, "loso.csv"))
loso_sum <- loso %>% group_by(effect) %>%
  summarise(min_est = min(estimate, na.rm = TRUE),
            max_est = max(estimate, na.rm = TRUE),
            sign_flips = sum(sign(estimate) != sign(median(estimate, na.rm = TRUE)),
                             na.rm = TRUE),
            max_p = max(p.value, na.rm = TRUE))
write_csv(loso_sum, file.path(OUT, "loso_summary.csv"))

# Cluster bootstrap (1,000)
boot <- imap_dfr(EFFECTS, function(spec, name) {
  sites <- unique(as.character(d$site))
  ests <- map_dbl(seq_len(1000), function(b) {
    samp <- sample(sites, length(sites), replace = TRUE)
    db <- map_dfr(seq_along(samp), function(i)
      d %>% filter(site == samp[i]) %>% mutate(site_b = paste0("s", i)))
    fit <- tryCatch(lmer(as.formula(str_replace(spec$f, "site\\)", "site_b)")),
                         data = db), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(term == spec$term)
    if (nrow(tt)) tt$estimate[1] else NA_real_
  })
  es <- ests[!is.na(ests)]
  tibble(effect = name, n_boot = length(es),
         ci_lo = quantile(es, 0.025), ci_hi = quantile(es, 0.975),
         excludes_zero = sign(quantile(es, 0.025)) == sign(quantile(es, 0.975)))
})
write_csv(boot, file.path(OUT, "cluster_bootstrap.csv"))

print(main %>% select(effect, check, estimate, p.value))
print(v1_extra %>% select(check, estimate, p.value))
print(loso_sum); print(boot)
