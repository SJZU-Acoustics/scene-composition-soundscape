# Analysis 14 — Opportunity-conditioned visual effects.
# Instead of asking every site to test every visual mechanism, condition the
# main visual effects on whether the relevant visual feature actually varies
# within the site. Coefficients keep the full-sample z scales from helpers.R.

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_14_opportunity_conditioned"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(
    comp_traffic_vis = comp_vehicles + ph_road,
    comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
    water_any = as.integer(comp_water > 0.005),
    flower_any = as.integer(replace_na(ph_flower, 0) > 0.001)
  )

site_opp <- d %>% group_by(site, site_label, space_class) %>%
  summarise(
    n = n(),
    veh_sd = sd(comp_vehicles), veh_prev = mean(comp_vehicles > 0.005),
    veh_max = max(comp_vehicles),
    traffic_vis_sd = sd(comp_traffic_vis),
    traffic_vis_prev = mean(comp_traffic_vis > 0.02),
    traffic_vis_max = max(comp_traffic_vis),
    water_sd = sd(comp_water), water_prev = mean(water_any),
    water_max = max(comp_water),
    flower_sd = sd(replace_na(ph_flower, 0)), flower_prev = mean(flower_any),
    flower_max = max(replace_na(ph_flower, 0)),
    .groups = "drop"
  ) %>%
  mutate(
    opp_vehicle = veh_prev >= 0.25 & veh_sd >= 0.007,
    opp_traffic_vis = traffic_vis_prev >= 0.30 & traffic_vis_sd >= 0.045,
    opp_water = water_prev >= 0.20 & water_sd >= 0.015,
    opp_flower = flower_prev >= 0.10 | flower_sd >= 0.005
  )
write_csv(site_opp, file.path(OUT, "site_opportunity_table.csv"))

d <- d %>% left_join(site_opp %>%
                      select(site, opp_vehicle, opp_traffic_vis,
                             opp_water, opp_flower),
                    by = "site")

EFFECTS <- tribble(
  ~effect, ~outcome, ~formula, ~fixed_formula, ~term, ~opp_flag,
  "vehicles_to_eventfulness", "isoe",
  "isoe_z ~ comp_vehicles_z + wave + (1 | site)",
  "isoe_z ~ comp_vehicles_z + wave + site", "comp_vehicles_z",
  "opp_vehicle",
  "traffic_seen_heard_to_oa", "overall_soundscape",
  "overall_soundscape_z ~ traffic_noise_z * comp_traffic_vis_z + wave + (1 | site)",
  "overall_soundscape_z ~ traffic_noise_z * comp_traffic_vis_z + wave + site",
  "traffic_noise_z:comp_traffic_vis_z", "opp_traffic_vis",
  "water_to_oa", "overall_soundscape",
  "overall_soundscape_z ~ water_any + wave + (1 | site)",
  "overall_soundscape_z ~ water_any + wave + site", "water_any",
  "opp_water",
  "flower_to_oa", "overall_soundscape",
  "overall_soundscape_z ~ flower_any + wave + (1 | site)",
  "overall_soundscape_z ~ flower_any + wave + site", "flower_any",
  "opp_flower"
)

fit_one <- function(f, term, dd) {
  fit <- suppressMessages(suppressWarnings(lmer(as.formula(f), data = dd)))
  tidy_term(fit, term)
}

models <- pmap_dfr(EFFECTS, function(effect, outcome, formula, fixed_formula,
                                     term, opp_flag) {
  dd_opp <- d %>% filter(.data[[opp_flag]])
  dd_low <- d %>% filter(!.data[[opp_flag]])
  fit_all <- suppressMessages(suppressWarnings(lmer(as.formula(formula), data = d)))
  fit_opp <- suppressMessages(suppressWarnings(lmer(as.formula(formula), data = dd_opp)))
  bind_rows(
    tidy_term(fit_all, term) %>%
      mutate(sample = "all", n = nobs(fit_all),
             n_sites = n_distinct(d$site)),
    tidy_term(fit_opp, term) %>%
      mutate(sample = "opportunity_sites", n = nobs(fit_opp),
             n_sites = n_distinct(dd_opp$site)),
    if (n_distinct(dd_low$site) >= 2) {
      fit_low <- suppressMessages(suppressWarnings(lmer(as.formula(formula), data = dd_low)))
      tidy_term(fit_low, term) %>%
        mutate(sample = "low_opportunity_sites",
               n = nobs(fit_low),
               n_sites = n_distinct(dd_low$site))
    } else {
      tibble()
    }
  ) %>%
    mutate(effect = effect, outcome = outcome, opportunity_flag = opp_flag,
           .before = 1)
})
write_csv(models, file.path(OUT, "opportunity_models.csv"))

loso <- pmap_dfr(EFFECTS, function(effect, outcome, formula, fixed_formula,
                                   term, opp_flag) {
  dd <- d %>% filter(.data[[opp_flag]])
  map_dfr(unique(as.character(dd$site)), function(s) {
    ddrop <- dd %>% filter(site != s)
    est <- tryCatch(fit_one(formula, term, ddrop),
                    error = function(e) tibble(estimate = NA_real_,
                                               p.value = NA_real_))
    tibble(effect = effect, dropped_site = s, estimate = est$estimate[1],
           p.value = est$p.value[1], n_sites = n_distinct(ddrop$site),
           n = nrow(ddrop))
  })
})
write_csv(loso, file.path(OUT, "opportunity_loso.csv"))

loso_summary <- loso %>% group_by(effect) %>%
  summarise(min_est = min(estimate, na.rm = TRUE),
            max_est = max(estimate, na.rm = TRUE),
            sign_flips = sum(sign(estimate) != sign(median(estimate, na.rm = TRUE)),
                             na.rm = TRUE),
            max_p = max(p.value, na.rm = TRUE),
            .groups = "drop")
write_csv(loso_summary, file.path(OUT, "opportunity_loso_summary.csv"))

boot <- pmap_dfr(EFFECTS, function(effect, outcome, formula, fixed_formula,
                                   term, opp_flag) {
  dd <- d %>% filter(.data[[opp_flag]])
  sites <- unique(as.character(dd$site))
  ests <- map_dbl(seq_len(1000), function(b) {
    samp <- sample(sites, length(sites), replace = TRUE)
    db <- map_dfr(seq_along(samp), function(i)
      dd %>% filter(site == samp[i]) %>% mutate(site_b = paste0("s", i)))
    f <- as.formula(str_replace(fixed_formula, "\\+ site$", "+ site_b"))
    fit <- tryCatch(lm(f, data = db), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    tt <- broom::tidy(fit) %>% filter(.data$term == !!term)
    if (nrow(tt)) tt$estimate[1] else NA_real_
  })
  es <- ests[!is.na(ests)]
  tibble(effect = effect, n_boot = length(es),
         ci_lo = quantile(es, 0.025), ci_hi = quantile(es, 0.975),
         excludes_zero = sign(quantile(es, 0.025)) == sign(quantile(es, 0.975)))
})
write_csv(boot, file.path(OUT, "opportunity_cluster_bootstrap.csv"))

print(site_opp %>% select(site_label, n, starts_with("opp_")), n = 20)
cat("\n--- Opportunity models ---\n")
print(models %>% select(effect, sample, n_sites, n, estimate, p.value), n = 40)
cat("\n--- LOSO summary ---\n")
print(loso_summary, n = 20)
cat("\n--- Cluster bootstrap ---\n")
print(boot, n = 20)
