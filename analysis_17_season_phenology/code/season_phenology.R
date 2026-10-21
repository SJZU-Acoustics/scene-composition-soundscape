# Analysis 17 - Season / phenology (Prof Zhang's hint, 2026-07-09): the photo
# variables measure visual QUANTITY (area shares), not quality; seasonal aspect
# (ji-xiang) proxies the missing quality axis. Questions:
#   (a) does greenery area share itself vary by season? (quantity-constancy)
#   (b) greenery x season interactions on the evaluation battery
#   (c) season-stratified slopes for greenery and the A12 key effects
#   (d) confound checks: same-year (2019) contrast, LOSO, cluster bootstrap
#       (1,000, site:wave batches), month-level exploratory phenology
# Identification: only 4 sites span >= 2 seasons (Zhongshan Park, Zhongshan
# Square, Dadong Square, Olympic Park); all season models run on that subset.

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_17_season_phenology"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(
    season = factor(sub(".* ", "", wave_label),
                    levels = c("autumn", "spring", "summer")),
    year = as.integer(sub(" .*", "", wave_label)),
    month = as.integer(format(parse_survey_date(survey_date), "%m")),
    # date-derived season (sensitivity). Under the settled 0717 workbook the
    # wave codes are date-consistent (the 2019 batch conflicts were fixed in
    # the 0716 round), so season and season_adj should now agree.
    season_adj = factor(case_when(
      !is.na(month) & month %in% c(3, 4, 5) ~ "spring",
      !is.na(month) & month %in% c(6, 7, 8) ~ "summer",
      !is.na(month) & month %in% c(9, 10, 11) ~ "autumn",
      TRUE ~ as.character(season)), levels = levels(season)),
    comp_traffic_vis = comp_vehicles + ph_road,
    comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
    water_any = as.integer(comp_water > 0.005),
    flower_any = as.integer(replace_na(ph_flower, 0) > 0.001))

season_tab <- table(as.character(d$site_label), d$season)
multi_sites <- rownames(season_tab)[rowSums(season_tab > 0) >= 2]
dm <- d %>% filter(site_label %in% multi_sites) %>% droplevels()
cat("Multi-season sites:", paste(multi_sites, collapse = ", "),
    "| n =", nrow(dm), "\n")

EVAL_OUTCOMES <- c(isop_z = "ISO pleasantness", isoe_z = "ISO eventfulness",
                   overall_soundscape_z = "Overall evaluation",
                   harmony_z = "Perceived harmony",
                   perceived_loudness_z = "Perceived loudness")

# (a) Quantity constancy: composites and element prevalence by season --------
quantity <- dm %>% group_by(season) %>%
  summarise(n = n(), across(all_of(COMPS8), mean),
            flower_any = mean(flower_any), water_any = mean(water_any),
            .groups = "drop")
write_csv(quantity, file.path(OUT, "quantity_by_season.csv"))

# (b) Greenery x season interactions on the evaluation battery ---------------
int_fit <- function(y, dd, season_var = "season") {
  f1 <- as.formula(paste0(y, " ~ comp_greenery_z * ", season_var,
                          " + (1 | site) + (1 | site:wave)"))
  f0 <- as.formula(paste0(y, " ~ comp_greenery_z + ", season_var,
                          " + (1 | site) + (1 | site:wave)"))
  m1 <- lmer(f1, data = dd, REML = FALSE)
  m0 <- lmer(f0, data = dd, REML = FALSE)
  cc <- broom.mixed::tidy(m1, effects = "fixed") %>%
    filter(str_detect(term, "comp_greenery_z"))
  tibble(outcome = y,
         b_autumn = cc$estimate[cc$term == "comp_greenery_z"],
         b_x_spring = cc$estimate[str_detect(cc$term, "spring")],
         p_x_spring = cc$p.value[str_detect(cc$term, "spring")],
         b_x_summer = cc$estimate[str_detect(cc$term, "summer")],
         p_x_summer = cc$p.value[str_detect(cc$term, "summer")],
         p_joint = anova(m0, m1)$`Pr(>Chisq)`[2])
}
interactions <- map_dfr(names(EVAL_OUTCOMES), int_fit, dd = dm) %>%
  mutate(q_joint = fdr(p_joint))
write_csv(interactions, file.path(OUT, "greenery_season_interactions.csv"))

# sensitivity: date-corrected season labels
interactions_adj <- map_dfr(names(EVAL_OUTCOMES), int_fit, dd = dm,
                            season_var = "season_adj") %>%
  mutate(q_joint = fdr(p_joint))
write_csv(interactions_adj,
          file.path(OUT, "greenery_season_interactions_dateadj.csv"))

# (c) Season-stratified slopes: greenery + the A12 key effects ---------------
strat_fit <- function(f, dd, term) {
  fit <- tryCatch(suppressMessages(lmer(as.formula(f), data = dd)),
                  error = function(e) NULL)
  if (is.null(fit)) return(tibble(estimate = NA_real_, std.error = NA_real_,
                                  p.value = NA_real_))
  tt <- broom.mixed::tidy(fit, effects = "fixed") %>% filter(term == !!term)
  tibble(estimate = tt$estimate[1], std.error = tt$std.error[1],
         p.value = tt$p.value[1],
         singular = isSingular(fit))
}

STRAT_EFFECTS <- tribble(
  ~effect, ~f, ~term,
  "greenery_isop", "isop_z ~ comp_greenery_z + (1 | site) + (1 | site:wave)", "comp_greenery_z",
  "greenery_oa", "overall_soundscape_z ~ comp_greenery_z + (1 | site) + (1 | site:wave)", "comp_greenery_z",
  "greenery_harmony", "harmony_z ~ comp_greenery_z + (1 | site) + (1 | site:wave)", "comp_greenery_z",
  "vehicles_isoe", "isoe_z ~ comp_vehicles_z + (1 | site) + (1 | site:wave)", "comp_vehicles_z",
  "traffic_congruence_oa", "overall_soundscape_z ~ traffic_noise_z * comp_traffic_vis_z + (1 | site) + (1 | site:wave)", "traffic_noise_z:comp_traffic_vis_z",
  "water_oa", "overall_soundscape_z ~ water_any + (1 | site) + (1 | site:wave)", "water_any",
  "flower_oa", "overall_soundscape_z ~ flower_any + (1 | site) + (1 | site:wave)", "flower_any")

stratified <- map_dfr(levels(dm$season), function(s) {
  ds <- filter(dm, season == s)
  STRAT_EFFECTS %>% rowwise() %>%
    mutate(res = list(strat_fit(f, ds, term))) %>% ungroup() %>%
    unnest(res) %>%
    transmute(season = s, n = nrow(ds), effect, estimate, std.error,
              p.value, singular)
})
write_csv(stratified, file.path(OUT, "stratified_slopes.csv"))

# (d1) Same-year check: 2019 spans spring/summer/autumn (waves 2/3/4) --------
d19 <- dm %>% filter(year == 2019) %>% droplevels()
same_year <- map_dfr(c("overall_soundscape_z", "isop_z"), int_fit, dd = d19,
                     season_var = "season_adj")
write_csv(same_year %>% mutate(n = nrow(d19)),
          file.path(OUT, "same_year_2019_check.csv"))

# (d2) LOSO over the 4 multi-season sites: spring greenery slopes ------------
dspring <- filter(dm, season == "spring")
loso <- map_dfr(multi_sites, function(s) {
  map_dfr(c(greenery_oa = "overall_soundscape_z", greenery_isop = "isop_z"),
          function(y) {
    strat_fit(paste0(y, " ~ comp_greenery_z + (1 | site) + (1 | site:wave)"),
              filter(dspring, site_label != s), "comp_greenery_z") %>%
      mutate(outcome = y)
  }) %>% mutate(dropped_site = s, .before = 1)
})
write_csv(loso, file.path(OUT, "loso_spring_slopes.csv"))

# (d3) Cluster bootstrap (1,000) resampling site:wave batches ----------------
# Only 4 sites span seasons, so site-level resampling is too coarse; the
# season contrast lives at the site:wave batch level (19 batches).
dm <- dm %>% mutate(batch = interaction(site, wave, drop = TRUE))
batches <- levels(dm$batch)
boot_one <- function(y) {
  ests <- map_dbl(seq_len(1000), function(b) {
    samp <- sample(batches, length(batches), replace = TRUE)
    db <- map_dfr(seq_along(samp), function(i)
      dm %>% filter(batch == samp[i]) %>% mutate(batch_b = paste0("b", i)))
    if (n_distinct(db$season) < 2) return(NA_real_)
    fit <- tryCatch(suppressMessages(lmer(
      as.formula(paste0(y, " ~ comp_greenery_z * season + (1 | site) + (1 | batch_b)")),
      data = db, REML = FALSE)), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(str_detect(term, "greenery") & str_detect(term, "spring"))
    if (nrow(tt)) tt$estimate[1] else NA_real_
  })
  es <- ests[!is.na(ests)]
  tibble(outcome = y, n_boot = length(es),
         ci_lo = quantile(es, 0.025), ci_hi = quantile(es, 0.975),
         excludes_zero = sign(quantile(es, 0.025)) == sign(quantile(es, 0.975)))
}
boot <- map_dfr(c("overall_soundscape_z", "isop_z"), boot_one)
write_csv(boot, file.path(OUT, "cluster_bootstrap_greenery_x_spring.csv"))

# (d4) Month-level exploratory phenology (rows with survey dates) ------------
by_month <- dm %>% filter(!is.na(month)) %>% group_by(month) %>%
  filter(n() >= 60) %>% group_modify(function(g, k) {
    bind_rows(
      strat_fit("overall_soundscape_z ~ comp_greenery_z + (1 | site) + (1 | site:wave)",
                g, "comp_greenery_z") %>% mutate(outcome = "overall_soundscape_z"),
      strat_fit("isop_z ~ comp_greenery_z + (1 | site) + (1 | site:wave)",
                g, "comp_greenery_z") %>% mutate(outcome = "isop_z")) %>%
      mutate(n = nrow(g))
  }) %>% ungroup()
write_csv(by_month, file.path(OUT, "month_exploratory_slopes.csv"))

print(quantity, width = Inf)
print(interactions); print(interactions_adj)
print(stratified %>% arrange(effect, season), n = 40)
print(same_year); print(loso, n = 20); print(boot)
print(by_month)
