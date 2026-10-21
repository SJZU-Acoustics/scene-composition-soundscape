# Analysis 20 - Auditory phenology: sound sources have seasonal aspect too
# (birdsong in spring, insects in summer, near-silence in autumn). Do the
# perceived-source slopes on evaluation vary by season, mirroring A17's visual
# test? Multi-season sites only; source x season interactions (FDR family) and
# season-stratified slopes. Sources: natural, traffic, human (4-cat battery).

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_20_auditory_phenology"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(season = factor(sub(".* ", "", wave_label),
                         levels = c("autumn", "spring", "summer")))
tab <- table(as.character(d$site_label), d$season)
multi <- rownames(tab)[rowSums(tab > 0) >= 2]
dm <- d %>% filter(site_label %in% multi) %>% droplevels()
cat("Multi-season n =", nrow(dm), "\n")

# mean natural-sound perception by season (auditory quantity-shift check)
season_means <- dm %>% group_by(season) %>%
  summarise(n = n(), across(all_of(SOURCES), \(x) mean(x, na.rm = TRUE)),
            .groups = "drop")
write_csv(season_means, file.path(OUT, "source_means_by_season.csv"))

PREDICTORS <- c("natural_sounds_z", "traffic_noise_z", "human_sounds_z")
OUTCOMES <- c("isop_z", "overall_soundscape_z", "harmony_z", "isoe_z")

# Source x season interactions (FDR = 12) -------------------------------------
interactions <- expand_grid(outcome = OUTCOMES, predictor = PREDICTORS) %>%
  rowwise() %>%
  mutate(res = list({
    f1 <- as.formula(paste0(outcome, " ~ ", predictor, " * season",
                            " + (1 | site) + (1 | site:wave)"))
    f0 <- as.formula(paste0(outcome, " ~ ", predictor, " + season",
                            " + (1 | site) + (1 | site:wave)"))
    m1 <- suppressMessages(lmer(f1, data = dm, REML = FALSE))
    m0 <- suppressMessages(lmer(f0, data = dm, REML = FALSE))
    cc <- broom.mixed::tidy(m1, effects = "fixed") %>%
      filter(str_detect(term, paste0(predictor, ":")))
    tibble(b_x_spring = cc$estimate[str_detect(cc$term, "spring")],
           b_x_summer = cc$estimate[str_detect(cc$term, "summer")],
           p_joint = anova(m0, m1)$`Pr(>Chisq)`[2])
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q_joint = fdr(p_joint)) %>% arrange(q_joint)
write_csv(interactions, file.path(OUT, "source_season_interactions.csv"))

# Season-stratified source slopes ----------------------------------------------
stratified <- map_dfr(levels(dm$season), function(s) {
  ds <- filter(dm, season == s)
  expand_grid(outcome = OUTCOMES, predictor = PREDICTORS) %>%
    rowwise() %>%
    mutate(res = list({
      f <- as.formula(paste0(outcome, " ~ ", predictor,
                             " + (1 | site) + (1 | site:wave)"))
      fit <- tryCatch(suppressMessages(lmer(f, data = ds)),
                      error = function(e) NULL)
      if (is.null(fit)) tibble(estimate = NA_real_, std.error = NA_real_,
                               p.value = NA_real_, singular = NA)
      else {
        tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
          filter(term == predictor)
        tibble(estimate = tt$estimate, std.error = tt$std.error,
               p.value = tt$p.value, singular = isSingular(fit))
      }
    })) %>% ungroup() %>% unnest(res) %>%
    mutate(season = s, n = nrow(ds), .before = 1)
})
write_csv(stratified, file.path(OUT, "stratified_source_slopes.csv"))

# Robustness of the natural x season survivors: same-year 2019, LOSO ---------
d19 <- dm %>% filter(str_detect(wave_label, "^2019")) %>% droplevels()
same_year <- map_dfr(c("overall_soundscape_z", "isop_z", "harmony_z"),
                     function(y) {
  f1 <- as.formula(paste0(y, " ~ natural_sounds_z * season + (1 | site) + (1 | site:wave)"))
  f0 <- as.formula(paste0(y, " ~ natural_sounds_z + season + (1 | site) + (1 | site:wave)"))
  m1 <- suppressMessages(lmer(f1, data = d19, REML = FALSE))
  m0 <- suppressMessages(lmer(f0, data = d19, REML = FALSE))
  cc <- broom.mixed::tidy(m1, effects = "fixed") %>%
    filter(str_detect(term, "natural_sounds_z:"))
  tibble(outcome = y, n = nrow(d19),
         b_x_spring = cc$estimate[str_detect(cc$term, "spring")],
         p_x_spring = cc$p.value[str_detect(cc$term, "spring")],
         p_joint = anova(m0, m1)$`Pr(>Chisq)`[2])
})
write_csv(same_year, file.path(OUT, "same_year_2019_check.csv"))

loso <- map_dfr(multi, function(s) {
  dl <- dm %>% filter(site_label != s) %>% droplevels()
  f1 <- as.formula("overall_soundscape_z ~ natural_sounds_z * season + (1 | site) + (1 | site:wave)")
  f0 <- as.formula("overall_soundscape_z ~ natural_sounds_z + season + (1 | site) + (1 | site:wave)")
  m1 <- suppressMessages(lmer(f1, data = dl, REML = FALSE))
  m0 <- suppressMessages(lmer(f0, data = dl, REML = FALSE))
  cc <- broom.mixed::tidy(m1, effects = "fixed") %>%
    filter(str_detect(term, "natural_sounds_z:") & str_detect(term, "spring"))
  tibble(dropped_site = s, b_x_spring = cc$estimate, p_x_spring = cc$p.value,
         p_joint = anova(m0, m1)$`Pr(>Chisq)`[2])
})
write_csv(loso, file.path(OUT, "loso_natural_oa.csv"))

print(season_means); print(interactions, n = 12)
print(stratified %>% arrange(predictor, outcome, season), n = 40)
print(same_year); print(loso)
