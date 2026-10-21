# Analysis 22 - Weather and haze as quality proxies for the sky (and as
# confound checks for the season findings). Date covariates fetched by
# fetch_date_covariates.py (meteostat Shenyang + quotsoft AQI/PM2.5; 32 dates).
# Families: (a) date-level main effects on evaluation; (b) sky x weather/haze
# quality interactions; (c) do the A17/A20 season interactions survive
# temperature + AQI adjustment (thermal-comfort / haze confound check)?

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_22_weather_aqi"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

wx <- read_csv(file.path(ANALYSIS_DIR, "data/date_covariates.csv"),
               show_col_types = FALSE) %>%
  mutate(rain_day = as.integer(replace_na(prcp, 0) > 0))

d <- load_vis_data() %>%
  mutate(date = parse_survey_date(survey_date),
         season = factor(sub(".* ", "", wave_label),
                         levels = c("autumn", "spring", "summer"))) %>%
  left_join(wx, by = "date") %>%
  mutate(temp_z = as.numeric(scale(temp)),
         rhum_z = as.numeric(scale(rhum)),
         aqi_z = as.numeric(scale(aqi_mean)),
         wspd_z = as.numeric(scale(wspd)))
cat("Photo-valid rows with date covariates:", sum(!is.na(d$aqi_mean)),
    "of", nrow(d), "| distinct dates:", n_distinct(na.omit(d$date)), "\n")
dd <- filter(d, !is.na(aqi_mean))

OUTCOMES <- c("isop_z", "overall_soundscape_z", "harmony_z", "isoe_z")
WXVARS <- c("temp_z", "rhum_z", "aqi_z", "wspd_z", "rain_day")

# (a) Date-level main effects (FDR = 20); date RE carries the day grain -------
mains <- expand_grid(outcome = OUTCOMES, wxvar = WXVARS) %>% rowwise() %>%
  mutate(res = list({
    f <- as.formula(paste0(outcome, " ~ ", wxvar,
                           " + (1 | site) + (1 | date)"))
    fit <- suppressMessages(lmer(f, data = dd))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>% filter(term == wxvar)
    tibble(estimate = tt$estimate, std.error = tt$std.error,
           p.value = tt$p.value, singular = isSingular(fit))
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q.value = fdr(p.value)) %>% arrange(q.value)
write_csv(mains, file.path(OUT, "date_covariate_mains.csv"))

# (b) Sky quality: sky share x weather/haze on ISO-P and OA (FDR = 8) --------
sky <- expand_grid(outcome = c("isop_z", "overall_soundscape_z"),
                   wxvar = c("temp_z", "rhum_z", "aqi_z", "rain_day")) %>%
  rowwise() %>%
  mutate(res = list({
    f <- as.formula(paste0(outcome, " ~ comp_sky_z * ", wxvar,
                           " + (1 | site) + (1 | date)"))
    fit <- suppressMessages(lmer(f, data = dd))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(str_detect(term, ":"))
    tibble(b_interaction = tt$estimate, p.value = tt$p.value,
           singular = isSingular(fit))
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q.value = fdr(p.value)) %>% arrange(q.value)
write_csv(sky, file.path(OUT, "sky_quality_interactions.csv"))

# cloud-cover variant on the 13 covered dates (exploratory, no FDR family)
dc <- filter(dd, !is.na(cldc)) %>% mutate(cldc_z = as.numeric(scale(cldc)))
cloud <- map_dfr(c("isop_z", "overall_soundscape_z"), function(y) {
  fit <- suppressMessages(lmer(as.formula(paste0(
    y, " ~ comp_sky_z * cldc_z + (1 | site) + (1 | date)")), data = dc))
  broom.mixed::tidy(fit, effects = "fixed") %>%
    filter(str_detect(term, ":")) %>%
    transmute(outcome = y, n = nrow(dc), term, estimate, p.value)
})
write_csv(cloud, file.path(OUT, "sky_cloud_exploratory.csv"))

# (c) Season interactions under weather adjustment (multi-season sites) ------
tab <- table(as.character(dd$site_label), dd$season)
multi <- rownames(tab)[rowSums(tab > 0) >= 2]
dm <- dd %>% filter(site_label %in% multi) %>% droplevels()
cat("Multi-season rows with covariates:", nrow(dm), "\n")

season_checks <- tribble(
  ~label, ~outcome, ~term,
  "greenery_x_season_oa", "overall_soundscape_z", "comp_greenery_z",
  "natural_x_season_oa", "overall_soundscape_z", "natural_sounds_z",
  "natural_x_season_isop", "isop_z", "natural_sounds_z",
  "natural_x_season_harmony", "harmony_z", "natural_sounds_z") %>%
  rowwise() %>%
  mutate(res = list({
    base <- " + (1 | site) + (1 | site:wave)"
    f1 <- paste0(outcome, " ~ ", term, " * season + temp_z + aqi_z", base)
    f0 <- paste0(outcome, " ~ ", term, " + season + temp_z + aqi_z", base)
    m1 <- suppressMessages(lmer(as.formula(f1), data = dm, REML = FALSE))
    m0 <- suppressMessages(lmer(as.formula(f0), data = dm, REML = FALSE))
    cc <- broom.mixed::tidy(m1, effects = "fixed") %>%
      filter(str_detect(term, ":seasonspring"))
    tibble(b_x_spring = cc$estimate, p_x_spring = cc$p.value,
           p_joint = anova(m0, m1)$`Pr(>Chisq)`[2])
  })) %>% ungroup() %>% unnest(res)
write_csv(season_checks, file.path(OUT, "season_weather_adjusted.csv"))

print(mains, n = 20); print(sky, n = 8); print(cloud); print(season_checks)
