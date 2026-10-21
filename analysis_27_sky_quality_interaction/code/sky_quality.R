# A27 - sky quality x sky share: the raw-photo route A22 asked for.
# A22's exploratory finding (13 cloud-covered dates): sky share helps ISO-P only
# under clear skies (comp_sky x cloud cover -0.165, p = 2e-4, wave-confounded).
# Here the moderator is the per-respondent photo-measured sky blueness
# (sky_blue_ratio = blue / (blue + white) among top-region sky pixels).
# Declared family: 2 interactions (ISO-P, OA) x photo moderator = 2 tests, BH.
# Meteorological validation: date-mean sky_blue_ratio vs cloud cover (cldc).
# Pre-scoping exploratory layer (raw-photo delivery 2026-08-18).
suppressPackageStartupMessages({library(tidyverse); library(lmerTest)})
source("analysis_23_photo_quality_metrics/code/photo_helpers.R")

ADIR <- "analysis_27_sky_quality_interaction"
OUT <- file.path(ADIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20260818)

wx <- read_csv("analysis_22_weather_aqi/data/date_covariates.csv",
               show_col_types = FALSE) %>% mutate(date = as.Date(date))

d <- load_quality_data() %>%
  mutate(equip = factor(equip_proxy),
         date = parse_survey_date(survey_date),
         sky_blue_ratio_z = as.numeric(scale(sky_blue_ratio))) %>%
  left_join(wx, by = "date")

# --- 1. meteorological validation: date-mean blueness vs cloud cover ---
dv <- d %>% filter(!is.na(cldc), !is.na(sky_blue_ratio)) %>%
  group_by(date) %>%
  summarise(blue = mean(sky_blue_ratio, na.rm = TRUE), cldc = first(cldc),
            n = n(), .groups = "drop")
r_val <- cor.test(dv$blue, dv$cldc, method = "spearman", exact = FALSE)
cat("Date-level validation: mean sky_blue_ratio vs cloud cover, n dates =", nrow(dv), "\n")
cat(sprintf("Spearman rho = %.3f, p = %.4f\n", r_val$estimate, r_val$p.value))
write_csv(dv, file.path(OUT, "date_blueness_vs_cloud.csv"))

# --- 2. declared interactions: comp_sky x sky_blue_ratio on ISO-P and OA ---
inter <- map_dfr(c("isop_z", "overall_soundscape_z"), function(oc) {
  dd <- d %>% drop_na(sky_blue_ratio, comp_sky)
  fit <- lmer(reformulate(c("comp_sky_z * sky_blue_ratio_z", "wave", "equip",
                            "(1 | site)", "(1 | panorama_group)"), oc), data = dd)
  co <- summary(fit)$coefficients
  tibble(outcome = oc, n = nrow(dd),
         b_sky = co["comp_sky_z", "Estimate"], p_sky = co["comp_sky_z", "Pr(>|t|)"],
         b_blue = co["sky_blue_ratio_z", "Estimate"], p_blue = co["sky_blue_ratio_z", "Pr(>|t|)"],
         b_inter = co["comp_sky_z:sky_blue_ratio_z", "Estimate"],
         se_inter = co["comp_sky_z:sky_blue_ratio_z", "Std. Error"],
         p_inter = co["comp_sky_z:sky_blue_ratio_z", "Pr(>|t|)"])
}) %>% mutate(q_inter = p.adjust(p_inter, "BH"))
write_csv(inter, file.path(OUT, "sky_share_x_blueness.csv"))
cat("\nDeclared interactions (family n = 2):\n"); print(inter, width = Inf)

# --- 3. stratified reading: sky-share slope by blueness tertile (descriptive) ---
strat <- d %>% drop_na(sky_blue_ratio) %>%
  mutate(blue_ter = cut(sky_blue_ratio, quantile(sky_blue_ratio, c(0, 1/3, 2/3, 1)),
                        labels = c("white/overcast", "mixed", "blue"), include.lowest = TRUE))
strat_res <- map_dfr(levels(strat$blue_ter), function(t) {
  dd <- strat %>% filter(blue_ter == t)
  fit <- lmer(isop_z ~ comp_sky_z + wave + equip + (1 | site) + (1 | panorama_group), data = dd)
  co <- summary(fit)$coefficients
  tibble(tertile = t, n = nrow(dd), b_sky = co["comp_sky_z", "Estimate"],
         se = co["comp_sky_z", "Std. Error"], p = co["comp_sky_z", "Pr(>|t|)"])
})
write_csv(strat_res, file.path(OUT, "sky_slope_by_blueness_tertile.csv"))
cat("\nISO-P sky-share slope by blueness tertile:\n"); print(strat_res)

# --- 4. cloud-cover replication on the joined sample (exploratory) ---
dc <- d %>% filter(!is.na(cldc)) %>% droplevels() %>%
  mutate(cldc_z = as.numeric(scale(cldc)))
cloud <- map_dfr(c("isop_z", "overall_soundscape_z"), function(oc) {
  rhs <- c("comp_sky_z * cldc_z",
           if (nlevels(dc$equip) > 1) "equip", "(1 | site)", "(1 | date)")
  fit <- lmer(reformulate(rhs, oc), data = dc)
  co <- summary(fit)$coefficients
  tibble(outcome = oc, n = nrow(dc), b_inter = co["comp_sky_z:cldc_z", "Estimate"],
         p_inter = co["comp_sky_z:cldc_z", "Pr(>|t|)"])
})
write_csv(cloud, file.path(OUT, "cloud_replication.csv"))
cat("\nCloud-cover variant (13-date subsample, exploratory):\n"); print(cloud)

# --- 5. sky blueness main-effect robustness recap on ISO-P family (context for A26) ---
fitb <- lmer(overall_soundscape_z ~ sky_blue_ratio_z + comp_sky_z + wave + equip +
               (1 | site) + (1 | panorama_group), data = d %>% drop_na(sky_blue_ratio))
cb <- summary(fitb)$coefficients
write_csv(as_tibble(cb, rownames = "term"), file.path(OUT, "blueness_sky_adjusted.csv"))
cat("\nOA ~ blueness (sky-share adjusted):\n")
print(cb[c("sky_blue_ratio_z", "comp_sky_z"), , drop = FALSE])

# --- 6. blueness main effect adjusted for daily weather (A26 robustness) ---
# Temperature and AQI from A22's date covariates, standardised over the joined
# sample; a survey-date random intercept absorbs the remaining day-level variation.
dw <- d %>% drop_na(sky_blue_ratio, temp, aqi_mean) %>%
  mutate(temp_z = as.numeric(scale(temp)), aqi_z = as.numeric(scale(aqi_mean)))
fitw <- lmer(overall_soundscape_z ~ sky_blue_ratio_z + comp_sky_z + temp_z + aqi_z + wave + equip +
               (1 | site) + (1 | date) + (1 | panorama_group), data = dw)
cw <- summary(fitw)$coefficients
write_csv(as_tibble(cw, rownames = "term"), file.path(OUT, "blueness_weather_adjusted.csv"))
cat("\nOA ~ blueness, weather-adjusted (n =", nrow(dw), "):\n")
print(cw["sky_blue_ratio_z", , drop = FALSE])
