# Analysis 19 - Loudness-conditioned visual effects (attention hypothesis):
# does vision matter only where sound is not capturing attention? Reverse of
# A07 (which tested vision moderating sound slopes): here LAeq moderates the
# VISUAL slopes. Continuous visual x LAeq interactions (FDR family) plus
# LAeq-tertile-stratified slopes; perceived loudness as sensitivity moderator.
# Sample: photo-valid rows with acoustics under the available-channel rule
# (2026-08-07: mono rows enter with their single valid channel; the
# channel-stratified sensitivity at the end separates the binaural and
# single-channel strata, which are confounded with survey wave).

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_19_loudness_conditioned_vision"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  filter(!is.na(laeq_a)) %>%
  mutate(comp_traffic_vis = comp_vehicles + ph_road,
         comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
         water_any = as.integer(comp_water > 0.005),
         flower_any = as.integer(replace_na(ph_flower, 0) > 0.001),
         laeq_tertile = cut(laeq_a, quantile(laeq_a, c(0, 1/3, 2/3, 1)),
                            labels = c("quiet", "mid", "loud"),
                            include.lowest = TRUE))
cat("Rows with acoustics:", nrow(d), "| tertile cuts:",
    paste(round(quantile(d$laeq_a, c(1/3, 2/3)), 1), collapse = "/"), "dB\n")

EFFECTS <- tribble(
  ~effect, ~outcome, ~term,
  "greenery_isop", "isop_z", "comp_greenery_z",
  "greenery_oa", "overall_soundscape_z", "comp_greenery_z",
  "flower_oa", "overall_soundscape_z", "flower_any",
  "water_oa", "overall_soundscape_z", "water_any",
  "vehicles_isoe", "isoe_z", "comp_vehicles_z",
  "sky_isop", "isop_z", "comp_sky_z",
  "people_oa", "overall_soundscape_z", "comp_people_z")

# Continuous moderation: visual x LAeq (FDR = 7) ------------------------------
moderation <- EFFECTS %>% rowwise() %>%
  mutate(res = list({
    f <- as.formula(paste0(outcome, " ~ ", term, " * laeq_a_z",
                           " + (1 | site) + (1 | site:wave)"))
    fit <- suppressMessages(lmer(f, data = d))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(str_detect(term, ":"))
    tibble(b_interaction = tt$estimate, se = tt$std.error, p.value = tt$p.value,
           singular = isSingular(fit))
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q.value = fdr(p.value)) %>% arrange(q.value)
write_csv(moderation, file.path(OUT, "laeq_moderation.csv"))

# sensitivity: perceived loudness as the moderator
moderation_pl <- EFFECTS %>% rowwise() %>%
  mutate(res = list({
    f <- as.formula(paste0(outcome, " ~ ", term, " * perceived_loudness_z",
                           " + (1 | site) + (1 | site:wave)"))
    fit <- suppressMessages(lmer(f, data = d))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(str_detect(term, ":"))
    tibble(b_interaction = tt$estimate, se = tt$std.error, p.value = tt$p.value)
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q.value = fdr(p.value)) %>% arrange(q.value)
write_csv(moderation_pl, file.path(OUT, "perceived_loudness_moderation.csv"))

# Tertile-stratified slopes ----------------------------------------------------
stratified <- map_dfr(levels(d$laeq_tertile), function(s) {
  ds <- filter(d, laeq_tertile == s)
  EFFECTS %>% rowwise() %>%
    mutate(res = list({
      tn <- term
      f <- as.formula(paste0(outcome, " ~ ", tn,
                             " + (1 | site) + (1 | site:wave)"))
      fit <- tryCatch(suppressMessages(lmer(f, data = ds)),
                      error = function(e) NULL)
      if (is.null(fit)) tibble(estimate = NA_real_, std.error = NA_real_,
                               p.value = NA_real_, singular = NA)
      else {
        tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
          filter(term == tn)
        tibble(estimate = tt$estimate, std.error = tt$std.error,
               p.value = tt$p.value, singular = isSingular(fit))
      }
    })) %>% ungroup() %>% unnest(res) %>%
    transmute(laeq_tertile = s, n = nrow(ds), effect, estimate, std.error,
              p.value, singular)
})
write_csv(stratified, file.path(OUT, "tertile_stratified_slopes.csv"))

# Channel-stratified sensitivity (available-channel rule, 2026-08-07): the
# LAeq-moderation family per acoustic_channels stratum. The strata are
# confounded with survey wave (waves 2/7 largely or wholly mono), so
# differences read as channel-cum-era heterogeneity, not channel noise alone.
channel_sens <- map_dfr(c(2, 1), function(ch) {
  dd <- filter(d, acoustic_channels == ch)
  EFFECTS %>% rowwise() %>%
    mutate(res = list({
      f <- as.formula(paste0(outcome, " ~ ", term, " * laeq_a_z",
                             " + (1 | site) + (1 | site:wave)"))
      fit <- suppressMessages(lmer(f, data = dd))
      tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
        filter(str_detect(term, ":"))
      tibble(b_interaction = tt$estimate, se = tt$std.error,
             p.value = tt$p.value, singular = isSingular(fit))
    })) %>% ungroup() %>% unnest(res) %>%
    mutate(q.value = fdr(p.value), channels = ch, n = nrow(dd)) %>%
    arrange(q.value)
})
write_csv(channel_sens, file.path(OUT, "channel_sensitivity.csv"))

print(moderation); print(moderation_pl)
print(stratified %>% arrange(effect, laeq_tertile), n = 30)
print(channel_sens %>% select(channels, n, effect, b_interaction, p.value, q.value), n = 14)
