# Analysis 21 - Familiarity subgroups: do visual effects depend on who is
# looking? Hypothesis: scene quality matters more for less-familiar visitors
# (tourists) than for locals who carry a stored image of the place.
# Visitor coding: round-3 harmonised visitor_group (local / tourist / other;
# built by build/build_core.py - old-era 1/2/3 = city resident/tourist/other,
# new-era 2023+ 1+2/3/4). "other" is semantically mixed and excluded.
# Activity-state subgroups use the 0716 codebook (状态 1 staying / 2 arriving /
# 3 leaving; waves 9/10 asked an open activity question instead -> missing).

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_21_familiarity_subgroups"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d0 <- load_vis_data() %>%
  mutate(comp_traffic_vis = comp_vehicles + ph_road,
         comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
         water_any = as.integer(comp_water > 0.005),
         flower_any = as.integer(replace_na(ph_flower, 0) > 0.001))

d <- d0 %>%
  mutate(visitor_status = fct_drop(visitor_group),
         non_local = case_when(
           visitor_group == "local" ~ 0L,
           visitor_group == "tourist" ~ 1L,
           TRUE ~ NA_integer_)) %>%
  filter(!is.na(non_local)) %>%
  mutate(visitor_status = fct_drop(visitor_status))
cat("n =", nrow(d), "| local:", sum(d$non_local == 0),
    "| non-local (tourist):", sum(d$non_local == 1), "\n")

EFFECTS <- tribble(
  ~effect, ~outcome, ~term,
  "greenery_isop", "isop_z", "comp_greenery_z",
  "greenery_oa", "overall_soundscape_z", "comp_greenery_z",
  "flower_oa", "overall_soundscape_z", "flower_any",
  "water_oa", "overall_soundscape_z", "water_any",
  "vehicles_isoe", "isoe_z", "comp_vehicles_z",
  "sky_isop", "isop_z", "comp_sky_z",
  "people_oa", "overall_soundscape_z", "comp_people_z")

# Visual x non-local interactions (FDR = 7) -----------------------------------
moderation <- EFFECTS %>% rowwise() %>%
  mutate(res = list({
    tn <- term
    f <- as.formula(paste0(outcome, " ~ ", tn, " * non_local + wave",
                           " + (1 | site)"))
    fit <- suppressMessages(lmer(f, data = d))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(term == paste0(tn, ":non_local"))
    tibble(b_interaction = tt$estimate, se = tt$std.error,
           p.value = tt$p.value, singular = isSingular(fit))
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q.value = fdr(p.value)) %>% arrange(q.value)
write_csv(moderation, file.path(OUT, "non_local_moderation.csv"))

# Stratified slopes: resident / city_other / tourist ---------------------------
stratified <- map_dfr(levels(d$visitor_status), function(g) {
  ds <- filter(d, visitor_status == g)
  EFFECTS %>% rowwise() %>%
    mutate(res = list({
      tn <- term
      f <- as.formula(paste0(outcome, " ~ ", tn, " + wave + (1 | site)"))
      fit <- tryCatch(suppressMessages(lmer(f, data = ds)),
                      error = function(e) NULL)
      if (is.null(fit)) tibble(estimate = NA_real_, std.error = NA_real_,
                               p.value = NA_real_, singular = NA)
      else {
        tt <- broom.mixed::tidy(fit, effects = "fixed") %>% filter(term == tn)
        tibble(estimate = tt$estimate, std.error = tt$std.error,
               p.value = tt$p.value, singular = isSingular(fit))
      }
    })) %>% ungroup() %>% unnest(res) %>%
    transmute(group = g, n = nrow(ds), effect, estimate, std.error,
              p.value, singular)
})
write_csv(stratified, file.path(OUT, "visitor_stratified_slopes.csv"))

# Activity-state subgroups (previously deferred; codebook confirmed 0716:
# 1 staying / 2 arriving / 3 leaving). staying vs in-transit moderation
# (FDR = 7) plus slopes stratified by the three states.
da <- d0 %>%
  mutate(staying = case_when(activity_state == 1 ~ 1L,
                             activity_state %in% c(2, 3) ~ 0L,
                             TRUE ~ NA_integer_),
         act_state = factor(activity_state, levels = 1:3,
                            labels = c("staying", "arriving", "leaving"))) %>%
  filter(!is.na(staying))
cat("activity subset n =", nrow(da), "| staying:", sum(da$staying == 1),
    "| in transit:", sum(da$staying == 0), "\n")

act_moderation <- EFFECTS %>% rowwise() %>%
  mutate(res = list({
    tn <- term
    f <- as.formula(paste0(outcome, " ~ ", tn, " * staying + wave",
                           " + (1 | site)"))
    fit <- suppressMessages(lmer(f, data = da))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(term == paste0(tn, ":staying"))
    tibble(b_interaction = tt$estimate, se = tt$std.error,
           p.value = tt$p.value, singular = isSingular(fit))
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q.value = fdr(p.value)) %>% arrange(q.value)
write_csv(act_moderation, file.path(OUT, "activity_state_moderation.csv"))

act_stratified <- map_dfr(levels(da$act_state), function(g) {
  ds <- filter(da, act_state == g)
  EFFECTS %>% rowwise() %>%
    mutate(res = list({
      tn <- term
      f <- as.formula(paste0(outcome, " ~ ", tn, " + wave + (1 | site)"))
      fit <- tryCatch(suppressMessages(lmer(f, data = ds)),
                      error = function(e) NULL)
      if (is.null(fit)) tibble(estimate = NA_real_, std.error = NA_real_,
                               p.value = NA_real_, singular = NA)
      else {
        tt <- broom.mixed::tidy(fit, effects = "fixed") %>% filter(term == tn)
        tibble(estimate = tt$estimate, std.error = tt$std.error,
               p.value = tt$p.value, singular = isSingular(fit))
      }
    })) %>% ungroup() %>% unnest(res) %>%
    transmute(group = g, n = nrow(ds), effect, estimate, std.error,
              p.value, singular)
})
write_csv(act_stratified, file.path(OUT, "activity_stratified_slopes.csv"))

print(moderation)
print(stratified %>% arrange(effect, group), n = 30)
print(act_moderation)
print(act_stratified %>% arrange(effect, group), n = 30)
