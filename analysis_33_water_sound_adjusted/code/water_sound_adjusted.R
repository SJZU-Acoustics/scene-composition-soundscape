# Analysis 33 - Water-in-view effect adjusted for the auditory measures
# (post-review robustness, 2026-09-18). The vehicle and sky-blueness effects
# already carry a refit with what is heard in the model; the water effect did
# not, so "visible water" could stand in for "audible water" (water sound is
# rated within the natural-sound saliency). Refit of the A11 water models with
# the four perceived source saliencies and perceived loudness added. No new
# test: robustness refit of an effect already significant in the element family.

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_33_water_sound_adjusted"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SOUND <- c("traffic_noise_z", "human_sounds_z", "natural_sounds_z",
           "other_noise_z", "perceived_loudness_z")
d <- load_vis_data() %>% mutate(water_any = as.integer(comp_water > 0.005))

fit_water <- function(outcome, covars, label) {
  dd <- d %>% drop_na(all_of(c(outcome, covars)))
  fit <- lmer(reformulate(c("water_any", covars, "wave", "(1 | site)"), outcome), data = dd)
  broom.mixed::tidy(fit, effects = "fixed") %>% filter(term == "water_any") %>%
    transmute(outcome = outcome, specification = label, estimate, se = std.error,
              p.value, n = nobs(fit))
}
res <- map_dfr(c("overall_soundscape_z", "isop_z"), function(oc) bind_rows(
  fit_water(oc, character(0), "as reported"),
  fit_water(oc, "natural_sounds_z", "natural-sound saliency adjusted"),
  fit_water(oc, SOUND, "all auditory measures adjusted")))
write_csv(res, file.path(OUT, "water_sound_adjusted.csv"))
print(res, width = Inf)

# The water x perceived-loudness moderation (A19 specification, acoustic
# subsample, site + site:wave random effects) with the source saliencies added.
da <- d %>% filter(!is.na(laeq_a))
SALIENCY <- setdiff(SOUND, "perceived_loudness_z")
mod <- map_dfr(list(`as reported` = character(0), `source saliencies adjusted` = SALIENCY),
               function(covars) {
  fit <- lmer(reformulate(c("water_any * perceived_loudness_z", covars,
                            "(1 | site)", "(1 | site:wave)"), "overall_soundscape_z"), data = da)
  broom.mixed::tidy(fit, effects = "fixed") %>%
    filter(term == "water_any:perceived_loudness_z") %>%
    transmute(estimate, se = std.error, p.value, n = nobs(fit))
}, .id = "specification")
write_csv(mod, file.path(OUT, "water_loudness_sound_adjusted.csv"))
print(mod, width = Inf)

# Does natural-sound saliency differ with water in view? (descriptive)
nat <- lmer(natural_sounds_z ~ water_any + wave + (1 | site), data = d)
nat_tab <- broom.mixed::tidy(nat, effects = "fixed") %>% filter(term == "water_any") %>%
  transmute(outcome = "natural_sounds_z", estimate, se = std.error, p.value, n = nobs(nat))
write_csv(nat_tab, file.path(OUT, "natural_sound_by_water.csv"))
print(nat_tab, width = Inf)
