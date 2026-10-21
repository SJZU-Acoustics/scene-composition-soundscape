# Analysis 05 — Visual composition -> soundscape evaluation.
# Outcomes: ISO-P, ISO-E, OA, perceived loudness, harmony. Screen of the 8
# composites (FDR family = 40), wave-adjusted, site RE; then acoustic
# adjustment (LAeq) and the sound-source-adjusted variant for the survivors.

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_05_visual_to_evaluation"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data()
OUTS <- c("isop", "isoe", "overall_soundscape", "perceived_loudness", "harmony")

screen <- map_dfr(OUTS, function(y) {
  map_dfr(COMPS8, function(v) {
    f <- paste0(y, "_z ~ ", v, "_z + wave + (1 | site)")
    fit <- lmer(as.formula(f), data = d)
    tidy_term(fit, paste0(v, "_z")) %>%
      mutate(outcome = y, composite = v, .before = 1)
  })
}) %>% mutate(q_fdr = fdr(p.value), family_size = n())
write_csv(screen, file.path(OUT, "visual_evaluation_screen.csv"))

surv <- screen %>% filter(q_fdr < 0.05)
write_csv(surv, file.path(OUT, "survivors.csv"))

# Acoustic adjustment (LAeq) for survivors, acoustic+photo subsample
da <- d %>% filter(!is.na(laeq_a))
adj_ac <- pmap_dfr(surv %>% select(outcome, composite), function(outcome, composite) {
  f <- paste0(outcome, "_z ~ ", composite, "_z + laeq_a_z + wave + (1 | site)")
  fit <- lmer(as.formula(f), data = da)
  tidy_term(fit, paste0(composite, "_z")) %>%
    mutate(outcome = outcome, composite = composite, adjustment = "laeq",
           n = nobs(fit), .before = 1)
})
write_csv(adj_ac, file.path(OUT, "survivors_laeq_adjusted.csv"))

# Sound-source adjustment (does vision matter beyond what people hear?)
adj_src <- pmap_dfr(surv %>% select(outcome, composite), function(outcome, composite) {
  f <- paste0(outcome, "_z ~ ", composite,
              "_z + traffic_noise_z + human_sounds_z + natural_sounds_z +",
              " other_noise_z + wave + (1 | site)")
  fit <- lmer(as.formula(f), data = d)
  tidy_term(fit, paste0(composite, "_z")) %>%
    mutate(outcome = outcome, composite = composite, adjustment = "sources",
           n = nobs(fit), .before = 1)
})
write_csv(adj_src, file.path(OUT, "survivors_source_adjusted.csv"))

print(surv %>% select(outcome, composite, estimate, q_fdr), n = 40)
cat("--- LAeq-adjusted survivors:\n")
print(adj_ac %>% select(outcome, composite, estimate, p.value), n = 40)
cat("--- source-adjusted survivors:\n")
print(adj_src %>% select(outcome, composite, estimate, p.value), n = 40)
