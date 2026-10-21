# Analysis 18 - Greenery internal composition as a direct quality measure
# (follow-up to A17's quantity-without-quality finding). Holding the AMOUNT of
# greenery fixed, does what the greenery is MADE OF (tree / grass / herb-layer
# shares, flower presence, strata layering, vegetation diversity) predict
# evaluation? comp_greenery = tree + grass + plant + palm + flower + field.
# Vegetated rows only (comp_greenery > 0.05). FDR across the quality battery.

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_18_greenery_composition"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

VEG <- c("ph_tree", "ph_grass", "ph_plant", "ph_palm", "ph_flower", "ph_field")

d <- load_vis_data() %>%
  mutate(across(all_of(VEG), ~ replace_na(.x, 0))) %>%
  filter(comp_greenery > 0.05) %>%
  mutate(
    veg_total = ph_tree + ph_grass + ph_plant + ph_palm + ph_flower + ph_field,
    tree_frac = ph_tree / veg_total,
    grass_frac = ph_grass / veg_total,
    herb_frac = (ph_plant + ph_flower) / veg_total,
    flower_in_green = ph_flower / veg_total,
    strata_n = (ph_tree > 0.01) + (ph_grass > 0.01) + (ph_plant + ph_flower > 0.01),
    layered = as.integer(strata_n >= 3),
    veg_shannon = {
      p <- cbind(ph_tree, ph_grass, ph_plant + ph_flower + ph_palm, ph_field) / veg_total
      -rowSums(ifelse(p > 0, p * log(p), 0))
    })
cat("Vegetated rows:", nrow(d), "of photo-valid core\n")

QUALITY <- c("tree_frac", "grass_frac", "herb_frac", "flower_in_green",
             "layered", "veg_shannon")
for (v in QUALITY) d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))

OUTCOMES <- c("isop_z", "overall_soundscape_z", "harmony_z", "isoe_z")

# Screen: each quality term, adjusted for greenery AMOUNT ---------------------
screen <- expand_grid(outcome = OUTCOMES, quality = QUALITY) %>%
  rowwise() %>%
  mutate(res = list({
    f <- as.formula(paste0(outcome, " ~ ", quality, "_z + comp_greenery_z",
                           " + (1 | site) + (1 | site:wave)"))
    fit <- suppressMessages(lmer(f, data = d))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(term == paste0(quality, "_z"))
    tibble(estimate = tt$estimate, std.error = tt$std.error,
           p.value = tt$p.value, singular = isSingular(fit))
  })) %>% ungroup() %>% unnest(res) %>%
  mutate(q.value = fdr(p.value)) %>% arrange(q.value)
write_csv(screen, file.path(OUT, "quality_screen.csv"))

# Joint quality block per outcome (LRT beyond amount) -------------------------
joint <- map_dfr(OUTCOMES, function(y) {
  f1 <- as.formula(paste0(y, " ~ ", paste0(QUALITY, "_z", collapse = " + "),
                          " + comp_greenery_z + (1 | site) + (1 | site:wave)"))
  f0 <- as.formula(paste0(y, " ~ comp_greenery_z + (1 | site) + (1 | site:wave)"))
  m1 <- suppressMessages(lmer(f1, data = d, REML = FALSE))
  m0 <- suppressMessages(lmer(f0, data = d, REML = FALSE))
  tibble(outcome = y, p_joint = anova(m0, m1)$`Pr(>Chisq)`[2])
})
write_csv(joint, file.path(OUT, "quality_joint_block.csv"))

# Survivors: LAeq / source adjustment + season interaction (A17 link) ---------
surv <- filter(screen, q.value < 0.05)
checks <- pmap_dfr(list(surv$outcome, surv$quality), function(y, v) {
  base <- paste0(" + comp_greenery_z + (1 | site) + (1 | site:wave)")
  fits <- list(
    laeq_adjusted = paste0(y, " ~ ", v, "_z + laeq_a_z", base),
    source_adjusted = paste0(y, " ~ ", v, "_z + traffic_noise_z + human_sounds_z + natural_sounds_z + other_noise_z", base))
  imap_dfr(fits, function(f, nm) {
    fit <- suppressMessages(lmer(as.formula(f), data = d))
    tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(term == paste0(v, "_z"))
    tibble(outcome = y, quality = v, check = nm, estimate = tt$estimate,
           p.value = tt$p.value)
  })
})
write_csv(checks, file.path(OUT, "survivor_checks.csv"))

print(screen, n = 30); print(joint); print(checks)
