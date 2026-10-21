# Analysis 13 — Broader context composites from semantic segmentation.
# Goal: test whether constructed context parameters (not only content shares)
# explain the evaluative layer beyond site/wave, and whether any signal remains
# beyond the existing 8 visual content composites or the sound battery.

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_13_context_composites"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d0 <- load_vis_data()
ph_cols <- names(d0)[stringr::str_starts(names(d0), "ph_")]

col0 <- function(dd, nm) {
  if (nm %in% names(dd)) replace_na(dd[[nm]], 0) else rep(0, nrow(dd))
}

ph <- d0[ph_cols] %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  as.matrix()
rs <- rowSums(ph)
phn <- ph / rs
phn[!is.finite(phn)] <- 0

entropy <- apply(phn, 1, function(x) {
  xp <- x[x > 0]
  -sum(xp * log(xp))
})
richness_001 <- rowSums(phn >= 0.01)
dominance <- apply(phn, 1, max)
effective_n <- exp(entropy)
evenness <- ifelse(richness_001 > 1, entropy / log(richness_001), NA_real_)

natural_mass <- d0$comp_greenery_rn + d0$comp_water_rn + d0$comp_terrain_rn
hard_mass <- d0$comp_paving_rn + d0$comp_built_rn + d0$comp_vehicles_rn
natural_hard_total <- natural_mass + hard_mass
natural_hard_mixedness <- ifelse(
  natural_hard_total > 0,
  1 - abs(natural_mass - hard_mass) / natural_hard_total,
  NA_real_
)

open_horizon <- col0(d0, "ph_sky") + col0(d0, "ph_mountain") +
  col0(d0, "ph_water") + col0(d0, "ph_sea") + col0(d0, "ph_lake") +
  col0(d0, "ph_river") - col0(d0, "ph_building") - col0(d0, "ph_wall") -
  col0(d0, "ph_fence") - col0(d0, "ph_railing")

mobility_infrastructure <- col0(d0, "ph_road") + d0$comp_vehicles_rn

social_amenity <- col0(d0, "ph_person") + col0(d0, "ph_bench") +
  col0(d0, "ph_chair") + col0(d0, "ph_seat") + col0(d0, "ph_table") +
  col0(d0, "ph_plaything")

restorative_cues <- d0$comp_water_rn + col0(d0, "ph_flower") +
  col0(d0, "ph_fountain") + col0(d0, "ph_bench")

d <- d0 %>%
  mutate(
    ctx_diversity = entropy,
    ctx_effective_n = effective_n,
    ctx_richness_001 = richness_001,
    ctx_evenness = evenness,
    ctx_dominance = dominance,
    ctx_natural_hard_mixedness = natural_hard_mixedness,
    ctx_open_horizon = open_horizon,
    ctx_mobility_infrastructure = mobility_infrastructure,
    ctx_social_amenity = social_amenity,
    ctx_restorative_cues = restorative_cues
  )

CTX_MAIN <- c("ctx_diversity", "ctx_dominance", "ctx_natural_hard_mixedness",
              "ctx_open_horizon", "ctx_mobility_infrastructure",
              "ctx_social_amenity", "ctx_restorative_cues")
CTX_DIAG <- c("ctx_diversity", "ctx_effective_n", "ctx_richness_001",
              "ctx_evenness", "ctx_dominance")
OUTS <- c("isop", "isoe", "overall_soundscape", "harmony", "raq_mean")
SOUND <- c("traffic_noise_z", "human_sounds_z", "natural_sounds_z",
           "other_noise_z", "perceived_loudness_z")
VISUAL <- paste0(COMPS8, "_z")

for (v in c(CTX_MAIN, CTX_DIAG)) {
  d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
}

definitions <- tribble(
  ~parameter, ~label, ~construction, ~intended_reading,
  "ctx_diversity", "Visual diversity",
  "Shannon entropy of the 150 row-normalised ADE20K shares",
  "Semantic richness / compositional complexity of the visual frame",
  "ctx_dominance", "Dominant-surface share",
  "Largest row-normalised ADE20K category share",
  "Visual simplicity or dominance by one surface/object class",
  "ctx_natural_hard_mixedness", "Natural-hard mixedness",
  "1 - |natural - hard| / (natural + hard), with natural = greenery + water + terrain and hard = paving + built + vehicles",
  "Balanced mosaic rather than one-sided natural or hardscape frame",
  "ctx_open_horizon", "Open horizon",
  "sky + mountain + water/sea/lake/river - building/wall/fence/railing",
  "Open visual field versus enclosed built edge",
  "ctx_mobility_infrastructure", "Mobility infrastructure",
  "road + vehicles",
  "Visible traffic/infrastructure cue, broader than vehicles alone",
  "ctx_social_amenity", "Social amenity",
  "person + bench + chair + seat + table + plaything",
  "Visible people and use-supporting amenities",
  "ctx_restorative_cues", "Restorative cues",
  "water + flower + fountain + bench",
  "Salient designed/restorative elements"
)
write_csv(definitions, file.path(OUT, "context_composite_definitions.csv"))

desc <- map_dfr(c(CTX_MAIN, CTX_DIAG), function(v) {
  tibble(parameter = v, n = sum(!is.na(d[[v]])),
         mean = mean(d[[v]], na.rm = TRUE),
         sd = sd(d[[v]], na.rm = TRUE),
         p05 = quantile(d[[v]], 0.05, na.rm = TRUE),
         median = median(d[[v]], na.rm = TRUE),
         p95 = quantile(d[[v]], 0.95, na.rm = TRUE))
})
write_csv(desc, file.path(OUT, "context_composite_descriptives.csv"))

write_csv(as_tibble(cor(d[CTX_DIAG], use = "pairwise"), rownames = "parameter"),
          file.path(OUT, "complexity_proxy_correlations.csv"))
write_csv(as_tibble(cor(d[CTX_MAIN], use = "pairwise"), rownames = "parameter"),
          file.path(OUT, "context_composite_correlations.csv"))

screen <- map_dfr(OUTS, function(y) {
  map_dfr(CTX_MAIN, function(v) {
    f <- paste0(y, "_z ~ ", v, "_z + wave + (1 | site)")
    fit <- lmer(as.formula(f), data = d)
    tidy_term(fit, paste0(v, "_z")) %>%
      mutate(outcome = y, parameter = v, .before = 1)
  })
}) %>% mutate(q_fdr = fdr(p.value), family_size = n())
write_csv(screen, file.path(OUT, "context_screen.csv"))

surv <- screen %>% filter(q_fdr < 0.05) %>%
  arrange(q_fdr, outcome, parameter)
write_csv(surv, file.path(OUT, "context_survivors.csv"))

fit_context_term <- function(y, v, extras = character()) {
  rhs <- c(paste0(v, "_z"), extras, "wave", "(1 | site)")
  f <- as.formula(paste0(y, "_z ~ ", paste(rhs, collapse = " + ")))
  fit <- lmer(f, data = d)
  tidy_term(fit, paste0(v, "_z"))
}

surv_adjusted <- if (nrow(surv) > 0) {
  pmap_dfr(surv %>% select(outcome, parameter), function(outcome, parameter) {
    bind_rows(
      fit_context_term(outcome, parameter, VISUAL) %>%
        mutate(adjustment = "existing_visual_composites"),
      fit_context_term(outcome, parameter, SOUND) %>%
        mutate(adjustment = "sound_battery"),
      fit_context_term(outcome, parameter, c(VISUAL, SOUND)) %>%
        mutate(adjustment = "visual_plus_sound")
    ) %>%
      mutate(outcome = outcome, parameter = parameter, .before = 1)
  })
} else {
  tibble(outcome = character(), parameter = character(), adjustment = character(),
         term = character(), estimate = numeric(), std.error = numeric(),
         statistic = numeric(), df = numeric(), p.value = numeric())
}
write_csv(surv_adjusted, file.path(OUT, "survivors_adjusted.csv"))

partial_r2 <- function(y, extra, dd) {
  f0 <- as.formula(paste0(y, "_z ~ site + wave"))
  f1 <- as.formula(paste0(y, "_z ~ site + wave + ",
                          paste(extra, collapse = " + ")))
  m0 <- lm(f0, data = dd); m1 <- lm(f1, data = dd)
  1 - sum(resid(m1)^2) / sum(resid(m0)^2)
}

joint <- map_dfr(OUTS, function(y) {
  dd <- d %>% drop_na(all_of(c(paste0(y, "_z"), paste0(CTX_MAIN, "_z"),
                              VISUAL, SOUND)))
  r_ctx <- partial_r2(y, paste0(CTX_MAIN, "_z"), dd)
  r_visual <- partial_r2(y, VISUAL, dd)
  r_sound <- partial_r2(y, SOUND, dd)
  r_visual_ctx <- partial_r2(y, c(VISUAL, paste0(CTX_MAIN, "_z")), dd)
  r_sound_ctx <- partial_r2(y, c(SOUND, paste0(CTX_MAIN, "_z")), dd)
  tibble(outcome = y, n = nrow(dd),
         context_total = r_ctx,
         existing_visual_total = r_visual,
         sound_total = r_sound,
         context_unique_beyond_visual = r_visual_ctx - r_visual,
         context_unique_beyond_sound = r_sound_ctx - r_sound)
})
write_csv(joint, file.path(OUT, "joint_context_increment.csv"))

joint_tests <- map_dfr(OUTS, function(y) {
  dd <- d %>% drop_na(all_of(c(paste0(y, "_z"), paste0(CTX_MAIN, "_z"),
                              VISUAL, SOUND)))
  yv <- paste0(y, "_z")
  base <- lm(as.formula(paste0(yv, " ~ site + wave")), data = dd)
  ctx <- lm(as.formula(paste0(yv, " ~ site + wave + ",
                              paste(paste0(CTX_MAIN, "_z"), collapse = " + "))),
            data = dd)
  visual <- lm(as.formula(paste0(yv, " ~ site + wave + ",
                                 paste(VISUAL, collapse = " + "))),
               data = dd)
  visual_ctx <- lm(as.formula(paste0(yv, " ~ site + wave + ",
                                     paste(c(VISUAL, paste0(CTX_MAIN, "_z")),
                                           collapse = " + "))),
                   data = dd)
  sound <- lm(as.formula(paste0(yv, " ~ site + wave + ",
                                paste(SOUND, collapse = " + "))),
              data = dd)
  sound_ctx <- lm(as.formula(paste0(yv, " ~ site + wave + ",
                                    paste(c(SOUND, paste0(CTX_MAIN, "_z")),
                                          collapse = " + "))),
                  data = dd)
  tibble(
    outcome = y,
    n = nrow(dd),
    context_vs_base_p = anova(base, ctx)$`Pr(>F)`[2],
    context_beyond_visual_p = anova(visual, visual_ctx)$`Pr(>F)`[2],
    context_beyond_sound_p = anova(sound, sound_ctx)$`Pr(>F)`[2]
  )
})
write_csv(joint_tests, file.path(OUT, "joint_context_tests.csv"))

# Robustness only for FDR survivors. Cluster bootstrap is site resampling.
robust_set <- surv %>% select(outcome, parameter)

loso <- if (nrow(robust_set) > 0) {
  pmap_dfr(robust_set, function(outcome, parameter) {
    map_dfr(unique(as.character(d$site)), function(s) {
      dd <- filter(d, site != s)
      f <- as.formula(paste0(outcome, "_z ~ ", parameter, "_z + wave + (1 | site)"))
      est <- tryCatch({
        fit <- lmer(f, data = dd)
        tidy_term(fit, paste0(parameter, "_z"))
      }, error = function(e) tibble(estimate = NA_real_, p.value = NA_real_))
      tibble(outcome = outcome, parameter = parameter, dropped_site = s,
             estimate = est$estimate[1], p.value = est$p.value[1])
    })
  })
} else {
  tibble()
}
write_csv(loso, file.path(OUT, "loso.csv"))

loso_sum <- if (nrow(loso) > 0) {
  loso %>% group_by(outcome, parameter) %>%
    summarise(min_est = min(estimate, na.rm = TRUE),
              max_est = max(estimate, na.rm = TRUE),
              sign_flips = sum(sign(estimate) != sign(median(estimate, na.rm = TRUE)),
                               na.rm = TRUE),
              max_p = max(p.value, na.rm = TRUE),
              .groups = "drop")
} else {
  tibble()
}
write_csv(loso_sum, file.path(OUT, "loso_summary.csv"))

boot <- if (nrow(robust_set) > 0) {
  pmap_dfr(robust_set, function(outcome, parameter) {
    sites <- unique(as.character(d$site))
    ests <- map_dbl(seq_len(1000), function(b) {
      samp <- sample(sites, length(sites), replace = TRUE)
      db <- map_dfr(seq_along(samp), function(i)
        d %>% filter(site == samp[i]) %>% mutate(site_b = paste0("s", i)))
      f <- as.formula(paste0(outcome, "_z ~ ", parameter,
                             "_z + wave + (1 | site_b)"))
      fit <- tryCatch(lmer(f, data = db), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      tt <- broom.mixed::tidy(fit, effects = "fixed") %>%
        filter(term == paste0(parameter, "_z"))
      if (nrow(tt)) tt$estimate[1] else NA_real_
    })
    es <- ests[!is.na(ests)]
    tibble(outcome = outcome, parameter = parameter, n_boot = length(es),
           ci_lo = quantile(es, 0.025), ci_hi = quantile(es, 0.975),
           excludes_zero = sign(quantile(es, 0.025)) == sign(quantile(es, 0.975)))
  })
} else {
  tibble()
}
write_csv(boot, file.path(OUT, "cluster_bootstrap.csv"))

print(definitions, n = nrow(definitions))
cat("\n--- FDR survivors ---\n")
print(surv %>% select(outcome, parameter, estimate, p.value, q_fdr), n = 50)
cat("\n--- adjusted survivors ---\n")
if (nrow(surv_adjusted) > 0) {
  print(surv_adjusted %>% select(outcome, parameter, adjustment, estimate, p.value),
        n = 80)
} else {
  cat("No FDR survivors; adjusted survivor checks not run.\n")
}
cat("\n--- joint increments ---\n")
print(joint, n = 10)
cat("\n--- joint F-tests ---\n")
print(joint_tests, n = 10)
cat("\n--- LOSO summary ---\n")
print(loso_sum, n = 80)
cat("\n--- cluster bootstrap ---\n")
print(boot, n = 80)
