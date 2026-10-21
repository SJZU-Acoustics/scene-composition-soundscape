# Analysis 15 — Explicit audio-visual match and conflict indices.
# Product indices are tested beyond their auditory and visual constituents, so
# an index effect means the combination matters, not just the two main effects.

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_15_congruence_conflict"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(
    comp_traffic_vis = comp_vehicles + ph_road,
    comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
    audio_activity = traffic_noise + human_sounds,
    visual_activity = comp_traffic_vis + comp_people,
    audio_activity_z = as.numeric(scale(audio_activity)),
    visual_activity_z = as.numeric(scale(visual_activity)),
    audio_nature_balance = natural_sounds - (traffic_noise + other_noise) / 2,
    visual_nature_balance = comp_nature - comp_traffic_vis,
    audio_nature_balance_z = as.numeric(scale(audio_nature_balance)),
    visual_nature_balance_z = as.numeric(scale(visual_nature_balance)),
    idx_traffic_match = traffic_noise_z * comp_traffic_vis_z,
    idx_nature_match = natural_sounds_z * comp_nature_z,
    idx_human_match = human_sounds_z * comp_people_z,
    idx_traffic_in_nature = traffic_noise_z * comp_nature_z,
    idx_nature_in_traffic = natural_sounds_z * comp_traffic_vis_z,
    idx_activity_match = audio_activity_z * visual_activity_z,
    idx_balance_match = audio_nature_balance_z * visual_nature_balance_z
  )

INDEX_SPECS <- tribble(
  ~index, ~label, ~controls, ~reading,
  "idx_traffic_match", "heard-and-seen traffic",
  list(c("traffic_noise_z", "comp_traffic_vis_z")),
  "traffic sound combined with visible traffic/road",
  "idx_nature_match", "heard-and-seen nature",
  list(c("natural_sounds_z", "comp_nature_z")),
  "natural sound combined with visible greenery/water",
  "idx_human_match", "heard-and-seen people",
  list(c("human_sounds_z", "comp_people_z")),
  "human sound combined with visible people",
  "idx_traffic_in_nature", "traffic-in-nature conflict",
  list(c("traffic_noise_z", "comp_nature_z")),
  "traffic sound in a visually natural frame",
  "idx_nature_in_traffic", "nature-in-traffic conflict",
  list(c("natural_sounds_z", "comp_traffic_vis_z")),
  "natural sound in a traffic/road visual frame",
  "idx_activity_match", "activity match",
  list(c("audio_activity_z", "visual_activity_z")),
  "heard activity combined with visible activity/infrastructure",
  "idx_balance_match", "nature-hard balance match",
  list(c("audio_nature_balance_z", "visual_nature_balance_z")),
  "alignment between auditory natural-vs-noise balance and visual natural-vs-traffic balance"
) %>%
  mutate(index_z = paste0(index, "_z"))

for (v in INDEX_SPECS$index) {
  d[[paste0(v, "_z")]] <- as.numeric(scale(d[[v]]))
}

defs <- INDEX_SPECS %>%
  mutate(controls = map_chr(controls, ~ paste(.x, collapse = " + ")))
write_csv(defs, file.path(OUT, "index_definitions.csv"))

OUTS <- c("isop", "isoe", "overall_soundscape", "harmony", "raq_mean")

screen <- map_dfr(OUTS, function(y) {
  pmap_dfr(INDEX_SPECS, function(index, label, controls, reading, index_z) {
    controls_vec <- unlist(controls)
    f <- as.formula(paste0(y, "_z ~ ", index_z, " + ",
                           paste(controls_vec, collapse = " + "),
                           " + wave + (1 | site)"))
    fit <- lmer(f, data = d)
    tidy_term(fit, index_z) %>%
      mutate(outcome = y, index = index, label = label,
             controls = paste(controls_vec, collapse = " + "),
             .before = 1)
  })
}) %>%
  mutate(q_fdr = fdr(p.value), family_size = n()) %>%
  arrange(q_fdr, p.value)
write_csv(screen, file.path(OUT, "match_conflict_screen.csv"))

surv <- screen %>% filter(q_fdr < 0.05)
write_csv(surv, file.path(OUT, "match_conflict_survivors.csv"))

# Targeted existing contrast for comparison with analysis 06: traffic match on OA.
target <- screen %>%
  filter(index == "idx_traffic_match",
         outcome %in% c("overall_soundscape", "isop", "isoe", "harmony"))
write_csv(target, file.path(OUT, "traffic_match_targeted.csv"))

loso <- if (nrow(surv) > 0) {
  pmap_dfr(surv %>% select(outcome, index, term, controls), function(outcome, index, term, controls) {
    controls_vec <- str_split(controls, " \\+ ", simplify = FALSE)[[1]]
    map_dfr(unique(as.character(d$site)), function(s) {
      dd <- filter(d, site != s)
      f <- as.formula(paste0(outcome, "_z ~ ", term, " + ",
                             paste(controls_vec, collapse = " + "),
                             " + wave + (1 | site)"))
      est <- tryCatch({
        fit <- lmer(f, data = dd)
        tidy_term(fit, term)
      }, error = function(e) tibble(estimate = NA_real_, p.value = NA_real_))
      tibble(outcome = outcome, index = index, dropped_site = s,
             estimate = est$estimate[1], p.value = est$p.value[1])
    })
  })
} else {
  tibble()
}
write_csv(loso, file.path(OUT, "loso.csv"))

loso_summary <- if (nrow(loso) > 0) {
  loso %>% group_by(outcome, index) %>%
    summarise(min_est = min(estimate, na.rm = TRUE),
              max_est = max(estimate, na.rm = TRUE),
              sign_flips = sum(sign(estimate) != sign(median(estimate, na.rm = TRUE)),
                               na.rm = TRUE),
              max_p = max(p.value, na.rm = TRUE),
              .groups = "drop")
} else {
  tibble()
}
write_csv(loso_summary, file.path(OUT, "loso_summary.csv"))

boot <- if (nrow(surv) > 0) {
  pmap_dfr(surv %>% select(outcome, index, term, controls), function(outcome, index, term, controls) {
    controls_vec <- str_split(controls, " \\+ ", simplify = FALSE)[[1]]
    sites <- unique(as.character(d$site))
    ests <- map_dbl(seq_len(1000), function(b) {
      samp <- sample(sites, length(sites), replace = TRUE)
      db <- map_dfr(seq_along(samp), function(i)
        d %>% filter(site == samp[i]) %>% mutate(site_b = paste0("s", i)))
      f <- as.formula(paste0(outcome, "_z ~ ", term, " + ",
                             paste(controls_vec, collapse = " + "),
                             " + wave + site_b"))
      fit <- tryCatch(lm(f, data = db), error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      tt <- broom::tidy(fit) %>% filter(.data$term == !!term)
      if (nrow(tt)) tt$estimate[1] else NA_real_
    })
    es <- ests[!is.na(ests)]
    tibble(outcome = outcome, index = index, n_boot = length(es),
           ci_lo = quantile(es, 0.025), ci_hi = quantile(es, 0.975),
           excludes_zero = sign(quantile(es, 0.025)) == sign(quantile(es, 0.975)))
  })
} else {
  tibble()
}
write_csv(boot, file.path(OUT, "cluster_bootstrap.csv"))

cat("--- top screen rows ---\n")
print(screen %>% select(outcome, label, estimate, p.value, q_fdr) %>% head(20),
      n = 20)
cat("\n--- survivors ---\n")
print(surv %>% select(outcome, label, estimate, p.value, q_fdr), n = 50)
cat("\n--- traffic match target ---\n")
print(target %>% select(outcome, estimate, p.value, q_fdr), n = 10)
cat("\n--- LOSO summary ---\n")
print(loso_summary, n = 20)
cat("\n--- bootstrap ---\n")
print(boot, n = 20)
