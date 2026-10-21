# Analysis 16 — Eventfulness pathway.
# Associational pathway test: do visual/cross-modal cues primarily affect
# ISO-E, and does ISO-E carry that signal onward to hedonic evaluation?
# Fixed site + wave effects are used for stable product-of-coefficients
# estimation; bootstrap resamples sites 1,000 times for indirect effects.

source("_shared/helpers.R")
set.seed(20260709)

ANALYSIS_DIR <- "analysis_16_eventfulness_pathway"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(
    comp_traffic_vis = comp_vehicles + ph_road,
    comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)),
    visual_activity = comp_traffic_vis + comp_people,
    visual_activity_z = as.numeric(scale(visual_activity)),
    social_amenity = replace_na(ph_person, 0) + replace_na(ph_bench, 0) +
      replace_na(ph_chair, 0) + replace_na(ph_seat, 0) +
      replace_na(ph_table, 0) + replace_na(ph_plaything, 0),
    social_amenity_z = as.numeric(scale(social_amenity)),
    audio_nature_balance = natural_sounds - (traffic_noise + other_noise) / 2,
    visual_nature_balance = comp_nature - comp_traffic_vis,
    audio_nature_balance_z = as.numeric(scale(audio_nature_balance)),
    visual_nature_balance_z = as.numeric(scale(visual_nature_balance)),
    idx_traffic_match = traffic_noise_z * comp_traffic_vis_z,
    idx_traffic_match_z = as.numeric(scale(idx_traffic_match)),
    idx_traffic_in_nature = traffic_noise_z * comp_nature_z,
    idx_traffic_in_nature_z = as.numeric(scale(idx_traffic_in_nature)),
    idx_balance_match = audio_nature_balance_z * visual_nature_balance_z,
    idx_balance_match_z = as.numeric(scale(idx_balance_match))
  )

CUES <- tribble(
  ~cue, ~label, ~controls, ~family,
  "comp_vehicles_z", "visible vehicles", list(character()), "visual",
  "comp_traffic_vis_z", "visible traffic/road", list(character()), "visual",
  "comp_people_z", "visible people", list(character()), "visual",
  "visual_activity_z", "visible activity", list(character()), "visual",
  "social_amenity_z", "social amenity", list(character()), "visual",
  "idx_traffic_match_z", "heard-and-seen traffic",
  list(c("traffic_noise_z", "comp_traffic_vis_z")), "cross_modal",
  "idx_traffic_in_nature_z", "traffic-in-nature conflict",
  list(c("traffic_noise_z", "comp_nature_z")), "cross_modal",
  "idx_balance_match_z", "nature-hard balance match",
  list(c("audio_nature_balance_z", "visual_nature_balance_z")), "cross_modal"
)

OUTS <- tribble(
  ~outcome, ~label,
  "isop_z", "ISO pleasantness",
  "overall_soundscape_z", "Overall evaluation",
  "harmony_z", "Harmony"
)

lm_coef <- function(formula, term, dd) {
  fit <- lm(as.formula(formula), data = dd)
  broom::tidy(fit) %>% filter(.data$term == !!term)
}

rhs <- function(primary, controls, extra = character(), site_var = "site") {
  vars <- c(primary, controls, extra, "wave", site_var)
  paste(vars, collapse = " + ")
}

a_paths <- pmap_dfr(CUES, function(cue, label, controls, family) {
  controls_vec <- unlist(controls)
  f <- paste0("isoe_z ~ ", rhs(cue, controls_vec))
  lm_coef(f, cue, d) %>%
    mutate(cue = cue, label = label, family = family,
           controls = paste(controls_vec, collapse = " + "), .before = 1)
}) %>%
  mutate(q_fdr = fdr(p.value), family_size = n()) %>%
  arrange(q_fdr, p.value)
write_csv(a_paths, file.path(OUT, "eventfulness_a_paths.csv"))

candidate_cues <- a_paths %>% filter(q_fdr < 0.05) %>% pull(cue)

pathways <- pmap_dfr(CUES, function(cue, label, controls, family) {
  controls_vec <- unlist(controls)
  map_dfr(seq_len(nrow(OUTS)), function(i) {
    y <- OUTS$outcome[i]
    ylab <- OUTS$label[i]
    total_f <- paste0(y, " ~ ", rhs(cue, controls_vec))
    a_f <- paste0("isoe_z ~ ", rhs(cue, controls_vec))
    med_f <- paste0(y, " ~ ", rhs(cue, controls_vec, extra = "isoe_z"))
    total_row <- lm_coef(total_f, cue, d)
    a_row <- lm_coef(a_f, cue, d)
    b_row <- lm_coef(med_f, "isoe_z", d)
    direct_row <- lm_coef(med_f, cue, d)
    tibble(cue = cue, label = label, family = family,
           controls = paste(controls_vec, collapse = " + "),
           outcome = y, outcome_label = ylab,
           total = total_row$estimate[1], total_p = total_row$p.value[1],
           a_eventfulness = a_row$estimate[1], a_p = a_row$p.value[1],
           b_eventfulness_to_outcome = b_row$estimate[1], b_p = b_row$p.value[1],
           direct = direct_row$estimate[1], direct_p = direct_row$p.value[1],
           indirect = a_row$estimate[1] * b_row$estimate[1])
  })
})
write_csv(pathways, file.path(OUT, "eventfulness_pathways.csv"))

# Bootstrap only the main overall-evaluation pathway for the strongest
# eventfulness cues. The full point-estimate table is still written for all
# cue-outcome combinations; bootstrapping every secondary combination is slow
# and adds little interpretive value.
boot_focus <- c("comp_vehicles_z", "idx_traffic_match_z", "idx_balance_match_z")
boot_targets <- pathways %>%
  filter(cue %in% intersect(candidate_cues, boot_focus),
         outcome == "overall_soundscape_z") %>%
  select(cue, label, controls, outcome, outcome_label)

boot <- if (nrow(boot_targets) > 0) {
  pmap_dfr(boot_targets, function(cue, label, controls, outcome, outcome_label) {
    controls_vec <- if (controls == "") character() else str_split(controls, " \\+ ")[[1]]
    sites <- unique(as.character(d$site))
    vals <- map_dbl(seq_len(1000), function(bi) {
      samp <- sample(sites, length(sites), replace = TRUE)
      db <- map_dfr(seq_along(samp), function(i)
        d %>% filter(site == samp[i]) %>% mutate(site_b = paste0("s", i)))
      a_f <- paste0("isoe_z ~ ", rhs(cue, controls_vec, site_var = "site_b"))
      med_f <- paste0(outcome, " ~ ", rhs(cue, controls_vec, extra = "isoe_z",
                                          site_var = "site_b"))
      a_fit <- tryCatch(lm(as.formula(a_f), data = db), error = function(e) NULL)
      b_fit <- tryCatch(lm(as.formula(med_f), data = db), error = function(e) NULL)
      if (is.null(a_fit) || is.null(b_fit)) return(NA_real_)
      a_est <- broom::tidy(a_fit) %>% filter(.data$term == !!cue) %>% pull(estimate)
      b_est <- broom::tidy(b_fit) %>% filter(.data$term == "isoe_z") %>% pull(estimate)
      if (!length(a_est) || !length(b_est)) return(NA_real_)
      a_est[1] * b_est[1]
    })
    vv <- vals[!is.na(vals)]
    tibble(cue = cue, label = label, outcome = outcome,
           outcome_label = outcome_label, n_boot = length(vv),
           indirect_ci_lo = quantile(vv, 0.025),
           indirect_ci_hi = quantile(vv, 0.975),
           indirect_excludes_zero = sign(quantile(vv, 0.025)) ==
             sign(quantile(vv, 0.975)))
  })
} else {
  tibble()
}
write_csv(boot, file.path(OUT, "eventfulness_indirect_bootstrap.csv"))

summary_table <- pathways %>%
  left_join(boot, by = c("cue", "label", "outcome", "outcome_label")) %>%
  arrange(outcome, desc(abs(a_eventfulness)))
write_csv(summary_table, file.path(OUT, "eventfulness_pathway_summary.csv"))

cat("--- Eventfulness a-paths ---\n")
print(a_paths %>% select(label, estimate, p.value, q_fdr), n = 20)
cat("\n--- Pathway summary for a-path survivors ---\n")
print(summary_table %>%
        filter(cue %in% candidate_cues) %>%
        select(label, outcome_label, a_eventfulness, b_eventfulness_to_outcome,
               indirect, indirect_ci_lo, indirect_ci_hi, direct, total),
      n = 50)
