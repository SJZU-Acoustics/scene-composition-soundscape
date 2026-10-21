# Analysis 34 - The mapping test: what does scene composition predict at a place
# the model has not seen? (2026-10-16, Prof Zhang: "do whatever analysis you
# would do"; targeted at the indicator question of the npj package.)
#
# Street-view soundscape maps train on surveyed places and predict evaluation at
# new ones. A28 tested composition within sites (site + wave adjusted inside the
# folds); this analysis tests the mapping use itself, between places:
#
# Part 1 - leave-one-site-out prediction (12 folds). Each block is fitted by OLS
#   with wave fixed effects on 11 sites and predicts the respondents of the
#   held-out site; the baseline is the training mean with the same wave
#   effects, so skill is what the block adds to knowing the survey wave.
#   Blocks: composition (8 composites), perceived auditory measures (4 source
#   saliencies + perceived loudness, as A28), measured acoustics (LAeq, L10-L90,
#   N5, sharpness, roughness, fluctuation strength, tonality TC; the acoustic
#   frame), and their unions. Read-outs per outcome:
#     person-level held-out R2 (1 - SSE / SSE of the wave-only baseline);
#     site-level held-out R2 and Spearman rho of predicted vs observed site means;
#     within-site R2 of the predicted deviations from the site's own mean.
#   Uncertainty: 1,000 site-cluster bootstrap resamples of the per-site
#   contributions (paired across blocks, so block differences get CIs).
#   Sensitivity: the same without wave.
# Part 2 - the proxy check. Variance share of the composition block (OLS R2
#   increment) pooled across places, pooled net of the perceived auditory
#   measures, pooled net of measured acoustics, and within sites (site fixed
#   effects) net of the auditory measures; site-cluster bootstrap CIs (1,000).
# Part 3 - context: share of each outcome's variance lying between sites
#   (null LMM ICC with wave fixed effects).

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_34_site_heldout_mapping"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SEED <- 20261016; N_BOOT <- 1000

OUTS <- c(overall_soundscape_z = "OA", isop_z = "ISO-P", isoe_z = "ISO-E", harmony_z = "Harmony")
VIS <- paste0(COMPS8, "_z")
SOUND <- c("traffic_noise_z", "human_sounds_z", "natural_sounds_z", "other_noise_z",
           "perceived_loudness_z")
MEAS <- c("laeq_a_z", "l10_l90_a_z", "loudness_n5_z", "sharpness_s_z", "roughness_r_z",
          "fluctuation_f_z", "tonality_tc_z")

base <- load_vis_data() %>%
  mutate(site = droplevels(factor(site)), wave = droplevels(factor(wave)))
d_per <- base %>% select(id, site, site_label, space_class, wave,
                         all_of(c(names(OUTS), VIS, SOUND))) %>% drop_na()
d_ac <- base %>% select(id, site, site_label, space_class, wave,
                        all_of(c(names(OUTS), VIS, SOUND, MEAS))) %>% drop_na()
cat("perception frame n =", nrow(d_per), "| acoustic frame n =", nrow(d_ac),
    "| sites =", n_distinct(d_per$site), "\n")

BLOCKS_PER <- list(composition = VIS, auditory = SOUND, auditory_composition = c(SOUND, VIS))
BLOCKS_AC <- list(composition = VIS, measured = MEAS, measured_composition = c(MEAS, VIS),
                  auditory = SOUND, auditory_measured = c(SOUND, MEAS),
                  all = c(SOUND, MEAS, VIS))

# ---- Part 1: leave-one-site-out ------------------------------------------------
# Per held-out site: SSE of the block model and of the wave-only baseline, plus
# observed and predicted site means and the within-site squared errors.
loso <- function(dat, cols, y, wave = TRUE) {
  sites <- levels(dat$site)
  map_dfr(sites, function(s) {
    tr <- dat[dat$site != s, ]; te <- dat[dat$site == s, ]
    rhs0 <- if (wave) "wave" else "1"
    f0 <- reformulate(rhs0, y); f1 <- reformulate(c(rhs0, cols), y)
    te_w <- te; tr$wave <- droplevels(tr$wave)
    stopifnot(all(as.character(te$wave) %in% levels(tr$wave)) || !wave)
    if (wave) te_w$wave <- factor(as.character(te$wave), levels = levels(tr$wave))
    m0 <- lm(f0, tr); m1 <- lm(f1, tr)
    p0 <- predict(m0, te_w); p1 <- predict(m1, te_w)
    yy <- te[[y]]
    tibble(site = s, n = length(yy), sse = sum((yy - p1)^2), sse0 = sum((yy - p0)^2),
           obs_mean = mean(yy), pred_mean = mean(p1), base_mean = mean(p0),
           sse_within = sum(((yy - mean(yy)) - (p1 - mean(p1)))^2),
           sst_within = sum((yy - mean(yy))^2))
  })
}
summ <- function(t) {
  tibble(r2_person = 1 - sum(t$sse) / sum(t$sse0),
         r2_site = 1 - sum((t$obs_mean - t$pred_mean)^2) / sum((t$obs_mean - t$base_mean)^2),
         rho_site = suppressWarnings(cor(t$obs_mean - t$base_mean, t$pred_mean - t$base_mean,
                                         method = "spearman")),
         r2_within = 1 - sum(t$sse_within) / sum(t$sst_within))
}

run_frame <- function(dat, blocks, frame, wave = TRUE) {
  per_site <- imap_dfr(blocks, function(cols, b) {
    map_dfr(names(OUTS), function(y) loso(dat, cols, y, wave) %>%
              mutate(block = b, outcome = y))
  }) %>% mutate(frame = frame, wave_adjusted = wave)
  per_site
}
ps <- bind_rows(run_frame(d_per, BLOCKS_PER, "perception"),
                run_frame(d_ac, BLOCKS_AC, "acoustic"),
                run_frame(d_per, BLOCKS_PER, "perception", wave = FALSE),
                run_frame(d_ac, BLOCKS_AC, "acoustic", wave = FALSE))
write_csv(ps, file.path(OUT, "loso_per_site.csv"))

point <- ps %>% group_by(frame, wave_adjusted, block, outcome) %>%
  group_modify(~ summ(.x)) %>% ungroup() %>% mutate(outcome_label = OUTS[outcome])

# Site-cluster bootstrap of the per-site contributions, paired across blocks.
set.seed(SEED)
site_ids <- levels(d_per$site)
boot_idx <- map(seq_len(N_BOOT), ~ sample(site_ids, length(site_ids), replace = TRUE))
boot <- ps %>% group_by(frame, wave_adjusted, block, outcome) %>%
  group_modify(function(t, k) {
    tt <- t %>% column_to_rownames("site")
    map_dfr(seq_len(N_BOOT), function(b) summ(tt[boot_idx[[b]], ] %>% as_tibble()) %>%
              mutate(boot = b))
  }) %>% ungroup()
ci <- boot %>% group_by(frame, wave_adjusted, block, outcome) %>%
  summarise(across(c(r2_person, r2_site, rho_site, r2_within),
                   list(lo = ~ quantile(.x, 0.025, na.rm = TRUE),
                        hi = ~ quantile(.x, 0.975, na.rm = TRUE))), .groups = "drop")
skill <- point %>% left_join(ci, by = c("frame", "wave_adjusted", "block", "outcome"))
write_csv(skill, file.path(OUT, "loso_skill.csv"))

# Paired block differences (same bootstrap sites).
DIFFS <- tribble(
  ~frame, ~a, ~b, ~label,
  "perception", "auditory_composition", "auditory", "composition added to auditory",
  "perception", "auditory", "composition", "auditory minus composition",
  "acoustic", "measured_composition", "measured", "composition added to measured acoustics",
  "acoustic", "measured", "composition", "measured acoustics minus composition",
  "acoustic", "all", "auditory_measured", "composition added to auditory + measured")
diff_tab <- pmap_dfr(DIFFS, function(frame, a, b, label) {
  map_dfr(c(TRUE, FALSE), function(w) {
    pa <- point %>% filter(.data$frame == !!frame, wave_adjusted == w, block == a)
    pb <- point %>% filter(.data$frame == !!frame, wave_adjusted == w, block == b)
    ba <- boot %>% filter(.data$frame == !!frame, wave_adjusted == w, block == a)
    bb <- boot %>% filter(.data$frame == !!frame, wave_adjusted == w, block == b)
    map_dfr(names(OUTS), function(y) {
      dd <- ba %>% filter(outcome == y) %>% select(boot, r2_person, r2_site) %>%
        left_join(bb %>% filter(outcome == y) %>% select(boot, r2_person, r2_site),
                  by = "boot", suffix = c("_a", "_b"))
      tibble(frame = frame, wave_adjusted = w, contrast = label, outcome = y,
             outcome_label = OUTS[y],
             d_person = pa$r2_person[pa$outcome == y] - pb$r2_person[pb$outcome == y],
             d_person_lo = quantile(dd$r2_person_a - dd$r2_person_b, 0.025),
             d_person_hi = quantile(dd$r2_person_a - dd$r2_person_b, 0.975),
             d_site = pa$r2_site[pa$outcome == y] - pb$r2_site[pb$outcome == y],
             d_site_lo = quantile(dd$r2_site_a - dd$r2_site_b, 0.025, na.rm = TRUE),
             d_site_hi = quantile(dd$r2_site_a - dd$r2_site_b, 0.975, na.rm = TRUE))
    })
  })
})
write_csv(diff_tab, file.path(OUT, "loso_block_differences.csv"))

# Site means for the display (perception frame, wave-adjusted, composition and
# auditory blocks), centred on the wave-only baseline.
site_means <- ps %>% filter(wave_adjusted, (frame == "perception" &
                                              block %in% c("composition", "auditory")) |
                              (frame == "acoustic" & block == "measured")) %>%
  mutate(obs_c = obs_mean - base_mean, pred_c = pred_mean - base_mean) %>%
  left_join(distinct(base, site = as.character(site), site_label, space_class), by = "site")
write_csv(site_means, file.path(OUT, "loso_site_means.csv"))

# ---- Part 2: proxy check - variance share of the composition block -----------
r2 <- function(f, dat) summary(lm(f, dat))$r.squared
share_rows <- function(dat, y) {
  w <- "wave"; s <- c("site", "wave")
  tibble(
    pooled = r2(reformulate(c(w, VIS), y), dat) - r2(reformulate(w, y), dat),
    pooled_net_auditory = r2(reformulate(c(w, SOUND, VIS), y), dat) -
      r2(reformulate(c(w, SOUND), y), dat),
    pooled_net_measured = if (all(MEAS %in% names(dat)))
      r2(reformulate(c(w, MEAS, VIS), y), dat) - r2(reformulate(c(w, MEAS), y), dat) else NA,
    within = r2(reformulate(c(s, VIS), y), dat) - r2(reformulate(s, y), dat),
    within_net_auditory = r2(reformulate(c(s, SOUND, VIS), y), dat) -
      r2(reformulate(c(s, SOUND), y), dat),
    auditory_pooled = r2(reformulate(c(w, SOUND), y), dat) - r2(reformulate(w, y), dat),
    measured_pooled = if (all(MEAS %in% names(dat)))
      r2(reformulate(c(w, MEAS), y), dat) - r2(reformulate(w, y), dat) else NA)
}
shares <- map_dfr(names(OUTS), ~ share_rows(d_ac, .x) %>% mutate(outcome = .x)) %>%
  mutate(frame = "acoustic")
set.seed(SEED + 1)
boot_sh <- map_dfr(seq_len(N_BOOT), function(b) {
  ss <- boot_idx[[b]]
  db <- map_dfr(seq_along(ss), ~ d_ac %>% filter(site == ss[.x]) %>%
                  mutate(site = factor(paste0("b", .x))))
  db$wave <- droplevels(db$wave)
  map_dfr(names(OUTS), ~ share_rows(db, .x) %>% mutate(outcome = .x, boot = b))
})
sh_long <- shares %>% pivot_longer(-c(outcome, frame), names_to = "quantity", values_to = "share")
sh_ci <- boot_sh %>% pivot_longer(-c(outcome, boot), names_to = "quantity", values_to = "share") %>%
  group_by(outcome, quantity) %>%
  summarise(lo = quantile(share, 0.025, na.rm = TRUE), hi = quantile(share, 0.975, na.rm = TRUE),
            .groups = "drop")
proxy <- sh_long %>% left_join(sh_ci, by = c("outcome", "quantity")) %>%
  mutate(outcome_label = OUTS[outcome], n = nrow(d_ac))
write_csv(proxy, file.path(OUT, "proxy_shares.csv"))

# ---- Part 3: between-site share of each outcome ------------------------------
icc <- map_dfr(names(OUTS), function(y) {
  fit <- lmer(reformulate(c("wave", "(1 | site)"), y), data = d_per, REML = TRUE)
  vc <- as.data.frame(VarCorr(fit))
  tibble(outcome = y, outcome_label = OUTS[y],
         icc_site = vc$vcov[vc$grp == "site"] / sum(vc$vcov))
})
write_csv(icc, file.path(OUT, "between_site_share.csv"))

# ---- console summary ------------------------------------------------------------
options(width = 200)
cat("\n== LOSO skill (wave-adjusted) ==\n")
print(skill %>% filter(wave_adjusted) %>%
        select(frame, block, outcome_label, r2_person, r2_person_lo, r2_person_hi,
               r2_site, r2_site_lo, r2_site_hi, rho_site, r2_within) %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))), n = 60)
cat("\n== Block differences (wave-adjusted) ==\n")
print(diff_tab %>% filter(wave_adjusted) %>% mutate(across(where(is.numeric), ~ round(.x, 3))),
      n = 40, width = Inf)
cat("\n== Proxy shares ==\n")
print(proxy %>% mutate(across(c(share, lo, hi), ~ round(.x, 4))), n = 40)
cat("\n== Between-site share ==\n"); print(icc)
