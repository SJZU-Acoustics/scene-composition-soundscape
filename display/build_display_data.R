# =============================================================================
# build_display_data.R -- Build the display-item data locks for the paper.
# Reads the analysis outcome CSVs plus the survey table, and writes:
#   output/fig*_*.csv        one data-lock CSV per figure panel
#   output/table1_data.csv   site table data
#   tables/table1_sample.tex LaTeX fragment (Table 1)
#   tables/table2_effects.tex LaTeX fragment (Table 2)
# No model is refitted here; every estimate is copied from the vetted
# exploration artefacts.
# =============================================================================

source("display/style.R")

# ---- Clean core: photo-valid analytic view ----------------------------------
core <- readr::read_csv(CLEAN_CSV, show_col_types = FALSE, guess_max = 3000) %>%
  filter(include_core == 1) %>%
  mutate(photo_valid = as.integer(photo_ok == 1 & photo_present == 1),
         season = word(wave_label, 2))

pv <- filter(core, photo_valid == 1)

stopifnot(nrow(pv) == 2262)                     # photo-valid core (A08/A05)
stopifnot(sum(!is.na(pv$laeq_a)) == 2229)       # available-channel acoustic rows

# ---- Table 1: sites ----------------------------------------------------------
season_abb <- c(spring = "Sp", summer = "Su", autumn = "Au")
tab1 <- pv %>%
  mutate(site_en = SITE_EN[site_key],
         class_en = recode(space_class, park = "Park",
                           civic_square = "Civic square",
                           residential = "Residential"),
         season_rank = match(word(wave_label, 2), c("spring", "summer", "autumn")),
         wave_ord = as.integer(word(wave_label, 1)) * 10 + season_rank,
         wave_abb = paste0(word(wave_label, 1), "~", season_abb[word(wave_label, 2)])) %>%
  group_by(site_key, site_en, class_en) %>%
  summarise(waves = paste(unique(wave_abb[order(wave_ord)]), collapse = ", "),
            n_photo = n(),
            n_ac = sum(!is.na(laeq_a)),
            greenery = mean(comp_greenery, na.rm = TRUE),
            sky = mean(comp_sky, na.rm = TRUE),
            oa = mean(overall_soundscape, na.rm = TRUE),
            laeq = mean(laeq_a, na.rm = TRUE), .groups = "drop") %>%
  arrange(factor(class_en, levels = c("Park", "Civic square", "Residential")),
          desc(n_photo))
write_csv(tab1, file.path(OUT_DIR, "table1_data.csv"))

fmt2 <- function(x) formatC(x, digits = 2, format = "f")
fmt1 <- function(x) formatC(x, digits = 1, format = "f")
tab1_tex <- c(
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}>{\\raggedright\\arraybackslash}p{2.2cm}>{\\raggedright\\arraybackslash}p{1.4cm}>{\\raggedright\\arraybackslash}p{2.7cm}rrrrrr@{}}",
  "\\toprule",
  paste0("Site & Space class & Survey waves & $n$ & $n_{\\mathrm{ac}}$ & ",
         "Greenery & Sky & OA & $L_\\mathrm{Aeq}$ \\\\"),
  "\\midrule",
  with(tab1, sprintf("%s & %s & %s & %d & %d & %s & %s & %s & %s \\\\",
                     site_en, class_en, waves, n_photo, n_ac,
                     fmt2(greenery), fmt2(sky), fmt2(oa), fmt1(laeq))),
  "\\midrule",
  sprintf("Total &  &  & %s & %s &  &  &  &  \\\\",
          format(sum(tab1$n_photo), big.mark = ","),
          format(sum(tab1$n_ac), big.mark = ",")),
  "\\bottomrule",
  "\\end{tabular*}")
writeLines(tab1_tex, file.path(TAB_DIR, "table1_sample.tex"))

# ---- Fig 1a: in-sample variance shares (A08) --------------------------------
fig1a <- read_csv(odir("analysis_08_variance_partitioning",
                       "variance_partitioning.csv"), show_col_types = FALSE) %>%
  filter(outcome != "raq_mean") %>%
  select(outcome, n, sound_total, visual_total, visual_to_sound_ratio)
write_csv(fig1a, file.path(OUT_DIR, "fig1_variance.csv"))

# ---- Fig 1b: out-of-sample increments and permutation nulls (A28) -----------
null_draws <- read_csv(odir("analysis_28_cv_readout",
                            "permutation_null_increments.csv"),
                       show_col_types = FALSE) %>%
  group_by(outcome) %>%
  summarise(null_lo05 = quantile(inc_lm, 0.05),
            null_hi95 = quantile(inc_lm, 0.95), .groups = "drop")
fig1b <- read_csv(odir("analysis_28_cv_readout", "permutation_increments.csv"),
                  show_col_types = FALSE) %>%
  left_join(null_draws, by = "outcome") %>%
  transmute(outcome = sub("_z$", "", outcome), outcome_label,
            inc = inc_lm, p_perm = p_inc_lm, null_lo05, null_hi95)
write_csv(fig1b, file.path(OUT_DIR, "fig1_cv.csv"))

# ---- Fig 1c: prediction at held-out sites (A34, leave-one-site-out) ---------
fig1c <- read_csv(odir("analysis_34_site_heldout_mapping", "loso_skill.csv"),
                  show_col_types = FALSE) %>%
  filter(frame == "perception", wave_adjusted,
         block %in% c("auditory", "composition")) %>%
  transmute(outcome = sub("_z$", "", outcome), outcome_label, block,
            r2_person, r2_person_lo, r2_person_hi)
stopifnot(nrow(fig1c) == 8)
write_csv(fig1c, file.path(OUT_DIR, "fig1_heldout.csv"))

# ---- Fig 2a: greenery share by season (A17) ---------------------------------
fig2a <- read_csv(odir("analysis_17_season_phenology", "quantity_by_season.csv"),
                  show_col_types = FALSE) %>%
  select(season, n, comp_greenery)
write_csv(fig2a, file.path(OUT_DIR, "fig2_quantity.csv"))

# ---- Fig 2b: greenery -> OA slope by season (A17) ---------------------------
fig2b <- read_csv(odir("analysis_17_season_phenology", "stratified_slopes.csv"),
                  show_col_types = FALSE) %>%
  filter(effect == "greenery_oa") %>%
  select(season, n, estimate, std.error, p.value)
write_csv(fig2b, file.path(OUT_DIR, "fig2_greenery_slopes.csv"))

# ---- Fig 2c: green pixels per unit segmented greenery (A31 subset frame; ----
# means and the A25 bootstrap contrast both on the multi-season subset)
fig2c <- read_csv(odir("analysis_31_panorama_clustering",
                       "quality_gap_by_season_subset.csv"),
                  show_col_types = FALSE)
boot <- read_csv(odir("analysis_25_chromatic_phenology",
                      "bootstrap_ratio_contrast.csv"), show_col_types = FALSE)
fig2c <- bind_cols(fig2c,
                   contrast_est = boot$est, contrast_mean = boot$mean,
                   contrast_lo = boot$lo, contrast_hi = boot$hi)
write_csv(fig2c, file.path(OUT_DIR, "fig2_quality_gap.csv"))

# ---- Fig 2d: source slopes on OA by season (A20) ----------------------------
fig2d <- read_csv(odir("analysis_20_auditory_phenology",
                       "stratified_source_slopes.csv"), show_col_types = FALSE) %>%
  filter(outcome == "overall_soundscape_z",
         predictor %in% c("natural_sounds_z", "traffic_noise_z")) %>%
  select(season, n, predictor, estimate, std.error, p.value)
write_csv(fig2d, file.path(OUT_DIR, "fig2_source_slopes.csv"))

# ---- Fig 3a: composites -> ISO-E screen (A05) -------------------------------
fig3a <- read_csv(odir("analysis_05_visual_to_evaluation",
                       "visual_evaluation_screen.csv"), show_col_types = FALSE) %>%
  filter(outcome == "isoe") %>%
  select(composite, estimate, std.error, p.value, q_fdr)
write_csv(fig3a, file.path(OUT_DIR, "fig3_isoe_screen.csv"))

# ---- Fig 3b: match/conflict indices -> ISO-E (A15) --------------------------
screen15 <- read_csv(odir("analysis_15_congruence_conflict",
                          "match_conflict_screen.csv"), show_col_types = FALSE)
boot15 <- read_csv(odir("analysis_15_congruence_conflict",
                        "cluster_bootstrap.csv"), show_col_types = FALSE)
fig3b <- screen15 %>%
  filter(outcome == "isoe") %>%
  select(index, label, estimate, std.error, p.value, q_fdr) %>%
  left_join(filter(boot15, outcome == "isoe") %>%
              select(index, ci_lo, ci_hi, excludes_zero), by = "index")
write_csv(fig3b, file.path(OUT_DIR, "fig3_indices.csv"))

# ---- Fig 4a: water -> OA by LAeq tertile (A19) ------------------------------
fig4a <- read_csv(odir("analysis_19_loudness_conditioned_vision",
                       "tertile_stratified_slopes.csv"), show_col_types = FALSE) %>%
  filter(effect == "water_oa") %>%
  select(laeq_tertile, n, estimate, std.error, p.value)
pl <- read_csv(odir("analysis_19_loudness_conditioned_vision",
                    "perceived_loudness_moderation.csv"), show_col_types = FALSE) %>%
  filter(effect == "water_oa")
fig4a <- fig4a %>%
  mutate(pl_interaction = pl$b_interaction, pl_p = pl$p.value, pl_q = pl$q.value)
write_csv(fig4a, file.path(OUT_DIR, "fig4_water_tertiles.csv"))

# ---- Fig 4b: blueness profile (A26/A29) -------------------------------------
scan26 <- read_csv(odir("analysis_26_quality_outcome_scan",
                        "quality_outcome_scan.csv"), show_col_types = FALSE) %>%
  filter(var == "sky_blue_ratio", outcome == "overall_soundscape_z")
adj29 <- read_csv(odir("analysis_29_blueness_integration",
                       "blueness_sound_adjusted.csv"), show_col_types = FALSE) %>%
  filter(outcome == "overall_soundscape_z")
sea29 <- read_csv(odir("analysis_29_blueness_integration",
                       "blueness_season_stratified_oa.csv"), show_col_types = FALSE)
cmp29 <- read_csv(odir("analysis_29_blueness_integration",
                       "blueness_loudness_compensation.csv"), show_col_types = FALSE) %>%
  filter(outcome == "overall_soundscape_z")
fig4b <- tibble(
  row = c("raw", "sound_adjusted", "spring", "autumn", "x_loudness"),
  label = c("Overall effect", "Sound-adjusted", "Spring subset", "Autumn subset",
            "× perceived loudness"),
  estimate = c(scan26$est, adj29$est, sea29$est[sea29$season == "spring"],
               sea29$est[sea29$season == "autumn"], cmp29$b_interaction),
  std.error = c(scan26$se, adj29$se, sea29$se[sea29$season == "spring"],
                sea29$se[sea29$season == "autumn"], cmp29$se_interaction),
  p.value = c(scan26$p, adj29$p, sea29$p[sea29$season == "spring"],
              sea29$p[sea29$season == "autumn"], cmp29$p_interaction),
  n = c(scan26$n, adj29$n, sea29$n[sea29$season == "spring"],
        sea29$n[sea29$season == "autumn"], cmp29$n),
  q.value = c(scan26$q, NA, NA, NA, cmp29$q_interaction))
write_csv(fig4b, file.path(OUT_DIR, "fig4_blueness.csv"))

# ---- Fig 4c: heard-and-seen traffic -> OA (A06/A17) -------------------------
cong <- read_csv(odir("analysis_06_congruence", "congruence_interactions.csv"),
                 show_col_types = FALSE) %>%
  filter(outcome == "overall_soundscape", sound == "traffic_noise") %>%
  transmute(stratum = "pooled", estimate, std.error, p.value, q_fdr)
cong_season <- read_csv(odir("analysis_17_season_phenology",
                             "stratified_slopes.csv"), show_col_types = FALSE) %>%
  filter(effect == "traffic_congruence_oa", season != "summer") %>%
  transmute(stratum = season, estimate, std.error, p.value, q_fdr = NA)
fig4c <- bind_rows(cong, cong_season)
write_csv(fig4c, file.path(OUT_DIR, "fig4_congruence.csv"))

# ---- Table 2: FDR-surviving effects -----------------------------------------
scr05 <- read_csv(odir("analysis_05_visual_to_evaluation",
                       "visual_evaluation_screen.csv"), show_col_types = FALSE)
adj05 <- read_csv(odir("analysis_05_visual_to_evaluation",
                       "survivors_laeq_adjusted.csv"), show_col_types = FALSE)
el11 <- read_csv(odir("analysis_11_element_candidates", "element_scan.csv"),
                 show_col_types = FALSE)
bal15 <- screen15 %>% filter(outcome == "isoe", index == "idx_balance_match")
adj33 <- read_csv(odir("analysis_33_water_sound_adjusted", "water_sound_adjusted.csv"),
                  show_col_types = FALSE) %>%
  filter(outcome == "overall_soundscape_z", specification == "all auditory measures adjusted")
mod33 <- read_csv(odir("analysis_33_water_sound_adjusted", "water_loudness_sound_adjusted.csv"),
                  show_col_types = FALSE) %>%
  filter(specification == "source saliencies adjusted")
stopifnot(nrow(adj33) == 1, nrow(mod33) == 1)

t2 <- tribble(
  ~group, ~effect, ~est, ~se, ~p, ~q, ~family, ~n,
  "Scene composition (per SD area share)",
    "Visible vehicles $\\rightarrow$ \\textit{ISO-E}",
    scr05$estimate[scr05$outcome == "isoe" & scr05$composite == "comp_vehicles"],
    scr05$std.error[scr05$outcome == "isoe" & scr05$composite == "comp_vehicles"],
    scr05$p.value[scr05$outcome == "isoe" & scr05$composite == "comp_vehicles"],
    scr05$q_fdr[scr05$outcome == "isoe" & scr05$composite == "comp_vehicles"], 40, 2262,
  "Scene composition (per SD area share)",
    "\\quad the same, $L_\\mathrm{Aeq}$-adjusted",
    adj05$estimate[1], adj05$std.error[1], adj05$p.value[1], NA, NA, 2229,
  "Scene condition (per SD condition metric)",
    "Sky blueness $\\rightarrow$ OA",
    scan26$est, scan26$se, scan26$p, scan26$q, 36, 2261,
  "Scene condition (per SD condition metric)",
    "\\quad the same, sound-adjusted",
    adj29$est, adj29$se, adj29$p, NA, NA, 2261,
  "Visual elements (present vs absent)",
    "Water in view $\\rightarrow$ OA",
    el11$estimate[el11$element == "water_any" & el11$outcome == "overall_soundscape"],
    el11$std.error[el11$element == "water_any" & el11$outcome == "overall_soundscape"],
    el11$p.value[el11$element == "water_any" & el11$outcome == "overall_soundscape"],
    el11$q_fdr[el11$element == "water_any" & el11$outcome == "overall_soundscape"], 8, 2262,
  "Visual elements (present vs absent)",
    "\\quad the same, sound-adjusted",
    adj33$estimate, adj33$se, adj33$p.value, NA, NA, 2262,
  "Visual elements (present vs absent)",
    "Water in view $\\rightarrow$ \\textit{ISO-P}",
    el11$estimate[el11$element == "water_any" & el11$outcome == "isop"],
    el11$std.error[el11$element == "water_any" & el11$outcome == "isop"],
    el11$p.value[el11$element == "water_any" & el11$outcome == "isop"],
    el11$q_fdr[el11$element == "water_any" & el11$outcome == "isop"], 8, 2262,
  "Visual elements (present vs absent)",
    "Flowers in view $\\rightarrow$ OA",
    el11$estimate[el11$element == "flower_any" & el11$outcome == "overall_soundscape"],
    el11$std.error[el11$element == "flower_any" & el11$outcome == "overall_soundscape"],
    el11$p.value[el11$element == "flower_any" & el11$outcome == "overall_soundscape"],
    el11$q_fdr[el11$element == "flower_any" & el11$outcome == "overall_soundscape"], 8, 2262,
  "Perceived sound characteristics: perceived loudness",
    "Water $\\times$ perceived loudness $\\rightarrow$ OA",
    pl$b_interaction, pl$se, pl$p.value, pl$q.value, 7, 2229,
  "Perceived sound characteristics: perceived loudness",
    "\\quad the same, source-adjusted",
    mod33$estimate, mod33$se, mod33$p.value, NA, NA, 2229,
  "Perceived sound characteristics: sound-source combination",
    "Heard-and-seen traffic $\\rightarrow$ OA",
    cong$estimate[1], cong$std.error[1], cong$p.value[1], cong$q_fdr[1], 9, 2262,
  "Perceived sound characteristics: sound-source combination",
    "Nature-hard balance match $\\rightarrow$ \\textit{ISO-E}",
    bal15$estimate[1], bal15$std.error[1], bal15$p.value[1], bal15$q_fdr[1],
    nrow(filter(screen15, !is.na(q_fdr))), 2262)
write_csv(t2, file.path(OUT_DIR, "table2_data.csv"))

fmt3 <- function(x) ifelse(is.na(x), "",
                    ifelse(x < 0, paste0("$-$", formatC(abs(x), digits = 3, format = "f")),
                           formatC(x, digits = 3, format = "f")))
fmt_pq <- function(x) ifelse(is.na(x), "--",
                      ifelse(x < 0.001, "$<$0.001", formatC(x, digits = 3, format = "f")))
t2_rows <- t2 %>%
  mutate(row = sprintf("%s & %s & %s & %s & %s \\\\", effect,
                       fmt3(est), fmt3(se), fmt_pq(p), fmt_pq(q)))
t2_tex <- c(
  "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrr@{}}",
  "\\toprule",
  "Effect & Estimate & SE & $p$ & $q$ \\\\",
  "\\midrule")
for (g in unique(t2_rows$group)) {
  t2_tex <- c(t2_tex,
              sprintf("\\multicolumn{5}{@{}l}{\\emph{%s}} \\\\", g),
              t2_rows$row[t2_rows$group == g])
}
t2_tex <- c(t2_tex, "\\bottomrule", "\\end{tabular*}")
writeLines(t2_tex, file.path(TAB_DIR, "table2_effects.tex"))

message("Display data locks written to ", OUT_DIR)
