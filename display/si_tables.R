# =============================================================================
# si_tables.R -- Generate the supplementary table BODIES as bare tabular
# fragments in tables/si/, embedding the FULL contents of the analysis outcome
# CSVs (every table in full, no truncation). Captions, floats and numbering
# live in the Supplementary Information document that \input{}s them.
# =============================================================================

source("display/style.R")
SI_DIR <- file.path(P22_ROOT, "tables", "si")
dir.create(SI_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_LAB <- c(overall_soundscape = "OA", overall_soundscape_z = "OA",
             isop = "ISO-P", isop_z = "ISO-P",
             isoe = "ISO-E", isoe_z = "ISO-E",
             harmony = "Harmony", harmony_z = "Harmony",
             perceived_loudness = "Perceived loudness",
             perceived_loudness_z = "Perceived loudness",
             raq_mean = "RAQ", raq_mean_z = "RAQ")
lab_out <- function(x) ifelse(x %in% names(OUT_LAB), OUT_LAB[x], x)
lab_comp <- function(x) ifelse(x %in% names(COMP_LABELS), COMP_LABELS[x],
                        gsub("_", " ", x))

ARROW <- " $\\rightarrow$ "
EFFECT_LAB <- c(
  greenery_isop = paste0("Greenery", ARROW, "ISO-P"),
  greenery_oa = paste0("Greenery", ARROW, "OA"),
  greenery_harmony = paste0("Greenery", ARROW, "Harmony"),
  sky_isop = paste0("Sky", ARROW, "ISO-P"),
  people_oa = paste0("People", ARROW, "OA"),
  vehicles_isoe = paste0("Visible vehicles", ARROW, "ISO-E"),
  vehicles_to_eventfulness = paste0("Visible vehicles", ARROW, "ISO-E"),
  traffic_congruence_oa = paste0("Heard-and-seen traffic", ARROW, "OA"),
  traffic_seen_heard_to_oa = paste0("Heard-and-seen traffic", ARROW, "OA"),
  balance_match_isoe = paste0("Nature-hard balance match", ARROW, "ISO-E"),
  water_oa = paste0("Water in view", ARROW, "OA"),
  water_to_oa = paste0("Water in view", ARROW, "OA"),
  water_isop = paste0("Water in view", ARROW, "ISO-P"),
  flower_oa = paste0("Flowers in view", ARROW, "OA"),
  flowers_oa = paste0("Flowers in view", ARROW, "OA"),
  flower_to_oa = paste0("Flowers in view", ARROW, "OA"),
  blueness_oa = paste0("Sky blueness", ARROW, "OA"),
  water_x_pl_oa = paste0("Water $\\times$ perceived loudness", ARROW, "OA"),
  greenery_x_season_oa = paste0("Greenery $\\times$ season", ARROW, "OA"),
  natural_x_season_oa = paste0("Natural sounds $\\times$ season", ARROW, "OA"),
  natural_x_season_isop = paste0("Natural sounds $\\times$ season", ARROW, "ISO-P"),
  natural_x_season_harmony = paste0("Natural sounds $\\times$ season", ARROW, "Harmony"))
lab_effect <- function(x) {
  stopifnot(all(x %in% names(EFFECT_LAB)))
  unname(EFFECT_LAB[x])
}
VAR_LAB <- c(
  ctx_diversity = "Scene diversity", ctx_dominance = "Dominance",
  ctx_natural_hard_mixedness = "Natural--hard mixedness", ctx_open_horizon = "Open horizon",
  ctx_mobility_infrastructure = "Mobility infrastructure", ctx_social_amenity = "Social amenity",
  ctx_restorative_cues = "Restorative cues",
  water_any = "Water in view", flower_any = "Flowers in view",
  sculpture_any = "Sculpture or fountain in view", person_share_z = "People share",
  veg_shannon = "Vegetation diversity", flower_in_green = "Flower content of greenery",
  tree_frac = "Tree fraction", grass_frac = "Grass fraction", herb_frac = "Herbaceous fraction",
  layered = "Layering", visual = "Visual", cross_modal = "Cross-modal")
lab_var <- function(x) {
  stopifnot(all(x %in% names(VAR_LAB)))
  unname(VAR_LAB[x])
}
WX_LAB <- c(temp_z = "Temperature", rhum_z = "Relative humidity", wspd_z = "Wind speed",
            rain_day = "Rain day", aqi_z = "Air quality index")
nfmt <- function(x) format(as.integer(x), big.mark = ",", trim = TRUE)
# one label per photograph metric across Supplementary Tables S18 and S19
METRIC_LAB <- c("Colourfulness (HS)" = "Colourfulness",
                "Vegetation hue yellow-shift" = "Vegetation hue (yellow shift)")

esc <- function(x) gsub("_", "\\\\_", as.character(x))
num3 <- function(x) ifelse(is.na(suppressWarnings(as.numeric(x))), esc(x),
                    ifelse(abs(as.numeric(x)) < 0.0005 & as.numeric(x) != 0,
                           formatC(as.numeric(x), format = "e", digits = 1),
                           formatC(as.numeric(x), digits = 3, format = "f")))
num4 <- function(x) formatC(as.numeric(x), digits = 4, format = "f")
pfmt <- function(x) ifelse(is.na(x), "--",
                    ifelse(x < 0.001, "$<$0.001",
                           formatC(x, digits = 3, format = "f")))

si_table <- function(df, file, align = NULL) {
  align <- align %||% paste0("l", strrep("r", ncol(df) - 1))
  hdr <- paste(colnames(df), collapse = " & ")
  rows <- apply(df, 1, function(r) paste(r, collapse = " & "))
  # negative numbers take a true minus sign, as in the main-text tables
  rows <- gsub("(?<![[:alnum:]$-])-(?=[0-9.])", "$-$", rows, perl = TRUE)
  rows <- gsub("([0-9.]+)e-0*([0-9]+)", "\\1$\\\\times10^{-\\2}$", rows)
  # ISO-P and ISO-E are set in italics wherever they print (registry form)
  ital_iso <- function(x) gsub("(?<!\\\\textit\\{)ISO-([PE])(?![}[:alnum:]])",
                               "\\\\textit{ISO-\\1}", x, perl = TRUE)
  hdr <- ital_iso(hdr)
  rows <- ital_iso(rows)
  lines <- c(sprintf("\\begin{tabular}{%s}", align),
             "\\toprule", paste0(hdr, " \\\\"), "\\midrule",
             paste0(rows, " \\\\"),
             "\\bottomrule", "\\end{tabular}")
  writeLines(lines, file.path(SI_DIR, file))
  message("  wrote ", file, " (", nrow(df), " rows)")
}

# Values that would print as the 0.050 threshold are shown to four places, or
# to five where four would still round onto it (e.g. 50/1001 = 0.04995).
thr_fmt <- function(x) ifelse(formatC(x, digits = 4, format = "f") == "0.0500" & x != 0.05,
                              formatC(x, digits = 5, format = "f"),
                              formatC(x, digits = 4, format = "f"))
qfmt <- function(x) ifelse(!is.na(x) & round(x, 3) == 0.05, thr_fmt(x), pfmt(x))

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- S1: composition screen (A05, 40 tests) ---------------------------------
s1 <- read_csv(odir("analysis_05_visual_to_evaluation",
                    "visual_evaluation_screen.csv"), show_col_types = FALSE) %>%
  transmute(Outcome = lab_out(outcome), Composite = lab_comp(composite),
            Estimate = num3(estimate), SE = num3(std.error),
            `$p$` = pfmt(p.value), `$q$` = pfmt(q_fdr))
si_table(s1, "si_composition_screen.tex")

# ---- S2: out-of-sample CV skill and permutation calibration (A28) -----------
OUTCOME_ORDER <- c("OA", "ISO-P", "Harmony", "ISO-E")
COND_MODELS <- c(auditory = "Linear, auditory",
                 auditory_composition = "Linear, auditory + composition",
                 auditory_condition = "Linear, auditory + condition",
                 auditory_composition_condition = "Linear, auditory + composition + condition")
s2a <- bind_rows(
  read_csv(odir("analysis_28_cv_readout", "cv_skill.csv"), show_col_types = FALSE) %>%
    transmute(frame = 1, Outcome = outcome_label,
              Model = recode(model, lm_sound = "Linear, auditory",
                             lm_visual = "Linear, composition",
                             lm_both = "Linear, auditory + composition",
                             rf_sound = "Forest, auditory",
                             rf_visual = "Forest, composition",
                             rf_both = "Forest, auditory + composition"),
              order = row_number(), r2_mean, r2_sd),
  read_csv(odir("analysis_28_cv_readout", "condition_block_skill.csv"),
           show_col_types = FALSE) %>%
    transmute(frame = 2, Outcome = outcome_label, Model = unname(COND_MODELS[model]),
              order = match(model, names(COND_MODELS)), r2_mean, r2_sd)) %>%
  arrange(frame, match(Outcome, OUTCOME_ORDER), order) %>%
  transmute(`$n$` = ifelse(frame == 1, "2,262", "2,261"), Outcome, Model,
            `CV $R^2$` = num4(r2_mean), SD = num4(r2_sd))
si_table(s2a, "si_cv_skill.tex", align = "lllrr")
perm <- read_csv(odir("analysis_28_cv_readout", "permutation_increments.csv"),
                 show_col_types = FALSE)
null_q <- read_csv(odir("analysis_28_cv_readout", "permutation_null_increments.csv"),
                   show_col_types = FALSE) %>% group_by(outcome) %>%
  summarise(lm_lo = quantile(inc_lm, 0.05, na.rm = TRUE), lm_hi = quantile(inc_lm, 0.95, na.rm = TRUE),
            rf_lo = quantile(inc_rf, 0.05, na.rm = TRUE), rf_hi = quantile(inc_rf, 0.95, na.rm = TRUE),
            .groups = "drop")
perm <- perm %>% left_join(null_q, by = "outcome")
INC_BLOCK <- c(composition_over_auditory = "Composition", condition_over_auditory = "Condition",
               condition_over_both = "Condition")
INC_BASE <- c(composition_over_auditory = "Auditory", condition_over_auditory = "Auditory",
              condition_over_both = "Auditory + composition")
s2b <- bind_rows(
  perm %>% transmute(frame = 1, Outcome = outcome_label, order = 1,
                     Added = "Composition", Base = "Auditory", Estimator = "Linear",
                     inc = inc_lm, lo = lm_lo, hi = lm_hi, p = p_inc_lm),
  perm %>% transmute(frame = 1, Outcome = outcome_label, order = 2,
                     Added = "Composition", Base = "Auditory", Estimator = "Forest",
                     inc = inc_rf, lo = rf_lo, hi = rf_hi, p = p_inc_rf),
  read_csv(odir("analysis_28_cv_readout", "condition_block_increments.csv"),
           show_col_types = FALSE) %>%
    transmute(frame = 2, Outcome = outcome_label, order = match(increment, names(INC_BLOCK)),
              Added = unname(INC_BLOCK[increment]), Base = unname(INC_BASE[increment]),
              Estimator = "Linear",
              inc, lo = null_lo05, hi = null_hi95, p = p_perm)) %>%
  arrange(frame, match(Outcome, OUTCOME_ORDER), order) %>%
  transmute(`$n$` = ifelse(frame == 1, "2,262", "2,261"), Outcome, `Block added` = Added,
            `Added to` = Base, Estimator, Gain = num4(inc),
            `Permuted, 5th--95th` = sprintf("[%s, %s]", num4(lo), num4(hi)),
            # a permutation p cannot fall below 1 / (permutations + 1), so it is
            # printed as computed and never as an inequality
            `$p_{\\mathrm{perm}}$` = ifelse(round(p, 3) == 0.05, thr_fmt(p),
                                             formatC(p, digits = 3, format = "f")))
si_table(s2b, "si_cv_permutation.tex", align = "lllllrrr")

# ---- S3: curvature check (A30) ----------------------------------------------
s3 <- read_csv(odir("analysis_30_curvature_check", "curvature_scan.csv"),
               show_col_types = FALSE) %>%
  transmute(Outcome = outcome_label, Composite = lab_comp(composite),
            `Quadratic estimate` = num3(b_quad), SE = num3(se_quad),
            `$p$` = pfmt(p_quad), `$q$` = pfmt(q_quad))
si_table(s3, "si_curvature.tex")

# ---- Prediction at held-out sites (A34, leave-one-site-out) -----------------
loso <- read_csv(odir("analysis_34_site_heldout_mapping", "loso_skill.csv"),
                 show_col_types = FALSE)
LOSO_BLOCKS <- tribble(
  ~frame, ~block, ~Block, ~n,
  "perception", "composition", "Composition", "2,262",
  "perception", "auditory", "Auditory", "2,262",
  "perception", "auditory_composition", "Auditory + composition", "2,262",
  "acoustic", "measured", "Measured acoustics", "2,228")
# R2 values and gains print to three decimals throughout (no exponent form,
# no negative zero)
r3 <- function(x) { x <- round(as.numeric(x), 3); x[x == 0] <- 0
  formatC(x, digits = 3, format = "f") }
ci3 <- function(lo, hi) sprintf("[%s, %s]", r3(lo), r3(hi))
s_loso <- loso %>% filter(wave_adjusted) %>%
  inner_join(LOSO_BLOCKS, by = c("frame", "block")) %>%
  arrange(match(outcome_label, OUTCOME_ORDER), match(block, LOSO_BLOCKS$block)) %>%
  transmute(`$n$` = n, Outcome = outcome_label, Block,
            `Person-level $R^2$` = r3(r2_person), `95\\% CI` = ci3(r2_person_lo, r2_person_hi),
            `Site-level $R^2$` = r3(r2_site), `95\\% CI ` = ci3(r2_site_lo, r2_site_hi))
stopifnot(nrow(s_loso) == 16)
si_table(s_loso, "si_site_heldout.tex", align = "lllrrrr")
gain <- read_csv(odir("analysis_34_site_heldout_mapping", "loso_block_differences.csv"),
                 show_col_types = FALSE) %>%
  filter(frame == "perception", contrast == "composition added to auditory")
icc <- read_csv(odir("analysis_34_site_heldout_mapping", "between_site_share.csv"),
                show_col_types = FALSE)
s_gain <- gain %>% filter(wave_adjusted) %>%
  left_join(gain %>% filter(!wave_adjusted) %>%
              select(outcome, nw = d_person, nw_lo = d_person_lo, nw_hi = d_person_hi),
            by = "outcome") %>%
  left_join(select(icc, outcome, icc_site), by = "outcome") %>%
  arrange(match(outcome_label, OUTCOME_ORDER)) %>%
  transmute(Outcome = outcome_label,
            `Between sites (\\%)` = formatC(100 * icc_site, digits = 1, format = "f"),
            `Person level` = paste(r3(d_person), ci3(d_person_lo, d_person_hi)),
            `Site level` = paste(r3(d_site), ci3(d_site_lo, d_site_hi)),
            `Person level, no wave terms` = paste(r3(nw), ci3(nw_lo, nw_hi)))
stopifnot(nrow(s_gain) == 4)
si_table(s_gain, "si_site_heldout_gain.tex", align = "lrrrr")

# ---- S4: match/conflict full screen (A15, 35 tests) -------------------------
s4 <- read_csv(odir("analysis_15_congruence_conflict",
                    "match_conflict_screen.csv"), show_col_types = FALSE) %>%
  transmute(Outcome = lab_out(outcome), Index = esc(label),
            Estimate = num3(estimate), SE = num3(std.error),
            `$p$` = pfmt(p.value), `$q$` = pfmt(q_fdr))
si_table(s4, "si_match_conflict.tex")

# ---- S5: eventfulness pathway (A16) -----------------------------------------
s5a <- read_csv(odir("analysis_16_eventfulness_pathway",
                     "eventfulness_a_paths.csv"), show_col_types = FALSE) %>%
  transmute(Cue = esc(label), Family = lab_var(family),
            `Slope on ISO-E` = num3(estimate), SE = num3(std.error),
            `$p$` = pfmt(p.value), `$q$` = pfmt(q_fdr))
si_table(s5a, "si_pathway_cues.tex")
boot5 <- read_csv(odir("analysis_16_eventfulness_pathway",
                       "eventfulness_indirect_bootstrap.csv"), show_col_types = FALSE)
s5b <- read_csv(odir("analysis_16_eventfulness_pathway", "eventfulness_pathways.csv"),
                show_col_types = FALSE) %>%
  filter(outcome == "overall_soundscape_z") %>%
  left_join(select(boot5, cue, indirect_ci_lo, indirect_ci_hi), by = "cue") %>%
  transmute(Cue = esc(label), `Indirect effect` = num3(indirect),
            `Bootstrap 95\\% CI` = ifelse(is.na(indirect_ci_lo), "--",
                                          sprintf("[%s, %s]", num3(indirect_ci_lo),
                                                  num3(indirect_ci_hi))))
stopifnot(nrow(s5b) == 8, sum(s5b$`Bootstrap 95\\% CI` != "--") == nrow(boot5))
si_table(s5b, "si_pathway_indirect.tex")

# ---- S6: greenery internal structure (A18) ----------------------------------
s6 <- read_csv(odir("analysis_18_greenery_composition", "quality_screen.csv"),
               show_col_types = FALSE) %>%
  transmute(Outcome = lab_out(outcome), Predictor = lab_var(quality),
            Estimate = num3(estimate), SE = num3(std.error),
            `$p$` = pfmt(p.value), `$q$` = pfmt(q.value))
si_table(s6, "si_greenery_structure.tex")

# ---- S7: constructed-context battery (A13) ----------------------------------
s7 <- read_csv(odir("analysis_13_context_composites", "context_screen.csv"),
               show_col_types = FALSE) %>%
  transmute(Outcome = lab_out(outcome), Parameter = lab_var(parameter),
            Estimate = num3(estimate), SE = num3(std.error),
            `$p$` = pfmt(p.value), `$q$` = pfmt(q_fdr))
si_table(s7, "si_context_battery.tex")

# ---- S8: photo-metric validation and chromatic phenology (A24/A25) ----------
s8a <- read_csv(odir("analysis_24_metric_validation",
                     "validation_correlations.csv"), show_col_types = FALSE) %>%
  transmute(Pair = esc(label), `$n$` = nfmt(n), `$r$` = num3(r_raw),
            `$r$ within sites` = num3(r_within))
si_table(s8a, "si_photometric_validation.tex")
s8b <- read_csv(odir("analysis_25_chromatic_phenology",
                     "season_effects_on_metrics.csv"), show_col_types = FALSE) %>%
  transmute(Metric = recode(esc(label), !!!METRIC_LAB), Type = str_to_sentence(type),
            `Spring $-$ autumn` = num3(est_spring_minus_autumn), SE = num3(se),
            `$p$` = num3(p), `Metric SD` = num3(sd_var),
            `In SD units` = num3(std_est), `$q$` = num3(q))
si_table(s8b, "si_chromatic_seasons.tex")

# ---- S9: seasonal interaction families (A17/A20/A22) ------------------------
s9a <- read_csv(odir("analysis_17_season_phenology",
                     "greenery_season_interactions.csv"),
                show_col_types = FALSE) %>%
  transmute(Outcome = lab_out(outcome), `Autumn slope` = num3(b_autumn),
            `$\\times$ spring` = num3(b_x_spring), `$p$` = pfmt(p_x_spring),
            `$\\times$ summer` = num3(b_x_summer),
            `Joint $p$` = pfmt(p_joint), `Joint $q$` = pfmt(q_joint))
si_table(s9a, "si_greenery_interactions.tex")
s9b <- read_csv(odir("analysis_20_auditory_phenology",
                     "source_season_interactions.csv"),
                show_col_types = FALSE) %>%
  transmute(Outcome = lab_out(outcome), Predictor = lab_comp(sub("_z$", "", predictor)),
            `$\\times$ spring` = num3(b_x_spring),
            `$\\times$ summer` = num3(b_x_summer),
            `Joint $p$` = pfmt(p_joint), `Joint $q$` = qfmt(q_joint))
si_table(s9b, "si_source_interactions.tex")
if (file.exists(odir("analysis_22_weather_aqi", "season_weather_adjusted.csv"))) {
  s9c <- read_csv(odir("analysis_22_weather_aqi", "season_weather_adjusted.csv"),
                  show_col_types = FALSE) %>%
    transmute(Interaction = lab_effect(label), `$\\times$ spring` = num3(b_x_spring),
              `$p$` = num3(p_x_spring), `Joint $p$` = num3(p_joint))
  si_table(s9c, "si_weather_adjusted.tex")
} else message("  analysis 22 outputs absent: si_weather_adjusted.tex not written")

# ---- S10: season-stratified slopes (A17/A20) --------------------------------
s10a <- read_csv(odir("analysis_17_season_phenology", "stratified_slopes.csv"),
                 show_col_types = FALSE) %>%
  transmute(Season = str_to_sentence(season), `$n$` = n, Effect = lab_effect(effect),
            Estimate = num3(estimate), SE = num3(std.error), `$p$` = pfmt(p.value))
si_table(s10a, "si_stratified_visual.tex", align = "lrlrrr")
s10b <- read_csv(odir("analysis_20_auditory_phenology",
                      "stratified_source_slopes.csv"), show_col_types = FALSE) %>%
  transmute(Season = str_to_sentence(season), `$n$` = n,
            Outcome = lab_out(outcome), Predictor = lab_comp(sub("_z$", "", predictor)),
            Estimate = num3(estimate), SE = num3(std.error), `$p$` = pfmt(p.value))
si_table(s10b, "si_stratified_sources.tex")

# ---- S11: condition-metric scan (A26) ---------------------------------------
s11 <- read_csv(odir("analysis_26_quality_outcome_scan",
                     "quality_outcome_scan.csv"), show_col_types = FALSE) %>%
  transmute(Metric = recode(esc(label), !!!METRIC_LAB), Block = esc(block),
            Outcome = lab_out(outcome),
            Estimate = num3(est), SE = num3(se),
            `$p$` = pfmt(p), `$q$` = pfmt(q))
si_table(s11, "si_condition_scan.tex")

# ---- S12: blueness integration (A29) ----------------------------------------
s12 <- bind_rows(
  read_csv(odir("analysis_29_blueness_integration", "blueness_sound_adjusted.csv"),
           show_col_types = FALSE) %>%
    transmute(Specification = paste0("Sound-adjusted, ", outcome_label),
              `$n$` = nfmt(n), Estimate = num3(est), SE = num3(se), `$p$` = pfmt(p)),
  read_csv(odir("analysis_29_blueness_integration",
                "blueness_season_stratified_oa.csv"), show_col_types = FALSE) %>%
    transmute(Specification = paste0(str_to_sentence(season), " subset, OA"),
              `$n$` = nfmt(n), Estimate = num3(est), SE = num3(se), `$p$` = pfmt(p)),
  read_csv(odir("analysis_29_blueness_integration",
                "blueness_loudness_compensation.csv"), show_col_types = FALSE) %>%
    transmute(Specification = paste0("$\\times$ perceived loudness, ", outcome_label),
              `$n$` = nfmt(n), Estimate = num3(b_interaction),
              SE = num3(se_interaction), `$p$` = pfmt(p_interaction)))
si_table(s12, "si_blueness_integration.tex")

# ---- S13: loudness moderation (A19) -----------------------------------------
s13a <- read_csv(odir("analysis_19_loudness_conditioned_vision",
                      "perceived_loudness_moderation.csv"),
                 show_col_types = FALSE) %>%
  transmute(Effect = lab_effect(effect), Interaction = num3(b_interaction),
            SE = num3(se), `$p$` = pfmt(p.value), `$q$` = pfmt(q.value))
si_table(s13a, "si_loudness_moderation.tex")
s13b <- read_csv(odir("analysis_19_loudness_conditioned_vision",
                      "laeq_moderation.csv"), show_col_types = FALSE) %>%
  transmute(Effect = lab_effect(effect), Interaction = num3(b_interaction),
            SE = num3(se), `$p$` = num3(p.value), `$q$` = num3(q.value))
si_table(s13b, "si_laeq_moderation.tex")
s13c <- read_csv(odir("analysis_19_loudness_conditioned_vision",
                      "tertile_stratified_slopes.csv"), show_col_types = FALSE) %>%
  transmute(Tertile = str_to_sentence(laeq_tertile), `$n$` = n,
            Effect = lab_effect(effect), Estimate = num3(estimate),
            SE = num3(std.error), `$p$` = pfmt(p.value))
si_table(s13c, "si_laeq_tertile_slopes.tex", align = "lrlrrr")

# ---- S14: elements and opportunity conditioning (A11/A14) -------------------
s14a <- read_csv(odir("analysis_11_element_candidates", "element_scan.csv"),
                 show_col_types = FALSE) %>%
  transmute(Element = lab_var(element), Outcome = lab_out(outcome),
            Estimate = num3(estimate), SE = num3(std.error),
            `$p$` = pfmt(p.value), `$q$` = pfmt(q_fdr))
si_table(s14a, "si_element_scan.tex")
s14b <- read_csv(odir("analysis_14_opportunity_conditioned",
                      "opportunity_cluster_bootstrap.csv"),
                 show_col_types = FALSE) %>%
  transmute(Effect = lab_effect(effect),
            `95\\% CI` = sprintf("[%s, %s]", num3(ci_lo), num3(ci_hi)),
            `Excludes zero` = ifelse(excludes_zero, "yes", "no"))
si_table(s14b, "si_opportunity_conditioning.tex")

# ---- S15: audience generality and date covariates (A21/A22) -----------------
moderation_table <- function(analysis, file) {
  read_csv(odir(analysis, file), show_col_types = FALSE) %>%
    transmute(Effect = lab_effect(effect), Interaction = num3(b_interaction),
              SE = num3(se), `$p$` = num3(p.value), `$q$` = num3(q.value))
}
si_table(moderation_table("analysis_21_familiarity_subgroups", "non_local_moderation.csv"),
         "si_familiarity.tex")
si_table(moderation_table("analysis_21_familiarity_subgroups", "activity_state_moderation.csv"),
         "si_activity.tex")
if (file.exists(odir("analysis_22_weather_aqi", "date_covariate_mains.csv"))) {
  s15c <- read_csv(odir("analysis_22_weather_aqi", "date_covariate_mains.csv"),
                   show_col_types = FALSE) %>%
    transmute(Outcome = lab_out(outcome), Covariate = unname(WX_LAB[wxvar]),
              Estimate = num3(estimate), SE = num3(std.error),
              `$p$` = num3(p.value), `$q$` = num3(q.value))
  si_table(s15c, "si_weather_mains.tex")
} else message("  analysis 22 outputs absent: si_weather_mains.tex not written")

# ---- S27: refits of the reported effects (A31, A26, A32, A33) ----------------
EFFECT_ORDER <- c("vehicles_isoe", "blueness_oa", "water_oa", "water_isop",
                  "flowers_oa", "water_x_pl_oa", "traffic_congruence_oa",
                  "balance_match_isoe")
SPEC_ORDER <- c("as reported", "plus panorama-group random effect",
                "within-site component", "between-site component",
                "natural-sound salience added", "all auditory measures added",
                "source salience ratings added")
refit_rows <- bind_rows(
  read_csv(odir("analysis_31_panorama_clustering", "panorama_sensitivity.csv"),
           show_col_types = FALSE) %>%
    transmute(effect = sub(" \\(A[0-9]+\\)$", "", effect),
              spec = recode(spec, `plus panorama RE` = "plus panorama-group random effect"),
              estimate, se, p.value, n),
  read_csv(odir("analysis_26_quality_outcome_scan", "quality_outcome_scan.csv"),
           show_col_types = FALSE) %>%
    filter(var == "sky_blue_ratio", outcome == "overall_soundscape_z") %>%
    transmute(effect = "blueness_oa", spec = "as reported", estimate = est, se, p.value = p, n),
  read_csv(odir("analysis_32_within_between", "within_between.csv"),
           show_col_types = FALSE) %>%
    transmute(effect, spec = paste(component, "component"), estimate, se, p.value, n),
  read_csv(odir("analysis_33_water_sound_adjusted", "water_sound_adjusted.csv"),
           show_col_types = FALSE) %>%
    filter(specification != "as reported") %>%
    transmute(effect = recode(outcome, overall_soundscape_z = "water_oa", isop_z = "water_isop"),
              spec = recode(specification,
                            `natural-sound saliency adjusted` = "natural-sound salience added",
                            `all auditory measures adjusted` = "all auditory measures added"),
              estimate, se, p.value, n),
  read_csv(odir("analysis_33_water_sound_adjusted", "water_loudness_sound_adjusted.csv"),
           show_col_types = FALSE) %>%
    filter(specification != "as reported") %>%
    transmute(effect = "water_x_pl_oa", spec = "source salience ratings added",
              estimate, se, p.value, n))
stopifnot(all(refit_rows$effect %in% EFFECT_ORDER), all(refit_rows$spec %in% SPEC_ORDER))
s16 <- refit_rows %>%
  arrange(match(effect, EFFECT_ORDER), match(spec, SPEC_ORDER)) %>%
  transmute(Effect = ifelse(duplicated(effect), "", lab_effect(effect)),
            Specification = str_to_sentence(spec),
            Estimate = num3(estimate), SE = num3(se),
            `$p$` = pfmt(p.value), `$n$` = nfmt(n))
si_table(s16, "si_effect_refits.tex", align = "llrrrr")

message("SI tables written to ", SI_DIR)
