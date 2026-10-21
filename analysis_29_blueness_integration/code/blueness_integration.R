# A29 - blueness integration batch: place the A26 survivor (sky blueness -> OA)
# on the paper's two moderating axes, and close its sound-confound gap.
#   Part 1: phenology - does sky blueness itself vary by season, and is the
#           blueness benefit season-specific? (A17/A25 machinery: multi-season
#           sites, autumn/spring, site + site:wave RE, equipment fixed.)
#   Part 2: compensation - does the blueness benefit grow under perceived
#           acoustic load like water's does (A19 machinery)?
#   Part 3: sound-adjusted robustness - A26's survivor model carried no sound
#           battery; refit with the full battery to close the confound gap.
# Pre-scoping exploratory layer (raw-photo delivery 2026-08-18): nothing here
# enters the manuscript storyline before Prof Zhang's scoping decision.
suppressPackageStartupMessages({library(tidyverse); library(lmerTest)})
source("analysis_23_photo_quality_metrics/code/photo_helpers.R")

ADIR <- "analysis_29_blueness_integration"
OUT <- file.path(ADIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(20260819)

d <- load_quality_data() %>% mutate(equip = factor(equip_proxy))

# ---- Part 1: phenology -------------------------------------------------------
season_tab <- table(as.character(d$site_label), d$season)
multi_sites <- rownames(season_tab)[rowSums(season_tab > 0) >= 2]
dm <- d %>% filter(site_label %in% multi_sites, season %in% c("autumn", "spring")) %>%
  droplevels()
cat("Part 1: multi-season sites:", paste(multi_sites, collapse = ", "),
    "| n =", nrow(dm), "\n")

# 1a. does sky blueness itself have phenology?
fit_q <- lmer(sky_blue_ratio_z ~ season + equip + (1 | site) + (1 | site:wave), data = dm)
co_q <- summary(fit_q)$coefficients
quantity <- tibble(term = "spring - autumn",
                   est = co_q["seasonspring", "Estimate"],
                   se = co_q["seasonspring", "Std. Error"],
                   p = co_q["seasonspring", "Pr(>|t|)"],
                   mean_autumn = mean(dm$sky_blue_ratio[dm$season == "autumn"]),
                   mean_spring = mean(dm$sky_blue_ratio[dm$season == "spring"]))
write_csv(quantity, file.path(OUT, "blueness_by_season.csv"))
cat("Blueness by season:\n"); print(quantity)

# 1b. blueness x season on the hedonic battery (declared family: 3 joint LRTs)
PHEN_OUTS <- c(overall_soundscape_z = "OA", isop_z = "ISO-P", harmony_z = "Harmony")
phen <- map_dfr(names(PHEN_OUTS), function(oc) {
  rhs_base <- c("sky_blue_ratio_z", "season", "comp_sky_z", "equip",
                "(1 | site)", "(1 | site:wave)")
  m0 <- lmer(reformulate(rhs_base, oc), data = dm, REML = FALSE)
  m1 <- lmer(reformulate(c("sky_blue_ratio_z * season", "comp_sky_z", "equip",
                           "(1 | site)", "(1 | site:wave)"), oc), data = dm, REML = FALSE)
  co <- summary(m1)$coefficients
  tibble(outcome = oc, outcome_label = PHEN_OUTS[oc], n = nrow(dm),
         b_autumn = co["sky_blue_ratio_z", "Estimate"],
         b_x_spring = co["sky_blue_ratio_z:seasonspring", "Estimate"],
         p_x_spring = co["sky_blue_ratio_z:seasonspring", "Pr(>|t|)"],
         p_joint = anova(m0, m1)$`Pr(>Chisq)`[2])
}) %>% mutate(q_joint = p.adjust(p_joint, "BH"))
write_csv(phen, file.path(OUT, "blueness_season_interactions.csv"))
cat("\nBlueness x season (joint LRTs, BH over 3):\n"); print(phen)

# 1c. season-stratified blueness -> OA slopes (descriptive ladder)
strat <- map_dfr(c("autumn", "spring"), function(s) {
  ds <- filter(dm, season == s)
  fit <- lmer(overall_soundscape_z ~ sky_blue_ratio_z + comp_sky_z + equip +
                (1 | site) + (1 | site:wave), data = ds)
  co <- summary(fit)$coefficients
  tibble(season = s, n = nrow(ds), est = co["sky_blue_ratio_z", "Estimate"],
         se = co["sky_blue_ratio_z", "Std. Error"],
         p = co["sky_blue_ratio_z", "Pr(>|t|)"])
})
write_csv(strat, file.path(OUT, "blueness_season_stratified_oa.csv"))
cat("\nStratified blueness -> OA:\n"); print(strat)

# ---- Part 2: compensation under perceived loudness ---------------------------
da <- d %>% filter(!is.na(laeq_a))
cat("\nPart 2: quality rows with acoustics, n =", nrow(da), "\n")
COMP_OUTS <- c(overall_soundscape_z = "OA", isop_z = "ISO-P")
comp <- map_dfr(names(COMP_OUTS), function(oc) {
  fit <- lmer(reformulate(c("sky_blue_ratio_z * perceived_loudness_z", "comp_sky_z",
                            "wave", "equip", "(1 | site)", "(1 | site:wave)"), oc),
              data = da)
  co <- summary(fit)$coefficients
  tibble(outcome = oc, outcome_label = COMP_OUTS[oc], n = nrow(da),
         b_main = co["sky_blue_ratio_z", "Estimate"],
         b_interaction = co["sky_blue_ratio_z:perceived_loudness_z", "Estimate"],
         se_interaction = co["sky_blue_ratio_z:perceived_loudness_z", "Std. Error"],
         p_interaction = co["sky_blue_ratio_z:perceived_loudness_z", "Pr(>|t|)"])
}) %>% mutate(q_interaction = p.adjust(p_interaction, "BH"))
write_csv(comp, file.path(OUT, "blueness_loudness_compensation.csv"))
cat("Blueness x perceived loudness (BH over 2):\n"); print(comp)

# 2b. PL-tertile stratified blueness -> OA slopes (descriptive ladder, A19-style)
da <- da %>% mutate(pl_tertile = cut(perceived_loudness,
                        quantile(perceived_loudness, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                        labels = c("quiet", "mid", "loud"), include.lowest = TRUE))
ladder <- map_dfr(levels(da$pl_tertile), function(s) {
  ds <- filter(da, pl_tertile == s)
  fit <- tryCatch(lmer(overall_soundscape_z ~ sky_blue_ratio_z + comp_sky_z +
                         wave + equip + (1 | site) + (1 | site:wave), data = ds),
                  error = function(e) NULL)
  if (is.null(fit)) return(tibble(pl_tertile = s, n = nrow(ds), est = NA_real_))
  co <- summary(fit)$coefficients
  tibble(pl_tertile = s, n = nrow(ds), est = co["sky_blue_ratio_z", "Estimate"],
         se = co["sky_blue_ratio_z", "Std. Error"],
         p = co["sky_blue_ratio_z", "Pr(>|t|)"])
})
write_csv(ladder, file.path(OUT, "blueness_pl_tertile_oa.csv"))
cat("\nBlueness -> OA by perceived-loudness tertile:\n"); print(ladder)

# ---- Part 3: sound-adjusted robustness of the A26 survivor -------------------
SOUND <- c("traffic_noise_z", "human_sounds_z", "natural_sounds_z",
           "other_noise_z", "perceived_loudness_z")
ADJ_OUTS <- c(overall_soundscape_z = "OA", isop_z = "ISO-P", harmony_z = "Harmony")
adj <- map_dfr(names(ADJ_OUTS), function(oc) {
  lab <- ADJ_OUTS[oc]
  dd <- d %>% drop_na(all_of(c(oc, "sky_blue_ratio_z")))
  fit <- lmer(reformulate(c("sky_blue_ratio_z", "comp_sky_z", SOUND, "wave", "equip",
                            "(1 | site)", "(1 | panorama_group)"), oc), data = dd)
  co <- summary(fit)$coefficients
  tibble(outcome = oc, outcome_label = lab, n = nrow(dd),
         est = co["sky_blue_ratio_z", "Estimate"],
         se = co["sky_blue_ratio_z", "Std. Error"],
         p = co["sky_blue_ratio_z", "Pr(>|t|)"])
})
write_csv(adj, file.path(OUT, "blueness_sound_adjusted.csv"))
cat("\nSound-adjusted blueness effects (A26 survivor robustness):\n"); print(adj)
