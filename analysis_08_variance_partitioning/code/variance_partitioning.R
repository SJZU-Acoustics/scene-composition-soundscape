# Analysis 08 — Relative contribution: auditory vs visual variance shares per
# outcome, net of site + wave (commonality decomposition on partial R2).

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_08_variance_partitioning"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data()

SOUND <- c("traffic_noise_z", "human_sounds_z", "natural_sounds_z",
           "other_noise_z", "perceived_loudness_z")
VISUAL <- paste0(COMPS8, "_z")
OUTS <- c("isop", "isoe", "overall_soundscape", "harmony", "raq_mean")

partial_r2 <- function(y, extra, dd) {
  f0 <- as.formula(paste0(y, "_z ~ site + wave"))
  f1 <- as.formula(paste0(y, "_z ~ site + wave + ",
                          paste(extra, collapse = " + ")))
  m0 <- lm(f0, data = dd); m1 <- lm(f1, data = dd)
  1 - sum(resid(m1)^2) / sum(resid(m0)^2)
}

res <- map_dfr(OUTS, function(y) {
  dd <- d %>% drop_na(all_of(c(paste0(y, "_z"), SOUND, VISUAL)))
  r_s <- partial_r2(y, SOUND, dd)
  r_v <- partial_r2(y, VISUAL, dd)
  r_sv <- partial_r2(y, c(SOUND, VISUAL), dd)
  tibble(outcome = y, n = nrow(dd),
         sound_total = r_s, visual_total = r_v, joint = r_sv,
         sound_unique = r_sv - r_v, visual_unique = r_sv - r_s,
         shared = r_s + r_v - r_sv,
         visual_to_sound_ratio = r_v / r_s)
})
write_csv(res, file.path(OUT, "variance_partitioning.csv"))
print(res, n = 6)
