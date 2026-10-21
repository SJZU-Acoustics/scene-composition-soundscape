# =============================================================================
# fig2_quality.R -- Figure 2: quality, not quantity.
# (a) Greenery area share by season (A17): the segmented amount barely moves.
# (b) Greenery -> overall evaluation slope by season (A17 stratified LMMs).
# (c) Green pixels per unit of segmented greenery by season (A24), with the
#     site-cluster bootstrap contrast (A25): the measured quality gap.
# (d) Natural-sound and traffic-noise slopes on overall evaluation by season
#     (A20 stratified LMMs): the auditory premium has the same phenology.
# =============================================================================

source("display/style.R")

season_scale <- function() list(
  scale_x_discrete(labels = SEASON_LABELS),
  theme_pub())

qty <- read_csv(file.path(OUT_DIR, "fig2_quantity.csv"), show_col_types = FALSE) %>%
  mutate(season = factor(season, levels = SEASON_LEVELS))
p_a <- ggplot(qty, aes(season, comp_greenery, fill = season)) +
  geom_col(width = 0.56, show.legend = FALSE) +
  scale_fill_manual(values = SEASON_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 0.4, 0.1)) +
  labs(x = NULL, y = "Greenery area share") +
  season_scale()

sl <- read_csv(file.path(OUT_DIR, "fig2_greenery_slopes.csv"),
               show_col_types = FALSE) %>%
  mutate(season = factor(season, levels = SEASON_LEVELS),
         sig = p.value < 0.05)
p_b <- ggplot(sl, aes(season, estimate, colour = season)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(ymin = ci_lo(estimate, std.error),
                     ymax = ci_hi(estimate, std.error)), linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 2.2, show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  scale_colour_manual(values = SEASON_COLS, guide = "none") +
  labs(x = NULL, y = "Greenery slope on OA (SD/SD)") +
  season_scale()

gap <- read_csv(file.path(OUT_DIR, "fig2_quality_gap.csv"),
                show_col_types = FALSE) %>%
  mutate(season = factor(season, levels = SEASON_LEVELS))
p_c <- ggplot(gap, aes(season, mean_ratio, fill = season)) +
  geom_col(width = 0.56, show.legend = FALSE) +
  scale_fill_manual(values = SEASON_COLS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 0.6, 0.2)) +
  labs(x = NULL, y = "Green pixels per unit greenery") +
  season_scale()

src <- read_csv(file.path(OUT_DIR, "fig2_source_slopes.csv"),
                show_col_types = FALSE) %>%
  mutate(season = factor(season, levels = SEASON_LEVELS),
         predictor = recode(predictor, natural_sounds_z = "Natural sounds",
                            traffic_noise_z = "Traffic noise"),
         sig = p.value < 0.05)
p_d <- ggplot(src, aes(season, estimate, colour = predictor,
                       group = predictor)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(ymin = ci_lo(estimate, std.error),
                     ymax = ci_hi(estimate, std.error)),
                 position = position_dodge(width = 0.45), linewidth = 0.5) +
  geom_point(aes(shape = sig),
             position = position_dodge(width = 0.45), size = 2.2) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_colour_manual(values = c(`Natural sounds` = "#009E73",
                                 `Traffic noise` = "#D55E00"), name = NULL) +
  labs(x = NULL, y = "Source slope on OA (SD/SD)") +
  season_scale() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.02, 0.98),
        legend.justification = c(0, 1))

fig2 <- add_tags((p_a | p_b) / (p_c | p_d))
save_fig(fig2, "fig2_quality.png", "double_column", aspect = 0.78)
