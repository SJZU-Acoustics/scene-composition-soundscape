# =============================================================================
# fig4_condition.R -- Figure 4: condition effects under acoustic load.
# (a) Water-in-view slope on overall evaluation by LAeq tertile (A19,
#     available-channel acoustic sample): compensation under load.
# (b) The sky-blueness effect on overall evaluation across specifications
#     (A26 survivor; A29 integration): robust, spring-tilted, load-independent.
# (c) Heard-and-seen traffic on overall evaluation, pooled and by season
#     (A06 congruence survivor; A17 stratified): the targeted penalty.
# Filled points: p < 0.05 (a, c strata) or the declared-family q < 0.05
# (pooled/overall rows).
# =============================================================================

source("display/style.R")

wt <- read_csv(file.path(OUT_DIR, "fig4_water_tertiles.csv"),
               show_col_types = FALSE) %>%
  mutate(laeq_tertile = factor(laeq_tertile, levels = c("quiet", "mid", "loud"),
                               labels = c("Quiet", "Mid", "Loud")),
         sig = p.value < 0.05)

p_a <- ggplot(wt, aes(laeq_tertile, estimate)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(ymin = ci_lo(estimate, std.error),
                     ymax = ci_hi(estimate, std.error)),
                 colour = "#0072B2", linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 2.2, colour = "#0072B2",
             show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(x = expression(italic(L)[plain("Aeq")]~"tertile"),
       y = "Water-in-view effect on OA (SD)") +
  theme_pub(base_size = 8, axis_title_size = 9)

bl <- read_csv(file.path(OUT_DIR, "fig4_blueness.csv"),
               show_col_types = FALSE) %>%
  mutate(label = factor(label, levels = rev(label)),
         sig = p.value < 0.05)

p_b <- ggplot(bl, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_lo(estimate, std.error),
                     xmax = ci_hi(estimate, std.error)),
                 colour = "#009E73", linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 2.1, colour = "#009E73",
             show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(x = "Sky-blueness effect on OA (SD/SD)", y = NULL) +
  theme_pub(base_size = 8, axis_title_size = 9)

cg <- read_csv(file.path(OUT_DIR, "fig4_congruence.csv"),
               show_col_types = FALSE) %>%
  mutate(stratum = factor(stratum, levels = c("pooled", "autumn", "spring"),
                          labels = c("Pooled", "Autumn", "Spring")),
         sig = p.value < 0.05)

p_c <- ggplot(cg, aes(stratum, estimate)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(ymin = ci_lo(estimate, std.error),
                     ymax = ci_hi(estimate, std.error)),
                 colour = "#D55E00", linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 2.2, colour = "#D55E00",
             show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(x = NULL, y = "Heard-and-seen traffic on OA (SD)") +
  theme_pub(base_size = 8, axis_title_size = 9)

fig4 <- add_tags(p_a | p_b | p_c)
save_fig(fig4, "fig4_condition.png", "double_column", aspect = 0.40)
