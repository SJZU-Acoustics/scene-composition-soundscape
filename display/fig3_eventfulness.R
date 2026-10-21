# =============================================================================
# fig3_eventfulness.R -- Figure 3: eventfulness is the audio-visual axis.
# (a) The eight scene composites on eventfulness (ISO-E) (A05 screen, wave-adjusted
#     site-RE LMMs, 40-test BH-FDR family): one survivor, visible vehicles.
# (b) The seven audio-visual match/conflict indices on ISO-E (A15,
#     constituent-controlled LMMs, declared BH-FDR family).
# Filled points survive their family's correction (q < 0.05).
# =============================================================================

source("display/style.R")

scr <- read_csv(file.path(OUT_DIR, "fig3_isoe_screen.csv"),
                show_col_types = FALSE) %>%
  mutate(label = COMP_LABELS[composite],
         sig = q_fdr < 0.05) %>%
  arrange(estimate) %>%
  mutate(label = factor(label, levels = label))

p_a <- ggplot(scr, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_lo(estimate, std.error),
                     xmax = ci_hi(estimate, std.error)),
                 colour = "#0072B2", linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 2.1, colour = "#0072B2",
             fill = "#0072B2", show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(x = expression("Slope on "*italic("ISO-E")*" (SD/SD)"), y = NULL) +
  theme_pub()

idx <- read_csv(file.path(OUT_DIR, "fig3_indices.csv"),
                show_col_types = FALSE) %>%
  mutate(label = str_to_sentence(str_remove(label, "^idx_")),
         sig = q_fdr < 0.05) %>%
  arrange(estimate) %>%
  mutate(label = factor(label, levels = label))

p_b <- ggplot(idx, aes(estimate, label)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(xmin = ci_lo(estimate, std.error),
                     xmax = ci_hi(estimate, std.error)),
                 colour = "#D55E00", linewidth = 0.5) +
  geom_point(aes(shape = sig), size = 2.1, colour = "#D55E00",
             fill = "#D55E00", show.legend = FALSE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
  labs(x = expression("Index slope on "*italic("ISO-E")*" (SD/SD)"), y = NULL) +
  theme_pub()

fig3 <- add_tags(p_a | p_b)
save_fig(fig3, "fig3_eventfulness.png", "double_column", aspect = 0.42)
