# =============================================================================
# fig1_boundary.R -- Figure 1: the visual-quantity boundary.
# (a) In-sample variance explained by the auditory measures vs the eight scene
#     composites per outcome (A08, joint LMM-frame R2 partition, n = 2,262).
# (b) Out-of-sample combined-over-sound linear increment per outcome with its
#     1,000-permutation null (A28, grouped 10-fold CV x 10 repeats).
# (c) Person-level R2 at held-out sites for the auditory measures and the
#     composites (A34, leave-one-site-out, 1,000 site-cluster bootstraps).
# Three panels in one row: ticks and legend 8 pt, axis titles 9 pt.
# =============================================================================

source("display/style.R")

BLOCK_COLS <- c(`Auditory measures` = "#0072B2", `Scene composites` = "#009E73")
theme_row3 <- function() theme_pub(base_size = 8, axis_title_size = 9)

v <- read_csv(file.path(OUT_DIR, "fig1_variance.csv"), show_col_types = FALSE) %>%
  pivot_longer(c(sound_total, visual_total), names_to = "block",
               values_to = "r2") %>%
  mutate(outcome = factor(outcome, levels = OUTCOME_LEVELS),
         block = recode(block, sound_total = "Auditory measures",
                        visual_total = "Scene composites"),
         block = factor(block, levels = names(BLOCK_COLS)))

# axis labels in outcome order; ISO-P and ISO-E italic (plotmath), hyphen kept
SHORT_LABELS <- expression("OA", italic("ISO-P"), italic("ISO-E"), "Harmony")

p_a <- ggplot(v, aes(outcome, r2, fill = block)) +
  geom_col(position = position_dodge(width = 0.62), width = 0.56) +
  scale_fill_manual(values = BLOCK_COLS, name = NULL) +
  scale_x_discrete(labels = SHORT_LABELS, expand = expansion(add = 0.25)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 0.20, 0.05)) +
  labs(x = NULL, y = expression("Explained variance"~(italic(R)^2))) +
  theme_row3() +
  theme(legend.position = "inside",
        legend.position.inside = c(0.98, 0.98),
        legend.justification = c(1, 1))

cv <- read_csv(file.path(OUT_DIR, "fig1_cv.csv"), show_col_types = FALSE) %>%
  mutate(outcome = factor(outcome, levels = OUTCOME_LEVELS),
         sig = p_perm < 0.05)

p_b <- ggplot(cv, aes(outcome, inc)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(ymin = null_lo05, ymax = null_hi95),
                 colour = "#E69F00", linewidth = 2.6, alpha = 0.55) +
  geom_point(aes(shape = sig), size = 2.0, colour = "black", fill = "black") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_x_discrete(labels = SHORT_LABELS, expand = expansion(add = 0.25)) +
  scale_y_continuous(breaks = seq(-0.01, 0.01, 0.005)) +
  labs(x = NULL,
       y = expression(Delta*italic(R)[plain("CV")]^2~"(visual increment)")) +
  theme_row3()

ho <- read_csv(file.path(OUT_DIR, "fig1_heldout.csv"), show_col_types = FALSE) %>%
  mutate(outcome = factor(outcome, levels = OUTCOME_LEVELS),
         block = recode(block, auditory = "Auditory measures",
                        composition = "Scene composites"),
         block = factor(block, levels = names(BLOCK_COLS)))

dodge <- position_dodge(width = 0.5)
p_c <- ggplot(ho, aes(outcome, r2_person, colour = block)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_linerange(aes(ymin = r2_person_lo, ymax = r2_person_hi),
                 position = dodge, linewidth = 0.5) +
  geom_point(position = dodge, size = 1.9) +
  scale_colour_manual(values = BLOCK_COLS, guide = "none") +
  scale_x_discrete(labels = SHORT_LABELS, expand = expansion(add = 0.25)) +
  scale_y_continuous(breaks = seq(0, 0.25, 0.05)) +
  labs(x = NULL, y = expression("Held-out-site"~italic(R)^2)) +
  theme_row3() +
  theme(plot.margin = margin(5.5, 9, 5.5, 5.5))

fig1 <- add_tags(p_a | p_b | p_c, size = 10)
save_fig(fig1, "fig1_boundary.png", "double_column", aspect = 0.40)
