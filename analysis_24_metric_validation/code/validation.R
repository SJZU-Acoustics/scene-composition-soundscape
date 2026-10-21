# A24 - convergent validation of the A23 colour-metric layer against the
# PSPNet area shares carried in the survey workbook.
# Questions: (a) do pixel-colour vegetation/sky metrics track the segmentation
# shares? (b) does the agreement move by season (the quality-blindness test:
# segmentation counts a bare canopy as tree; colour does not)?
# Pre-scoping exploratory layer (raw-photo delivery 2026-08-18).
suppressPackageStartupMessages({library(tidyverse); library(lmerTest)})
source("analysis_23_photo_quality_metrics/code/photo_helpers.R")

ADIR <- "analysis_24_metric_validation"
OUT <- file.path(ADIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- load_quality_data() %>% mutate(equip = factor(equip_proxy))

pairs <- tribble(
  ~metric,        ~share,          ~label,
  "veg_frac",     "comp_greenery", "Green-pixel fraction vs greenery share",
  "exg_mean",     "comp_greenery", "Excess green vs greenery share",
  "sky_frac_top", "comp_sky",      "Sky-pixel fraction (top) vs sky share",
  "lum_mean",     "comp_sky",      "Brightness vs sky share",
  "edge_density", "comp_built",    "Edge density vs built share",
  "veg_frac",     "comp_nature",   "Green-pixel fraction vs nature share",
)

res <- pmap_dfr(pairs, function(metric, share, label) {
  dd <- d %>% select(all_of(c(metric, share)), site, wave, equip) %>% drop_na()
  r_raw <- cor(dd[[metric]], dd[[share]])
  # within-site, wave+equipment-adjusted partial correlation
  rm_ <- resid(lm(reformulate(c("site", "wave", "equip"), metric), dd))
  rs_ <- resid(lm(reformulate(c("site", "wave", "equip"), share), dd))
  tibble(label = label, metric = metric, share = share, n = nrow(dd),
         r_raw = r_raw, r_within = cor(rm_, rs_))
})
write_csv(res, file.path(OUT, "validation_correlations.csv"))
cat("Convergent validation:\n"); print(res)

# --- seasonal agreement: greenery share vs green pixels by season ---
by_season <- d %>% filter(season %in% c("autumn", "spring")) %>%
  group_by(season) %>%
  group_modify(function(g, ...) {
    rm_ <- resid(lm(veg_frac ~ site + equip, g))
    rs_ <- resid(lm(comp_greenery ~ site + equip, g))
    tibble(n = nrow(g), r_within = cor(rm_, rs_),
           veg_frac_mean = mean(g$veg_frac), greenery_share_mean = mean(g$comp_greenery),
           veg_sat_mean = mean(g$veg_sat, na.rm = TRUE))
  })
write_csv(by_season, file.path(OUT, "validation_by_season.csv"))
cat("\nSeasonal agreement (autumn vs spring):\n"); print(by_season)

# --- the quality-gap variable: green pixels per unit of segmented greenery ---
# ratio-style discrepancy: high greenery share with few green pixels = leafless/
# brown canopy. Defined on rows with meaningful greenery (share > 0.05).
gap <- d %>% filter(comp_greenery > 0.05) %>%
  mutate(green_per_greenery = veg_frac / comp_greenery)
gap_season <- gap %>% filter(season %in% c("autumn", "spring")) %>%
  group_by(season) %>%
  summarise(n = n(), mean_ratio = mean(green_per_greenery),
            sd = sd(green_per_greenery), .groups = "drop")
fit_gap <- lmer(green_per_greenery ~ season + equip + (1 | site) + (1 | site:wave),
                data = gap %>% filter(season %in% c("autumn", "spring")))
co <- summary(fit_gap)$coefficients
write_csv(gap_season, file.path(OUT, "green_per_greenery_by_season.csv"))
cat("\nGreen pixels per unit greenery share (>5% greenery rows):\n"); print(gap_season)
cat("season contrast (spring - autumn), equip-adjusted, site + site:wave RE:\n")
print(co["seasonspring", , drop = FALSE])

# --- figure: veg_frac vs comp_greenery coloured by season (single column) ---
WIDTH_MM <- c(single_column = 85, double_column = 178)
mm2in <- function(mm) mm / 25.4
okabe_ito <- c(autumn = "#E69F00", spring = "#009E73", summer = "#56B4E9")
theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = "Helvetica") %+replace%
    theme(axis.line = element_line(colour = "black", linewidth = 0.5),
          axis.ticks = element_line(colour = "black", linewidth = 0.4),
          axis.ticks.length = unit(0.10, "cm"),
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = base_size, colour = "black"),
          legend.text = element_text(size = base_size), legend.title = element_blank(),
          legend.position = "inside", legend.position.inside = c(0.02, 0.98),
          legend.justification = c(0, 1), legend.background = element_blank(),
          legend.key = element_blank(), panel.grid = element_blank(),
          plot.title = element_blank(),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.background = element_rect(fill = "white", colour = NA))
}
p <- ggplot(d %>% filter(season %in% c("autumn", "spring")),
            aes(comp_greenery, veg_frac, colour = season)) +
  geom_point(size = 0.5, alpha = 0.35, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_colour_manual(values = okabe_ito) +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Segmented greenery share (proportion)",
       y = "Green-pixel fraction (proportion)") +
  theme_pub()
ggsave(file.path(OUT, "fig_veg_vs_greenery_by_season.png"), p,
       width = mm2in(85), height = mm2in(85), dpi = 600, bg = "white")
cat("\nfigure written\n")
