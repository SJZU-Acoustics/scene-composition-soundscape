# =============================================================================
# style.R -- Shared publication-plot style for the display scripts: A4 journal,
# single/double column, Helvetica 8-10 pt at print scale, left+bottom spines
# only, outward ticks, no grid, Okabe-Ito palette, PNG @ 600 dpi, opaque white
# background. In-plot text kept minimal (stats annotations + direct data labels
# only); descriptive text lives in the captions.
# Every estimate is read from the analysis outcome CSVs -- this file never
# refits models.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(ggtext)
  library(ragg)
})

# ---- Paths -------------------------------------------------------------------
P22_ROOT  <- "."  # repository root (run from there)
EXPLORE   <- "."
FIG_DIR   <- file.path(P22_ROOT, "figures")
TAB_DIR   <- file.path(P22_ROOT, "tables")
OUT_DIR   <- file.path(P22_ROOT, "output")
CLEAN_CSV <- file.path(EXPLORE,
  "analysis_01_clean_core/outcome/clean_data/clean_core_analysis.csv")
odir <- function(analysis, f) file.path(EXPLORE, analysis, "outcome", f)
for (d in c(FIG_DIR, TAB_DIR, OUT_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Width constants (mm -> in) ----------------------------------------------
WIDTH_MM  <- c(single_column = 85, double_column = 178)
mm2in     <- function(mm) as.numeric(mm) / 25.4
fig_width <- function(mode = "single_column") mm2in(WIDTH_MM[[mode]])

BASE_FAMILY <- "Helvetica"

# ---- Palettes ------------------------------------------------------------------
okabe_ito <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9", "#CC79A7",
               "#F0E442", "#000000")

# Outcomes of this paper (evaluation layer; RAQ belongs to project 21).
OUTCOME_LEVELS <- c("overall_soundscape", "isop", "isoe", "harmony")
OUTCOME_LABELS <- c(overall_soundscape = "Overall evaluation",
                    isop = "Pleasantness",
                    isoe = "Eventfulness",
                    harmony = "Harmony")
OUTCOME_COLS   <- c(overall_soundscape = "#0072B2", isop = "#56B4E9",
                    isoe = "#D55E00", harmony = "#009E73")

# Eight scene composites (A03 reduction; raw shares, the A02 decision).
COMP_LEVELS <- c("comp_greenery", "comp_sky", "comp_water", "comp_paving",
                 "comp_built", "comp_people", "comp_vehicles", "comp_furniture")
COMP_LABELS <- c(comp_greenery = "Greenery", comp_sky = "Sky",
                 comp_water = "Water", comp_paving = "Paving",
                 comp_built = "Built surface", comp_people = "People",
                 comp_vehicles = "Vehicles", comp_furniture = "Street furniture")

SEASON_LEVELS <- c("autumn", "spring", "summer")
SEASON_LABELS <- c(autumn = "Autumn", spring = "Spring", summer = "Summer")
SEASON_COLS   <- c(autumn = "#E69F00", spring = "#009E73", summer = "#56B4E9")

# ---- Site label map (English; survey-table site keys) ---------------------------
SITE_EN <- c(
  "CN_01_中山公园"       = "Zhongshan Park",
  "CN_02_奥林匹克公园"   = "Olympic Park",
  "CN_03_中山广场"       = "Zhongshan Square",
  "CN_04_大东广场"       = "Dadong Square",
  "CN_05_劝工园"         = "Quan Park",
  "CN_06_邻里公园"       = "Linli Park",
  "CN_07_长白岛森林公园" = "Changbaidao Park",
  "CN_08_工人新村"       = "Gongren Xincun",
  "CN_09_劳模园"         = "Laomo Park",
  "CN_10_罗氏圈生态公园" = "Luoshiquan Ecological Park",
  "CN_11_方塔园"         = "Fangta Park",
  "CN_12_南湖公园"       = "Nanhu Park",
  "CN_13_北陵公园"       = "Beiling Park")

# ---- Reusable theme -------------------------------------------------------------
theme_pub <- function(base_size = 9, axis_title_size = 10) {
  theme_classic(base_size = base_size, base_family = BASE_FAMILY) %+replace%
    theme(
      axis.line         = element_line(colour = "black", linewidth = 0.4),
      axis.ticks        = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.09, "cm"),
      axis.title        = element_text(size = axis_title_size, colour = "black"),
      axis.title.x      = element_text(margin = margin(t = 2.5)),
      axis.title.y      = element_text(margin = margin(r = 2.5), angle = 90),
      axis.text         = element_text(size = base_size, colour = "black"),
      legend.text       = element_text(size = base_size),
      legend.title      = element_text(size = base_size, colour = "black"),
      legend.key        = element_blank(),
      legend.key.size   = unit(0.32, "cm"),
      legend.background = element_blank(),
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.spacing = unit(2, "pt"),
      panel.grid        = element_blank(),
      plot.title        = element_blank(),
      plot.background   = element_rect(fill = "white", colour = NA),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.tag          = element_text(size = axis_title_size + 1, face = "bold",
                                       family = BASE_FAMILY, hjust = 0),
      strip.background  = element_blank(),
      strip.text        = element_text(size = base_size, colour = "black")
    )
}

add_tags <- function(p, size = 11) {
  p + plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = size, face = "bold",
                                  family = BASE_FAMILY, hjust = 0))
}

# ---- Saver ----------------------------------------------------------------------
save_fig <- function(plot, file, width_mode = "double_column",
                     height_in = NULL, aspect = 0.62, dir = FIG_DIR) {
  w <- fig_width(width_mode)
  h <- if (is.null(height_in)) w * aspect else height_in
  path <- file.path(dir, file)
  ggsave(path, plot, width = w, height = h, units = "in", dpi = 600,
         device = ragg::agg_png, bg = "white")
  message(sprintf("  wrote %-32s (%.2f x %.2f in)", file, w, h))
  invisible(path)
}

# ---- Small helpers ----------------------------------------------------------------
fmt_p <- function(p) ifelse(p < 0.001, "italic(p) < 0.001",
                            paste0("italic(p) == ", formatC(p, digits = 3, format = "f")))
ci_lo <- function(est, se) est - 1.96 * se
ci_hi <- function(est, se) est + 1.96 * se

message("style.R loaded; FIG_DIR = ", FIG_DIR)
