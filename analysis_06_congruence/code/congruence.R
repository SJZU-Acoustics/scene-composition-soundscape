# Analysis 06 — Audio-visual congruence: is a sound evaluated better when its
# visual source is in view? Natural sound x visual nature; traffic sound x
# visible vehicles/roads; human sound x visible people.
# Outcomes: ISO-P, OA, harmony. FDR family over the 9 congruence interactions.

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_06_congruence"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(comp_traffic_vis = comp_vehicles + ph_road,
         comp_traffic_vis_z = as.numeric(scale(comp_traffic_vis)))

PAIRS <- tribble(
  ~sound, ~visual,
  "natural_sounds", "comp_nature",
  "traffic_noise", "comp_traffic_vis",
  "human_sounds", "comp_people",
)
OUTS <- c("isop", "overall_soundscape", "harmony")

cong <- pmap_dfr(PAIRS, function(sound, visual) {
  map_dfr(OUTS, function(y) {
    f <- paste0(y, "_z ~ ", sound, "_z * ", visual, "_z + wave + (1 | site)")
    fit <- lmer(as.formula(f), data = d)
    broom.mixed::tidy(fit, effects = "fixed") %>%
      filter(str_detect(term, ":")) %>%
      mutate(outcome = y, sound = sound, visual = visual, .before = 1)
  })
}) %>% mutate(q_fdr = fdr(p.value), family_size = n())
write_csv(cong, file.path(OUT, "congruence_interactions.csv"))

# Full coefficients for the natural x nature models (main-effect context)
for (y in OUTS) {
  f <- paste0(y, "_z ~ natural_sounds_z * comp_nature_z + wave + (1 | site)")
  fit <- lmer(as.formula(f), data = d)
  write_csv(broom.mixed::tidy(fit, effects = "fixed"),
            file.path(OUT, paste0("natural_x_nature_", y, ".csv")))
}

print(cong %>% select(outcome, sound, visual, estimate, p.value, q_fdr), n = 9)
