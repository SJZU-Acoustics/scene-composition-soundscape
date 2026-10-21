# Analysis 11 — Element-level candidates: water in view, sculpture, flowers,
# visual crowding (person share). Presence/absence framing for sparse elements
# (analysis 04 gate). FDR family = 8 (4 elements x ISO-P/OA), wave + site RE.
# Plus: visual crowding vs heard human sound (joint model).

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_11_element_candidates"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(
    water_any = as.integer(comp_water > 0.005),
    sculpture_any = as.integer(replace_na(ph_sculpture, 0) +
                                 replace_na(ph_fountain, 0) > 0.001),
    flower_any = as.integer(replace_na(ph_flower, 0) > 0.001),
    person_share_z = as.numeric(scale(ph_person)))

ELEMENTS <- c("water_any", "sculpture_any", "flower_any", "person_share_z")
OUTS <- c("isop", "overall_soundscape")

scan <- map_dfr(ELEMENTS, function(e) {
  map_dfr(OUTS, function(y) {
    f <- paste0(y, "_z ~ ", e, " + wave + (1 | site)")
    fit <- lmer(as.formula(f), data = d)
    tidy_term(fit, e) %>% mutate(element = e, outcome = y, .before = 1)
  })
}) %>% mutate(q_fdr = fdr(p.value), family_size = n())
write_csv(scan, file.path(OUT, "element_scan.csv"))

write_csv(d %>% summarise(across(all_of(c("water_any", "sculpture_any",
                                          "flower_any")), mean)),
          file.path(OUT, "element_prevalence.csv"))

# Visual crowding vs heard human sound
m_crowd <- lmer(isop_z ~ person_share_z + human_sounds_z + wave + (1 | site),
                data = d)
write_csv(broom.mixed::tidy(m_crowd, effects = "fixed"),
          file.path(OUT, "crowding_vs_heard_human_isop.csv"))
m_crowd2 <- lmer(harmony_z ~ person_share_z + human_sounds_z + wave + (1 | site),
                 data = d)
write_csv(broom.mixed::tidy(m_crowd2, effects = "fixed"),
          file.path(OUT, "crowding_vs_heard_human_harmony.csv"))

print(scan %>% select(element, outcome, estimate, p.value, q_fdr), n = 8)
print(broom.mixed::tidy(m_crowd, effects = "fixed") %>%
        select(term, estimate, p.value))
