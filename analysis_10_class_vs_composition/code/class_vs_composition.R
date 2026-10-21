# Analysis 10 — Space class vs visual composition: is the typology just a
# visual composite in disguise, and which representation carries the signal?

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_10_class_vs_composition"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data() %>%
  mutate(is_park = as.integer(space_class == "park"))

# --- (a) Can the class be read off the photo? -------------------------------------
m_log <- glm(is_park ~ comp_greenery_z + comp_built_z + comp_sky_z,
             data = d, family = binomial)
pred <- as.integer(predict(m_log, type = "response") > 0.5)
acc <- mean(pred == d$is_park)
auc_num <- suppressWarnings({
  r <- rank(predict(m_log, type = "response"))
  n1 <- sum(d$is_park); n0 <- sum(1 - d$is_park)
  (sum(r[d$is_park == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
})
write_csv(tibble(statistic = c("accuracy", "auc"), value = c(acc, auc_num)),
          file.path(OUT, "photo_predicts_class.csv"))
write_csv(broom::tidy(m_log), file.path(OUT, "class_logit_coefficients.csv"))

# --- (b) Photo-level: greenery slope on isop with vs without class/site -------------
specs <- list(
  greenery_only = "isop_z ~ comp_greenery_z + wave",
  plus_class = "isop_z ~ comp_greenery_z + space_class + wave",
  plus_site_re = "isop_z ~ comp_greenery_z + wave + (1 | site)"
)
sl <- imap_dfr(specs, function(f, name) {
  fit <- if (str_detect(f, "\\|")) lmer(as.formula(f), data = d)
         else lm(as.formula(f), data = d)
  est <- if (inherits(fit, "lmerMod") || inherits(fit, "lmerModLmerTest"))
    tidy_term(fit, "comp_greenery_z") else
    broom::tidy(fit) %>% filter(term == "comp_greenery_z") %>%
      mutate(df = NA) %>% select(term, estimate, std.error, statistic, df, p.value)
  est %>% mutate(model = name, .before = 1)
})
write_csv(sl, file.path(OUT, "greenery_isop_by_specification.csv"))

# --- (c) Site level (12 sites, descriptive): class vs greenery as predictors of
# site-mean outcomes -----------------------------------------------------------------
site <- d %>% group_by(site_label, space_class) %>%
  summarise(across(c(comp_greenery, comp_built, isop, overall_soundscape,
                     raq_mean, natural_share), mean), n = n(), .groups = "drop")
write_csv(site, file.path(OUT, "site_level_table.csv"))
site_cor <- tibble(
  pair = c("greenery_vs_park", "isop_vs_greenery", "isop_vs_park",
           "natural_share_vs_greenery"),
  r = c(cor(site$comp_greenery, as.integer(site$space_class == "park")),
        cor(site$isop, site$comp_greenery),
        cor(site$isop, as.integer(site$space_class == "park")),
        cor(site$natural_share, site$comp_greenery)))
write_csv(site_cor, file.path(OUT, "site_level_correlations.csv"))

cat("accuracy:", round(acc, 3), " AUC:", round(auc_num, 3), "\n")
print(sl %>% select(model, estimate, p.value))
print(site_cor)
