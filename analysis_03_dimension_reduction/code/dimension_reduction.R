# Analysis 03 — Dimension reduction of the 150 ADE20K categories.
# (a) Sparsity inventory; (b) theory-driven composites coverage; (c) data-driven
# PCA cross-check on the prevalent categories; (d) composite descriptives.

source("_shared/helpers.R")

ANALYSIS_DIR <- "analysis_03_dimension_reduction"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE)

d <- load_vis_data()
ph_cols <- names(d)[str_starts(names(d), "ph_")]

# --- (a) Sparsity inventory ---------------------------------------------------------
inv <- map_dfr(ph_cols, function(v) {
  x <- d[[v]]
  tibble(category = str_remove(v, "ph_"), mean_share = mean(x, na.rm = TRUE),
         pct_nonzero = mean(x > 0, na.rm = TRUE) * 100,
         max_share = max(x, na.rm = TRUE))
}) %>% arrange(desc(mean_share))
write_csv(inv, file.path(OUT, "category_sparsity_inventory.csv"))

summ <- tibble(
  n_categories = length(ph_cols),
  n_mean_share_ge_001 = sum(inv$mean_share >= 0.01),
  n_nonzero_ge_10pct = sum(inv$pct_nonzero >= 10),
  n_nonzero_lt_1pct = sum(inv$pct_nonzero < 1),
  share_top10 = sum(inv$mean_share[1:10]),
  share_top20 = sum(inv$mean_share[1:20]))
write_csv(summ, file.path(OUT, "sparsity_summary.csv"))

# --- (b) Theory-driven composite coverage --------------------------------------------
comp_cols <- c(COMPS8, "comp_terrain", "comp_animal", "comp_other")
cov <- d %>%
  summarise(across(all_of(comp_cols), mean)) %>%
  pivot_longer(everything(), names_to = "composite", values_to = "mean_share") %>%
  arrange(desc(mean_share)) %>%
  mutate(pct_of_total = mean_share / sum(mean_share) * 100)
write_csv(cov, file.path(OUT, "composite_coverage.csv"))

# --- (c) PCA cross-check on prevalent categories (mean share >= 0.5%) ----------------
prev <- inv %>% filter(mean_share >= 0.005) %>% pull(category)
X <- d[paste0("ph_", prev)] %>%
  mutate(across(everything(), ~ replace_na(.x, 0))) %>%
  select(where(~ sd(.x) > 0)) %>%
  mutate(across(everything(), ~ as.numeric(scale(.x))))
pca <- prcomp(X)
eig <- tibble(component = seq_along(pca$sdev), eigenvalue = pca$sdev^2,
              var_explained = pca$sdev^2 / sum(pca$sdev^2))
write_csv(eig %>% slice(1:10), file.path(OUT, "pca_eigenvalues.csv"))
loads <- as_tibble(pca$rotation[, 1:4], rownames = "category")
write_csv(loads, file.path(OUT, "pca_loadings_pc1_4.csv"))

# --- (d) Composite descriptives and intercorrelations ---------------------------------
write_csv(as_tibble(cor(d[COMPS8], use = "pairwise"), rownames = "composite"),
          file.path(OUT, "composite_correlations.csv"))
by_site <- d %>% group_by(site_label, space_class) %>%
  summarise(across(all_of(COMPS8), mean), n = n(), .groups = "drop")
write_csv(by_site, file.path(OUT, "composite_site_means.csv"))

print(summ); print(cov, n = 12)
print(inv %>% slice(1:15), n = 15)
print(eig %>% slice(1:5))
