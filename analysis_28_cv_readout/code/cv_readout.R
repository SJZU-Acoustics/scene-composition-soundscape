# Analysis 28 - Cross-validated multivariate readout (estimator-family change).
# Every positive claim so far rests on per-predictor LMM slopes (mass-univariate)
# and the 30:1 boundary claim (A08) is an in-sample partial-R2 partition. Here
# the same within-site estimand is read OUT OF SAMPLE in one pipeline:
#   - grouped 10-fold CV (folds keep shared panoramas together) x 10 repeats;
#   - outcomes and predictors residualised on site + wave INSIDE each training
#     fold, the held-out fold adjusted with the training coefficients;
#   - six models per outcome: auditory / visual / combined x linear / forest,
#     forest settings (mtry, min.node.size) chosen per outcome and block from a
#     small grid by the same grouped CV;
#   - the combined-over-auditory increment calibrated by permuting the visual
#     block across panorama groups within site-by-wave cells, every permuted
#     data set scored on the same 100 folds (1,000 permutations linear, 200
#     forest).
# Condition block (2026-09-18): the nine photograph condition metrics of the
# A26 family enter the same pipeline as a third block (linear; nuisance set
# site + wave + capture system; n = 2,261), so that composition and condition
# are compared like for like: increment over auditory, and over auditory +
# composition, each against 1,000 block permutations.

source("analysis_23_photo_quality_metrics/code/photo_helpers.R")
suppressPackageStartupMessages(library(ranger))

ANALYSIS_DIR <- "analysis_28_cv_readout"
OUT <- file.path(ANALYSIS_DIR, "outcome"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SOUND <- c("traffic_noise_z", "human_sounds_z", "natural_sounds_z",
           "other_noise_z", "perceived_loudness_z")
VIS <- paste0(COMPS8, "_z")
COND <- paste0(QUALITY_VARS$var, "_z")
OUTS <- c(overall_soundscape_z = "OA", isop_z = "ISO-P",
          isoe_z = "ISO-E", harmony_z = "Harmony")
K_FOLDS <- 10; R_REPEATS <- 10
N_PERM_LM <- 1000; N_PERM_RF <- 200; N_TREES <- 500
N_THREADS <- max(1L, parallel::detectCores() - 1L)
SEED <- 20260819

# ---- data --------------------------------------------------------------------
pg <- load_photo_metrics() %>% select(id, panorama_group)
d <- load_vis_data() %>% left_join(pg, by = "id") %>%
  mutate(panorama_group = ifelse(is.na(panorama_group), paste0("solo_", id),
                                 panorama_group)) %>%
  select(id, site, wave, panorama_group, all_of(c(names(OUTS), SOUND, VIS))) %>%
  drop_na() %>% mutate(site = droplevels(factor(site)), wave = droplevels(factor(wave)))
dq <- load_quality_data() %>%
  mutate(equip = factor(equip_proxy)) %>%
  select(id, site, wave, equip, panorama_group,
         all_of(c(names(OUTS), SOUND, VIS, COND))) %>%
  drop_na() %>% mutate(site = droplevels(factor(site)), wave = droplevels(factor(wave)))
cat("Primary frame n =", nrow(d), "| sites:", n_distinct(d$site), "| waves:",
    n_distinct(d$wave), "| panorama groups:", n_distinct(d$panorama_group), "\n")
cat("Condition frame n =", nrow(dq), "| panorama groups:",
    n_distinct(dq$panorama_group), "\n")

# ---- machinery ---------------------------------------------------------------
make_folds <- function(groups, k = K_FOLDS) {
  ug <- sample(unique(groups))
  fm <- setNames(rep(seq_len(k), length.out = length(ug)), ug)
  as.integer(fm[groups])
}

# One object per (repeat, fold): training/test rows, the QR of the training
# nuisance design and the test nuisance design.
build_folds <- function(dat, nuis) {
  X <- model.matrix(reformulate(nuis), dat)
  flatten(map(seq_len(R_REPEATS), function(r) {
    set.seed(SEED + r)
    f <- make_folds(dat$panorama_group)
    map(seq_len(K_FOLDS), function(k) {
      tr <- which(f != k); te <- which(f == k)
      list(rep = r, fold = k, tr = tr, te = te,
           qrX = qr(X[tr, , drop = FALSE]), Xte = X[te, , drop = FALSE])
    })
  }))
}

# Residualise the columns of M on the nuisance design, training rows only;
# held-out rows are adjusted with the training coefficients.
resid_fold <- function(M, fo) {
  B <- qr.coef(fo$qrX, M[fo$tr, , drop = FALSE]); B[is.na(B)] <- 0
  list(tr = qr.resid(fo$qrX, M[fo$tr, , drop = FALSE]),
       te = M[fo$te, , drop = FALSE] - fo$Xte %*% B)
}

lm_predict <- function(Ptr, Ytr, Pte) {
  B <- qr.coef(qr(cbind(1, Ptr)), Ytr); B[is.na(B)] <- 0
  cbind(1, Pte) %*% B
}

# Pooled out-of-fold R2 per repeat and outcome from a list of per-fold
# squared-error and outcome sums.
r2_by_repeat <- function(folds, Yte, pred) {
  map_dfr(seq_len(R_REPEATS), function(r) {
    idx <- which(map_int(folds, "rep") == r)
    y <- do.call(rbind, Yte[idx]); p <- do.call(rbind, pred[idx])
    tibble(repeat_id = r, outcome = colnames(y),
           r2 = 1 - colSums((y - p)^2) / colSums(sweep(y, 2, colMeans(y))^2))
  })
}

# Permute a predictor block across panorama groups within strata: every group
# receives the block profile of a donor group from the same stratum, so shared
# panoramas stay shared and the nuisance structure of the block is preserved.
make_permuter <- function(dat, block, strata_vars) {
  g <- dat %>% mutate(.row = row_number()) %>%
    group_by(panorama_group) %>% slice(1) %>% ungroup() %>%
    mutate(.stratum = interaction(!!!syms(strata_vars), drop = TRUE))
  M <- as.matrix(dat[, block]); gi <- match(dat$panorama_group, g$panorama_group)
  by_stratum <- split(seq_len(nrow(g)), g$.stratum)
  function() {
    donor <- seq_len(nrow(g))
    for (s in by_stratum) if (length(s) > 1) donor[s] <- s[sample.int(length(s))]
    M[g$.row[donor[gi]], , drop = FALSE]
  }
}

# ---- primary frame: fold objects and residualised blocks ---------------------
folds <- build_folds(d, c("site", "wave"))
Y <- as.matrix(d[, names(OUTS)]); A <- as.matrix(d[, SOUND]); V <- as.matrix(d[, VIS])
Yr <- map(folds, ~ resid_fold(Y, .x)); Ar <- map(folds, ~ resid_fold(A, .x))
Vr <- map(folds, ~ resid_fold(V, .x))
Yte <- map(Yr, "te")

# ---- linear models ------------------------------------------------------------
lm_block <- function(Pr_list) {
  pred <- map2(Pr_list, Yr, ~ lm_predict(.x$tr, .y$tr, .x$te))
  r2_by_repeat(folds, Yte, pred)
}
bind_blocks <- function(...) {
  bl <- list(...)
  map(seq_along(bl[[1]]), function(i)
    list(tr = do.call(cbind, map(bl, ~ .x[[i]]$tr)),
         te = do.call(cbind, map(bl, ~ .x[[i]]$te))))
}
cv_lm <- bind_rows(lm_block(Ar) %>% mutate(model = "lm_sound"),
                   lm_block(Vr) %>% mutate(model = "lm_visual"),
                   lm_block(bind_blocks(Ar, Vr)) %>% mutate(model = "lm_both"))

# ---- forest: settings chosen by the same grouped CV (first repeat) -----------
# Forest fits run one per core (single-threaded, explicitly seeded), so results
# do not depend on scheduling.
rf_fit_predict <- function(Ptr, ytr, Pte, mtry, mns, seed) {
  fit <- ranger(x = as.data.frame(Ptr), y = ytr, num.trees = N_TREES, mtry = mtry,
                min.node.size = mns, num.threads = 1, seed = seed)
  predict(fit, data = as.data.frame(Pte), num.threads = 1)$predictions
}
par_map <- function(x, f) parallel::mclapply(x, f, mc.cores = N_THREADS)
BLOCKS <- list(sound = Ar, visual = Vr, both = bind_blocks(Ar, Vr))
rep1 <- which(map_int(folds, "rep") == 1)
grid_for <- function(p) crossing(mtry = unique(c(floor(sqrt(p)), ceiling(p / 2), p)),
                                 mns = c(5, 20, 50, 100, 200, 400))
tuning <- imap_dfr(BLOCKS, function(Pr, bname) {
  map_dfr(names(OUTS), function(oc) {
    pmap_dfr(grid_for(ncol(Pr[[1]]$tr)), function(mtry, mns) {
      pred <- par_map(rep1, function(i)
        rf_fit_predict(Pr[[i]]$tr, Yr[[i]]$tr[, oc], Pr[[i]]$te, mtry, mns, SEED + i))
      y <- unlist(map(rep1, ~ Yte[[.x]][, oc])); p <- unlist(pred)
      tibble(block = bname, outcome = oc, mtry = mtry, min_node_size = mns,
             r2 = 1 - sum((y - p)^2) / sum((y - mean(y))^2))
    })
  })
})
settings <- tuning %>% group_by(block, outcome) %>% slice_max(r2, n = 1, with_ties = FALSE) %>%
  ungroup()
write_csv(tuning, file.path(OUT, "rf_tuning.csv"))
write_csv(settings, file.path(OUT, "rf_settings.csv"))
cat("\nForest settings chosen:\n"); print(settings, n = 20)

rf_block <- function(Pr, bname, seed_offset = 0) {
  st <- settings %>% filter(block == bname)
  st <- split(st, st$outcome)
  pred <- par_map(seq_along(folds), function(i) {
    sapply(names(OUTS), function(oc)
      rf_fit_predict(Pr[[i]]$tr, Yr[[i]]$tr[, oc], Pr[[i]]$te, st[[oc]]$mtry,
                     st[[oc]]$min_node_size, SEED + seed_offset + i))
  })
  r2_by_repeat(folds, Yte, pred)
}
cv_rf <- imap_dfr(BLOCKS, ~ rf_block(.x, .y) %>% mutate(model = paste0("rf_", .y)))

cv <- bind_rows(cv_lm, cv_rf)
MODEL_ORDER <- c("lm_sound", "lm_visual", "lm_both", "rf_sound", "rf_visual", "rf_both")
cv_sum <- cv %>% group_by(outcome, model) %>%
  summarise(r2_mean = mean(r2), r2_sd = sd(r2), .groups = "drop") %>%
  mutate(outcome_label = OUTS[outcome]) %>%
  arrange(outcome, match(model, MODEL_ORDER))
write_csv(cv, file.path(OUT, "cv_repeats.csv"))
write_csv(cv_sum, file.path(OUT, "cv_skill.csv"))
print(cv_sum, n = 30)

mean_r2 <- function(tab, m) tab %>% filter(model == m) %>% group_by(outcome) %>%
  summarise(r2 = mean(r2), .groups = "drop") %>% deframe()
obs_lm <- mean_r2(cv, "lm_both") - mean_r2(cv, "lm_sound")
obs_rf <- mean_r2(cv, "rf_both") - mean_r2(cv, "rf_sound")

# ---- permutation null: visual block, all 100 folds per permutation -----------
perm_vis <- make_permuter(d, VIS, c("site", "wave"))
set.seed(SEED)
t0 <- Sys.time()
null_lm <- map_dfr(seq_len(N_PERM_LM), function(b) {
  Vp <- perm_vis()
  Vpr <- map(folds, ~ resid_fold(Vp, .x))
  r2 <- lm_block(bind_blocks(Ar, Vpr)) %>% group_by(outcome) %>%
    summarise(r2 = mean(r2), .groups = "drop") %>% deframe()
  tibble(perm = b, outcome = names(r2), inc_lm = r2 - mean_r2(cv, "lm_sound")[names(r2)])
})
cat("linear permutations done:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
write_csv(null_lm, file.path(OUT, "permutation_null_increments.csv"))

# ---- condition block (linear), same pipeline on the condition frame ----------
folds_q <- build_folds(dq, c("site", "wave", "equip"))
Yq <- as.matrix(dq[, names(OUTS)])
res_q <- function(M) map(folds_q, ~ resid_fold(M, .x))
Yqr <- res_q(Yq); Aqr <- res_q(as.matrix(dq[, SOUND]))
Cqr <- res_q(as.matrix(dq[, VIS])); Qqr <- res_q(as.matrix(dq[, COND]))
Yqte <- map(Yqr, "te")
lm_block_q <- function(Pr_list) {
  pred <- map2(Pr_list, Yqr, ~ lm_predict(.x$tr, .y$tr, .x$te))
  r2_by_repeat(folds_q, Yqte, pred)
}
bind_q <- function(...) {
  bl <- list(...)
  map(seq_along(bl[[1]]), function(i)
    list(tr = do.call(cbind, map(bl, ~ .x[[i]]$tr)),
         te = do.call(cbind, map(bl, ~ .x[[i]]$te))))
}
cvq <- bind_rows(lm_block_q(Aqr) %>% mutate(model = "auditory"),
                 lm_block_q(bind_q(Aqr, Cqr)) %>% mutate(model = "auditory_composition"),
                 lm_block_q(bind_q(Aqr, Qqr)) %>% mutate(model = "auditory_condition"),
                 lm_block_q(bind_q(Aqr, Cqr, Qqr)) %>% mutate(model = "auditory_composition_condition"))
cvq_sum <- cvq %>% group_by(outcome, model) %>%
  summarise(r2_mean = mean(r2), r2_sd = sd(r2), .groups = "drop") %>%
  mutate(outcome_label = OUTS[outcome])
write_csv(cvq_sum, file.path(OUT, "condition_block_skill.csv"))
print(cvq_sum, n = 20)

mq <- function(m) mean_r2(cvq, m)
obs_q <- tibble(outcome = names(OUTS),
                composition_over_auditory = (mq("auditory_composition") - mq("auditory"))[names(OUTS)],
                condition_over_auditory = (mq("auditory_condition") - mq("auditory"))[names(OUTS)],
                condition_over_both = (mq("auditory_composition_condition") -
                                         mq("auditory_composition"))[names(OUTS)])
perm_c <- make_permuter(dq, VIS, c("site", "wave", "equip"))
perm_q <- make_permuter(dq, COND, c("site", "wave", "equip"))
t0 <- Sys.time()
null_q <- map_dfr(seq_len(N_PERM_LM), function(b) {
  Cp <- res_q(perm_c()); Qp <- res_q(perm_q())
  r2m <- function(Pr) lm_block_q(Pr) %>% group_by(outcome) %>%
    summarise(r2 = mean(r2), .groups = "drop") %>% deframe()
  tibble(perm = b, outcome = names(OUTS),
         composition_over_auditory = (r2m(bind_q(Aqr, Cp)) - mq("auditory"))[names(OUTS)],
         condition_over_auditory = (r2m(bind_q(Aqr, Qp)) - mq("auditory"))[names(OUTS)],
         condition_over_both = (r2m(bind_q(Aqr, Cqr, Qp)) -
                                  mq("auditory_composition"))[names(OUTS)])
})
cat("condition-block permutations done:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
INCS <- c("composition_over_auditory", "condition_over_auditory", "condition_over_both")
cond_inc <- map_dfr(INCS, function(v) {
  map_dfr(names(OUTS), function(oc) {
    nl <- null_q[[v]][null_q$outcome == oc]; ob <- obs_q[[v]][obs_q$outcome == oc]
    tibble(outcome = oc, outcome_label = OUTS[oc], increment = v, n = nrow(dq), inc = ob,
           null_mean = mean(nl), null_lo05 = quantile(nl, 0.05), null_hi95 = quantile(nl, 0.95),
           p_perm = (1 + sum(nl >= ob)) / (N_PERM_LM + 1))
  })
})
write_csv(null_q, file.path(OUT, "condition_block_null.csv"))
write_csv(cond_inc, file.path(OUT, "condition_block_increments.csv"))
print(cond_inc, n = 20, width = Inf)

# ---- forest permutation null (combined model refitted on every permuted set) --
set.seed(SEED + 1000)
CKPT <- file.path(OUT, "null_rf_checkpoint.rds")
t0 <- Sys.time()
# A checkpoint from an interrupted run is resumed; the permutation draws are
# replayed so the random stream is the same as in an uninterrupted run.
null_rf <- if (file.exists(CKPT)) readRDS(CKPT) else list()
for (b in seq_len(N_PERM_RF)) {
  Vp <- perm_vis()
  if (b <= length(null_rf)) next
  Vpr <- map(folds, ~ resid_fold(Vp, .x))
  r2 <- rf_block(bind_blocks(Ar, Vpr), "both", seed_offset = 1000 * b) %>%
    group_by(outcome) %>% summarise(r2 = mean(r2), .groups = "drop") %>% deframe()
  null_rf[[b]] <- tibble(perm = b, outcome = names(r2),
                         inc_rf = r2 - mean_r2(cv, "rf_sound")[names(r2)])
  if (b %% 10 == 0) {
    saveRDS(null_rf, CKPT)
    cat("rf perm", b, "/", N_PERM_RF, "|", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
  }
}
null_rf <- bind_rows(null_rf)

null <- full_join(null_lm, null_rf, by = c("outcome", "perm"))
perm <- null %>% group_by(outcome) %>%
  summarise(p_inc_lm = (1 + sum(inc_lm >= obs_lm[outcome[1]], na.rm = TRUE)) / (N_PERM_LM + 1),
            p_inc_rf = (1 + sum(inc_rf >= obs_rf[outcome[1]], na.rm = TRUE)) / (N_PERM_RF + 1),
            null_mean_lm = mean(inc_lm, na.rm = TRUE), null_sd_lm = sd(inc_lm, na.rm = TRUE),
            null_mean_rf = mean(inc_rf, na.rm = TRUE), null_sd_rf = sd(inc_rf, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(inc_lm = obs_lm[outcome], inc_rf = obs_rf[outcome], outcome_label = OUTS[outcome])
write_csv(null, file.path(OUT, "permutation_null_increments.csv"))
write_csv(perm, file.path(OUT, "permutation_increments.csv"))
file.remove(CKPT)
print(perm, width = Inf)
