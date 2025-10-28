# --- Libraries ---------------------------------------------------------------
library(mice)
library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(ggplot2)

# ==== Lean EPDS_9 imputation setup (keep IDs, exclude from imputation) =========


# --- Inputs (as you had) -------------------------------------------------------
trt_col   <- "medgroup"
ref_label <- "Sulfadoxine-pryimethamine dose"
id_vars   <- c("midm","mide")  # keep but exclude from imputation

# Vars used to impute EPDS_9 under MAR
predictors_for_EPDS9 <- c(
  trt_col, "v01bas_age_2", "v01bas_canread_4", "v01bas_canwrite_5",
  "v01bas_livingchd_cat", "v01bas_weight_14", "v01bas_gestage_17"
)

# --- 1) Trim dataset to only what you need ------------------------------------
dat <- tipping_maternal_db %>%
  select(EPDS_9, all_of(trt_col), all_of(predictors_for_EPDS9), any_of(id_vars))

# Types + reference
dat[[trt_col]]           <- relevel(factor(dat[[trt_col]]), ref = ref_label)
dat$v01bas_canread_4     <- factor(dat$v01bas_canread_4)
dat$v01bas_canwrite_5    <- factor(dat$v01bas_canwrite_5)
dat$v01bas_livingchd_cat <- ordered(dat$v01bas_livingchd_cat)   # ordered factor
dat$EPDS_9               <- as.factor(dat$EPDS_9)
dat$was_missing_EPDS9    <- is.na(dat$EPDS_9)

# --- 2) mice setup: impute EPDS_9 + the one missing predictor ------------------
tmp  <- mice(dat, m = 1, maxit = 0)
meth <- tmp$method
pred <- tmp$predictorMatrix

# start clean
meth[]  <- ""
pred[,] <- 0

# Outcome: EPDS_9 via logistic regression, predicted by your chosen covariates
meth["EPDS_9"] <- "logreg"
pred["EPDS_9", intersect(predictors_for_EPDS9, colnames(dat))] <- 1

# Also impute the ONLY predictor with missingness: v01bas_livingchd_cat (ordered)
meth["v01bas_livingchd_cat"] <- "polr"
# sensible predictors for it (use the others, which are complete)
pred["v01bas_livingchd_cat", c(trt_col, "v01bas_age_2", "v01bas_canread_4",
                               "v01bas_canwrite_5", "v01bas_weight_14",
                               "v01bas_gestage_17")] <- 1

# IDs: keep but exclude from imputation/prediction
for (v in intersect(id_vars, colnames(pred))) {
  meth[v]   <- ""
  pred[v, ] <- 0
  pred[, v] <- 0
}

# Never use the missingness flag inside mice
if ("was_missing_EPDS9" %in% colnames(pred)) {
  meth["was_missing_EPDS9"]   <- ""
  pred["was_missing_EPDS9", ] <- 0
  pred[, "was_missing_EPDS9"] <- 0
}

# (Optional) limit imputation to rows where a var is missing
where <- matrix(FALSE, nrow(dat), ncol(dat), dimnames = list(NULL, names(dat)))
where[, "EPDS_9"]               <- is.na(dat$EPDS_9)
where[, "v01bas_livingchd_cat"] <- is.na(dat$v01bas_livingchd_cat)

# Run mice
m <- 20; set.seed(817)
imp <- mice(dat, m = m, method = meth, predictorMatrix = pred, where = where, maxit = 20, seed = 817)

# Inspect
print(imp$loggedEvents)

# (Optional) sanity: quick NA counts
colSums(is.na(dat[, c("EPDS_9","v01bas_livingchd_cat","v01bas_age_2","v01bas_canread_4",
                      "v01bas_canwrite_5","v01bas_weight_14","v01bas_gestage_17")]))


# --- 2) δ grid ---------------------------------------------------------------
delta_grid <- seq(-2, 2, by = 0.25)

# --- 3) Models & helpers -----------------------------------------------------
# Ensure reference level in any data.frame before fitting
set_ref <- function(d) {
  d[[trt_col]] <- stats::relevel(factor(d[[trt_col]]), ref = ref_label)
  d
}
form_pred <- reformulate(
  c(trt_col, "v01bas_age_2","v01bas_canread_4","v01bas_canwrite_5",
    "v01bas_livingchd_cat","v01bas_weight_14","v01bas_gestage_17"),
  response = "EPDS_9"
)
form_bin_unadj <- reformulate(trt_col, response = "EPDS_9")

# Apply a common δ shift to imputed outcomes in NON-reference arms
adjust_imputed_binary <- function(d, delta) {
  d[[trt_col]] <- stats::relevel(factor(d[[trt_col]]), ref = ref_label)
  d_obs <- d[!d$was_missing_EPDS9, , drop = FALSE]
  fit_pred <- suppressWarnings(glm(form_pred, data = d_obs, family = binomial()))
  eta_hat  <- suppressWarnings(predict(fit_pred, newdata = d, type = "link"))
  p_del    <- plogis(eta_hat + delta)
  p_del    <- pmin(pmax(p_del, 1e-8), 1 - 1e-8)
  idx <- which(d$was_missing_EPDS9 & d[[trt_col]] != ref_label & is.finite(p_del) & !is.na(p_del))
  if (length(idx)) d$EPDS_9[idx] <- rbinom(length(idx), 1, p_del[idx])
  d
}

extract_contrasts <- function(fit) {
  broom::tidy(fit) %>%
    filter(grepl(paste0("^", trt_col), term)) %>%
    transmute(
      contrast = sub(paste0("^", trt_col), "", term),
      est_vec  = estimate,
      se_vec   = std.error
    )
}

rubin_pool <- function(ests) {
  ests %>%
    summarise(
      m_eff = n(),
      Qbar  = mean(est_vec),
      W     = mean(se_vec^2),
      B     = var(est_vec),
      Tvar  = W + (1 + 1/m_eff) * B,
      se    = sqrt(pmax(Tvar, 0)),
      lwr   = Qbar - qnorm(0.975) * se,
      upr   = Qbar + qnorm(0.975) * se,
      pval  = 2 * pnorm(-abs(Qbar / se)),
      .groups = "drop"
    )
}

# Run δ for SP-referenced contrasts (Full vs SP; Split vs SP)
run_one_delta <- function(delta) {
  ests <- map_dfr(seq_len(m), function(i) {
    d_i <- complete(imp, i); d_i$was_missing_EPDS9 <- dat$was_missing_EPDS9
    d_i <- adjust_imputed_binary(d_i, delta)
    fit <- suppressWarnings(glm(form_bin_unadj, data = d_i, family = binomial()))
    out <- extract_contrasts(fit); out$.imp <- i; out
  }) %>% filter(is.finite(est_vec), is.finite(se_vec))
  if (!nrow(ests)) return(tibble())
  ests %>% group_by(contrast) %>%
    rubin_pool() %>%
    transmute(contrast, m_eff, estimate = Qbar, W, B, Tvar, se, lwr, upr, pval, delta = delta)
}

tp_bin_results <- map_dfr(delta_grid, run_one_delta) %>% arrange(contrast, delta)

# --- 4) Pairwise Full vs Split (shift one arm at a time) ---------------------
pair_levels <- c("Mefloquine full dose", "Mefloquine split dose")

adjust_imputed_binary_perarm <- function(d, delta_map) {
  d[[trt_col]] <- stats::relevel(factor(d[[trt_col]]), ref = ref_label)
  d_obs <- d[!d$was_missing_EPDS9, , drop = FALSE]
  fit   <- suppressWarnings(glm(form_pred, data = d_obs, family = binomial()))
  eta0  <- suppressWarnings(predict(fit, newdata = d, type = "link"))
  p0    <- plogis(eta0); p_del <- p0
  for (arm in names(delta_map)) {
    idx <- d$was_missing_EPDS9 & d[[trt_col]] == arm
    if (any(idx)) p_del[idx] <- plogis(qlogis(p0[idx]) + delta_map[[arm]])
  }
  p_del <- pmin(pmax(p_del, 1e-8), 1 - 1e-8)
  idx_all <- which(d$was_missing_EPDS9 & d[[trt_col]] %in% names(delta_map) & !is.na(p_del))
  if (length(idx_all)) d$EPDS_9[idx_all] <- rbinom(length(idx_all), 1, p_del[idx_all])
  d
}

extract_pair_contrast <- function(fit, A = pair_levels[1], B = pair_levels[2]) {
  b <- coef(fit); V <- vcov(fit)
  nA <- paste0(trt_col, A); nB <- paste0(trt_col, B)
  if (!(nA %in% names(b)) || !(nB %in% names(b)))
    return(tibble(contrast = paste(A, "vs", B), est_vec = NA_real_, se_vec = NA_real_))
  est <- b[nA] - b[nB]
  var <- V[nA, nA] + V[nB, nB] - 2 * V[nA, nB]
  tibble(contrast = paste(A, "vs", B), est_vec = est, se_vec = sqrt(max(var, 0)))
}

run_one_delta_pair <- function(delta, shift_arm = pair_levels[1]) {
  ests <- map_dfr(seq_len(m), function(i) {
    d_i <- complete(imp, i); d_i$was_missing_EPDS9 <- dat$was_missing_EPDS9
    d_i <- adjust_imputed_binary_perarm(d_i, delta_map = setNames(delta, shift_arm))
    fit <- suppressWarnings(glm(form_bin_unadj, data = d_i, family = binomial()))
    out <- extract_pair_contrast(fit); out$.imp <- i; out
  }) %>% filter(is.finite(est_vec), is.finite(se_vec))
  if (!nrow(ests)) return(tibble())
  ests %>% rubin_pool() %>%
    transmute(contrast = paste(pair_levels, collapse = " vs "),
              shifted_arm = shift_arm,
              m_eff, estimate = Qbar, W, B, Tvar, se, lwr, upr, pval, delta)
}

delta_grid_pair <- seq(-2, 2, by = 0.25)
tp_pair_results_full  <- map_dfr(delta_grid_pair, ~ run_one_delta_pair(.x, shift_arm = "Mefloquine full dose"))
tp_pair_results_split <- map_dfr(delta_grid_pair, ~ run_one_delta_pair(.x, shift_arm = "Mefloquine split dose"))

# --- 5) Build one long DF with 4 panels --------------------------------------
df_full_vs_sp  <- tp_bin_results %>%
  filter(contrast == "Mefloquine full dose") %>%
  mutate(panel = "Full vs SP (δ on MQ imputed)")

df_split_vs_sp <- tp_bin_results %>%
  filter(contrast == "Mefloquine split dose") %>%
  mutate(panel = "Split vs SP (δ on MQ imputed)")

df_pair_full_shift  <- tp_pair_results_full  %>% mutate(panel = "Full vs Split (shift Full)")
df_pair_split_shift <- tp_pair_results_split %>% mutate(panel = "Full vs Split (shift Split)")

df_all <- bind_rows(
  df_full_vs_sp  %>% select(delta, estimate, lwr, upr, pval, panel),
  df_split_vs_sp %>% select(delta, estimate, lwr, upr, pval, panel),
  df_pair_full_shift  %>% select(delta, estimate, lwr, upr, pval, panel),
  df_pair_split_shift %>% select(delta, estimate, lwr, upr, pval, panel)
) %>%
  mutate(
    sig   = pval < 0.05,
    panel = factor(panel, levels = c(
      "Full vs SP (δ on MQ imputed)",
      "Split vs SP (δ on MQ imputed)",
      "Full vs Split (shift Full)",
      "Full vs Split (shift Split)"
    ))
  )

# --- 6) Tipping points per panel (first δ where 95% CI excludes 0) -----------
# Helper: return first positive and first negative tipping deltas for one panel
# Auto-detect delta and CI columns, then find first tipping delta in +/– directions
tip_one <- function(d) {
  nm <- names(d)
  
  # detect delta column
  delta_candidates <- c("delta","delta_shift","d","delta_v","shift","logodds_shift")
  delta_col <- intersect(delta_candidates, nm)[1]
  if (is.na(delta_col)) stop("No delta column found. Tried: ", paste(delta_candidates, collapse=", "))
  
  # detect CI columns
  lower_candidates <- c("lcl","ci_low","lo","lower","ymin","lwr","low","cil","ci_l")
  upper_candidates <- c("ucl","ci_high","hi","upper","ymax","upr","high","ciu","ci_u")
  lo_col <- intersect(lower_candidates, nm)[1]
  hi_col <- intersect(upper_candidates, nm)[1]
  if (is.na(lo_col) || is.na(hi_col)) {
    stop("No CI columns found. Tried lower: ",
         paste(lower_candidates, collapse=", "),
         " | upper: ", paste(upper_candidates, collapse=", "))
  }
  
  d <- d %>%
    mutate(
      ..delta = as.numeric(.data[[delta_col]]),
      ..lo    = as.numeric(.data[[lo_col]]),
      ..hi    = as.numeric(.data[[hi_col]])
    ) %>%
    arrange(..delta)
  
  excl0 <- (d$..lo > 0) | (d$..hi < 0)
  
  # first δ ≥ 0 where CI excludes 0
  pos <- d %>%
    filter(..delta >= 0, excl0) %>%
    summarize(delta_pos = if (n() > 0) min(..delta) else NA_real_) %>%
    pull(delta_pos)
  
  # first δ ≤ 0 (closest to 0 from the right)
  neg <- d %>%
    filter(..delta <= 0, excl0) %>%
    summarize(delta_neg = if (n() > 0) max(..delta) else NA_real_) %>%
    pull(delta_neg)
  
  tibble(delta_pos = pos, delta_neg = neg)
}

# --- your pipeline (keep as-is) -----------------------------------------------
tips_4 <- df_all %>%
  group_by(panel) %>%
  group_modify(~ tip_one(.x)) %>%
  ungroup() %>%
  mutate(
    first_tip_delta = case_when(
      !is.na(delta_pos) & !is.na(delta_neg) ~ if_else(abs(delta_pos) <= abs(delta_neg), delta_pos, delta_neg),
      !is.na(delta_pos) ~ delta_pos,
      !is.na(delta_neg) ~ delta_neg,
      TRUE ~ NA_real_
    ),
    first_tip_dir = case_when(
      is.na(first_tip_delta) ~ "No tipping in range",
      first_tip_delta > 0    ~ "CI > 0",
      TRUE                   ~ "CI < 0"
    )
  )

print(tips_4 %>% mutate(across(c(delta_pos, delta_neg, first_tip_delta), ~ round(., 2))))

tp_vlines4 <- tips_4 %>%
  pivot_longer(c(delta_pos, delta_neg), names_to = "type", values_to = "delta_v") %>%
  filter(!is.na(delta_v))

# --- 6b) Display both positive and negative δ tipping points ------------------
tips_4_summary <- tips_4 %>%
  dplyr::mutate(
    delta_pos = round(delta_pos, 2),
    delta_neg = round(delta_neg, 2),
    first_tip_delta = round(first_tip_delta, 2),
    odds_multiplier_pos = exp(delta_pos),
    odds_multiplier_neg = exp(delta_neg),
    odds_multiplier_first = exp(first_tip_delta)
  ) %>%
  dplyr::select(
    panel,
    delta_neg, odds_multiplier_neg,
    delta_pos, odds_multiplier_pos,
    first_tip_delta, odds_multiplier_first, first_tip_dir
  )

print(tips_4_summary, n = Inf)


# ==== Tipping-point interpretability: % of missing depressed at the tip =======
# Needs: imp, dat, trt_col, ref_label, form_pred, tips_4


# --- 7) Single ggplot with 4 panels ------------------------------------------
# Build per-panel y-top for label placement
# Create a y_top per panel (for placing labels)

# --- Harmonize panel names between df_all and tips_4 --------------------------
panel_map <- c(
  "Full vs SP (δ on MQ imputed)"  = "Full vs SP (δ on MQ imputed: Full & Split shifted)",
  "Split vs SP (δ on MQ imputed)" = "Split vs SP (δ on MQ imputed: Full & Split shifted)",
  "Full vs Split (shift Full)"    = "Full vs Split (shift Full)",
  "Full vs Split (shift Split)"   = "Full vs Split (shift Split)"
)

# Recode tips_4 panels to match df_all's labels
tips_4 <- tips_4 %>%
  dplyr::mutate(panel = panel_map[as.character(panel)]) %>%
  dplyr::mutate(panel = factor(panel,
                               levels = c(
                                 "Full vs SP (δ on MQ imputed: Full & Split shifted)",
                                 "Split vs SP (δ on MQ imputed: Full & Split shifted)",
                                 "Full vs Split (shift Full)",
                                 "Full vs Split (shift Split)"
                               )))

# Also ensure df_all uses the same factor levels
df_all <- df_all %>%
  dplyr::mutate(panel = factor(as.character(panel), levels = levels(tips_4$panel)))


y_tops <- df_all %>%
  dplyr::group_by(panel) %>%
  dplyr::summarise(
    y_top = if (all(is.infinite(upr) | is.na(upr))) 0 else max(upr, na.rm = TRUE),
    .groups = "drop"
  )

tp_vlines4 <- tips_4 %>%
  dplyr::select(panel, delta_pos, delta_neg) %>%
  tidyr::pivot_longer(c(delta_pos, delta_neg),
                      names_to = "type", values_to = "delta_v") %>%
  dplyr::filter(!is.na(delta_v)) %>%
  dplyr::left_join(y_tops, by = "panel") %>%
  dplyr::mutate(lbl = paste0("δ=", round(delta_v, 2)))


ylim_all <- range(c(df_all$lwr, df_all$upr), na.rm = TRUE)

p4 <- ggplot(df_all, aes(x = delta, y = estimate)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15) +
  geom_line() +
  geom_point(aes(shape = sig), size = 2) +
  geom_vline(data = tp_vlines4,
             aes(xintercept = delta_v),
             linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
  geom_text(data = tp_vlines4,
            aes(x = delta_v, y = y_top, label = lbl),
            vjust = -0.35, size = 3) +
  facet_wrap(~ panel, ncol = 2) +
  coord_cartesian(ylim = ylim_all) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  theme_bw() +
  labs(
    title    = "Tipping-point sensitivity to missing outcomes — EPDS ≥ 9 (binary)",
    subtitle = "δ represents a log-odds (≈ SD-unit) shift applied to imputed outcomes; models are unadjusted.",
    x = expression(delta~"(log-odds shift ≈ SD-unit)"),
    y = "Pooled log-odds ratio",
    shape = "p < 0.05",
    caption = "Dashed vertical lines mark the first δ (positive and negative) where the 95% CI excludes 0."
  )


print(p4)

ggplot2::ggsave("tp_epds9_four_panel_faceted.png", p4, width = 11, height = 9, dpi = 300)


# look at my pooled results with no tipping point, just imputed

# Pool results across imputations (unadjusted model)
mi_unadj <- purrr::map_dfr(seq_len(m), function(i) {
  d_i <- complete(imp, i) |> set_ref()
  fit <- suppressWarnings(glm(form_bin_unadj, data = d_i, family = binomial()))
  extract_contrasts(fit) |> dplyr::mutate(.imp = i)
}) |>
  dplyr::group_by(contrast) |>
  rubin_pool() |>
  # rename Qbar -> estimate and keep the essentials
  dplyr::transmute(
    contrast,
    estimate = Qbar, se, lwr, upr, pval,
    model = "Unadjusted (δ = 0)"
  )


# View (log-odds). If you want ORs, see below.
mi_unadj |>
  dplyr::mutate(dplyr::across(c(estimate, se, lwr, upr, pval), ~ round(., 3))) |>
  print(n = Inf)

mi_unadj_or <- mi_unadj |>
  dplyr::transmute(
    contrast, 
    OR = exp(estimate),
    CI_low = exp(lwr),
    CI_high = exp(upr),
    pval
  )

mi_unadj_or |>
  dplyr::mutate(dplyr::across(c(OR, CI_low, CI_high, pval), ~ round(., 3))) |>
  print(n = Inf)


# 1) Get MAR baseline probabilities p0 for EVERYONE by pooling across completed predictor datasets
get_p0_EPDS9_from_imp <- function(imp_obj, form_pred, trt_col, ref_label) {
  m <- imp_obj$m
  n <- nrow(mice::complete(imp_obj, 1))
  pmat <- matrix(NA_real_, nrow = n, ncol = m)
  
  for (i in seq_len(m)) {
    d_i <- mice::complete(imp_obj, i) |> set_ref()
    # Fit MAR predictor on observed EPDS_9 (predictors are complete in d_i)
    fit_i <- suppressWarnings(glm(form_pred,
                                  data = d_i[!is.na(d_i$EPDS_9), , drop = FALSE],
                                  family = binomial()))
    eta_i <- suppressWarnings(predict(fit_i, newdata = d_i, type = "link"))
    pmat[, i] <- plogis(eta_i)
  }
  list(p0 = rowMeans(pmat, na.rm = TRUE))
}

# 2) Summarize, among MISSING in a given arm, what % would be depressed at the tipping delta*
summarize_missing_arm_all <- function(d, p0_list, arm, delta_star) {
  d2 <- set_ref(d)
  p0 <- p0_list$p0
  
  idx_miss <- d2$was_missing_EPDS9 & d2[[trt_col]] == arm
  M <- sum(idx_miss)
  if (M == 0) {
    return(dplyr::tibble(
      arm = arm, M_missing = 0,
      p0_mean = NA_real_, p_tip_mean = NA_real_,
      pct_needed = NA_real_, extra_cases = NA_real_
    ))
  }
  # clip to avoid qlogis(0) or qlogis(1)
  p0_clip <- pmin(pmax(p0[idx_miss], 1e-8), 1 - 1e-8)
  p_tip   <- plogis(qlogis(p0_clip) + delta_star)
  
  dplyr::tibble(
    arm         = arm,
    M_missing   = M,
    p0_mean     = mean(p0_clip),
    p_tip_mean  = mean(p_tip),
    pct_needed  = 100 * mean(p_tip),
    extra_cases = sum(p_tip) - sum(p0_clip)
  )
}

# 3) Wrapper across your four panels, applying delta* and the correct arm(s)
# --- Fix: consistent panel names in calc_tip_requirements_all() ---------------
calc_tip_requirements_all <- function(original_dat, tips_tbl) {
  d0 <- set_ref(original_dat)
  p0_list <- get_p0_EPDS9_from_imp(imp, form_pred, trt_col, ref_label)
  
  levs <- levels(d0[[trt_col]])
  nonref_arms <- setdiff(levs, ref_label)
  full_arm    <- "Mefloquine full dose"
  split_arm   <- "Mefloquine split dose"
  
  purrr::map_dfr(seq_len(nrow(tips_tbl)), function(i) {
    panel <- tips_tbl$panel[i]
    delta_star <- tips_tbl$first_tip_delta[i]
    if (is.na(delta_star)) return(NULL)
    
    # ✅ FIX: match updated labels used in df_all and tips_4
    arms_to_shift <-
      if (panel %in% c(
        "Full vs SP (δ on MQ imputed: Full & Split shifted)",
        "Split vs SP (δ on MQ imputed: Full & Split shifted)"
      )) {
        nonref_arms
      } else if (panel == "Full vs Split (shift Full)") {
        full_arm
      } else if (panel == "Full vs Split (shift Split)") {
        split_arm
      } else character(0)
    
    purrr::map_dfr(arms_to_shift, ~ summarize_missing_arm_all(d0, p0_list, .x, delta_star)) |>
      dplyr::mutate(
        panel = panel,
        delta_star = delta_star,
        odds_multiplier = exp(delta_star)
      ) |>
      dplyr::select(
        panel, delta_star, odds_multiplier, arm, M_missing,
        p0_mean, p_tip_mean, pct_needed, extra_cases
      )
  })
}

# ---- Run + pretty print -------------------------------------------------------
tip_requirements <- calc_tip_requirements_all(dat, tips_4)

tip_requirements |>
  dplyr::mutate(dplyr::across(
    c(delta_star, odds_multiplier, p0_mean, p_tip_mean, pct_needed, extra_cases),
    ~ round(., 3)
  )) |>
  dplyr::arrange(panel, arm) |>
  print(n = Inf)

# === Clean and labeled tipping-point summary table for EPDS ≥ 9 ==============

tip_requirements_epds9 <- tip_requirements %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      delta_star > 0  ~ "Positive δ (higher odds of depression)",
      delta_star < 0  ~ "Negative δ (lower odds of depression)",
      TRUE            ~ "No tipping"
    )
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(delta_star, odds_multiplier, p0_mean, p_tip_mean, pct_needed, extra_cases),
      ~ round(., 3)
    )
  ) %>%
  dplyr::arrange(panel, direction, arm) %>%
  dplyr::select(
    Panel = panel,
    Direction = direction,
    Arm = arm,
    `δ (log-odds)` = delta_star,
    `Odds ×` = odds_multiplier,
    `Mean p(MAR)` = p0_mean,
    `Mean p(δ*)` = p_tip_mean,
    `% Depressed (δ*)` = pct_needed,
    `Δ Cases (δ* - MAR)` = extra_cases,
    `N Missing` = M_missing
  )

# Print table
print(tip_requirements_epds9, n = Inf)

# Combine the two MQ arms into one row for SP-referenced panels (optional view)
combined_view <- tip_requirements_epds9 %>%
  dplyr::mutate(is_sp_panel = grepl("^Full vs SP|^Split vs SP", Panel)) %>%
  dplyr::group_by(Panel, Direction, `δ (log-odds)`, `Odds ×`, is_sp_panel) %>%
  dplyr::summarise(
    `N Missing` = sum(`N Missing`),
    # weighted means by # missing
    `Mean p(MAR)` = sum(`Mean p(MAR)` * `N Missing`) / sum(`N Missing`),
    `Mean p(δ*)`  = sum(`Mean p(δ*)`  * `N Missing`) / sum(`N Missing`),
    `% Depressed (δ*)` = 100 * `Mean p(δ*)`,
    `Δ Cases (δ* - MAR)` = sum(`Δ Cases (δ* - MAR)`),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Panel, Direction)

print(combined_view, n = Inf)










# Export to Excel
library(openxlsx)
write.xlsx(
  tip_requirements_epds9,
  "epds9_tipping_point_results_table.xlsx",
  rowNames = FALSE
)



