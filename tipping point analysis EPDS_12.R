# --- Libraries ---------------------------------------------------------------
library(mice)
library(dplyr)
library(purrr)
library(broom)
library(tidyr)
library(ggplot2)

# --- 0) Prep -----------------------------------------------------------------
dat <- tipping_maternal_db
trt_col   <- "medgroup"
ref_label <- "Sulfadoxine-pryimethamine dose"

# Keep ONLY variables used in the MI and analysis
dat <- dat %>%
  dplyr::select(
    EPDS_12, !!trt_col,
    v01bas_age_2, v01bas_canread_4, v01bas_canwrite_5,
    v01bas_livingchd_cat, v01bas_weight_14, v01bas_gestage_17
  )

# Types
dat[[trt_col]]              <- as.factor(dat[[trt_col]])
dat$v01bas_canread_4        <- as.factor(dat$v01bas_canread_4)
dat$v01bas_canwrite_5       <- as.factor(dat$v01bas_canwrite_5)
dat$v01bas_livingchd_cat    <- ordered(dat$v01bas_livingchd_cat)
dat$EPDS_12                 <- as.factor(dat$EPDS_12)

# Flag who was missing EPDS_12 pre-imputation (used for tipping shifts later)
dat$was_missing_epds12 <- is.na(dat$EPDS_12)

# --- 1) Multiple imputation (MAR for EPDS_12 + impute missing predictor) -----
m <- 20; set.seed(817)

# Method & predictor matrix scaffolding
tmp  <- mice(dat, m = 1, maxit = 0)
meth <- tmp$method
pred <- tmp$predictorMatrix

# Start with 'no imputation' then turn on for variables that have missingness & are needed
meth[] <- ""

# Binary outcome by logistic regression
meth["EPDS_12"] <- "logreg"

# Ordered categorical predictor has 15 missings → impute with proportional odds model
meth["v01bas_livingchd_cat"] <- "polr"

# Predictors used for EPDS_12 imputation (MAR model; unadjusted analysis comes later)
pred["EPDS_12", ] <- 0
pred["EPDS_12", c(trt_col, "v01bas_age_2","v01bas_canread_4","v01bas_canwrite_5",
                  "v01bas_livingchd_cat","v01bas_weight_14","v01bas_gestage_17")] <- 1

# Predictors used for v01bas_livingchd_cat imputation (rich but reasonable)
pred["v01bas_livingchd_cat", ] <- 0
pred["v01bas_livingchd_cat", c(trt_col, "v01bas_age_2","v01bas_canread_4","v01bas_canwrite_5",
                               "v01bas_weight_14","v01bas_gestage_17","EPDS_12")] <- 1

# Restrict imputation to variables that actually have NA
where <- matrix(FALSE, nrow(dat), ncol(dat), dimnames = list(NULL, names(dat)))
where[, "EPDS_12"]              <- is.na(dat$EPDS_12)
where[, "v01bas_livingchd_cat"] <- is.na(dat$v01bas_livingchd_cat)

# Run mice
imp <- mice(dat, m = m, method = meth, predictorMatrix = pred, where = where,
            maxit = 20, seed = 817)

# Quick checks
plot(imp)
densityplot(imp, ~ EPDS_12)

# Helper to enforce SP as reference later
set_ref <- function(d) {
  d[[trt_col]] <- stats::relevel(factor(d[[trt_col]]), ref = ref_label)
  d
}

# --- 2) δ grid ---------------------------------------------------------------
delta_grid <- seq(-2, 2, by = 0.25)

# --- 3) Models & helpers -----------------------------------------------------
# MAR predictor (used ONLY to generate shifted probabilities for missing EPDS_12)
form_pred <- reformulate(
  c(trt_col, "v01bas_age_2","v01bas_canread_4","v01bas_canwrite_5",
    "v01bas_livingchd_cat","v01bas_weight_14","v01bas_gestage_17"),
  response = "EPDS_12"
)

# Unadjusted analysis model (this is what gets pooled and plotted)
form_bin_unadj <- reformulate(trt_col, response = "EPDS_12")

# Apply a common δ shift to imputed outcomes in NON-reference arms
adjust_imputed_binary <- function(d, delta) {
  d <- set_ref(d)
  d_obs <- d[!d$was_missing_epds12, , drop = FALSE]
  fit_pred <- suppressWarnings(glm(form_pred, data = d_obs, family = binomial()))
  eta_hat  <- suppressWarnings(predict(fit_pred, newdata = d, type = "link"))
  p_del    <- plogis(eta_hat + delta)
  p_del    <- pmin(pmax(p_del, 1e-8), 1 - 1e-8)
  idx <- which(d$was_missing_epds12 & d[[trt_col]] != ref_label & is.finite(p_del) & !is.na(p_del))
  if (length(idx)) d$EPDS_12[idx] <- rbinom(length(idx), 1, p_del[idx])
  d
}

extract_contrasts <- function(fit) {
  broom::tidy(fit) %>%
    dplyr::filter(grepl(paste0("^", trt_col), term)) %>%
    dplyr::transmute(
      contrast = sub(paste0("^", trt_col), "", term),
      est_vec  = estimate,
      se_vec   = std.error
    )
}

rubin_pool <- function(ests) {
  ests %>%
    dplyr::summarise(
      m_eff = dplyr::n(),
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
  ests <- purrr::map_dfr(seq_len(m), function(i) {
    d_i <- complete(imp, i)
    # restore the original missingness flag (from the raw data)
    d_i$was_missing_epds12 <- dat$was_missing_epds12
    d_i <- adjust_imputed_binary(d_i, delta)
    fit <- suppressWarnings(glm(form_bin_unadj, data = set_ref(d_i), family = binomial()))
    out <- extract_contrasts(fit); out$.imp <- i; out
  }) %>% dplyr::filter(is.finite(est_vec), is.finite(se_vec))
  if (!nrow(ests)) return(tibble())
  ests %>% dplyr::group_by(contrast) %>%
    rubin_pool() %>%
    dplyr::transmute(contrast, m_eff, estimate = Qbar, W, B, Tvar, se, lwr, upr, pval, delta = delta)
}

tp_bin_results <- purrr::map_dfr(delta_grid, run_one_delta) %>% dplyr::arrange(contrast, delta)

# --- 4) Pairwise Full vs Split (shift one arm at a time) ---------------------
pair_levels <- c("Mefloquine full dose", "Mefloquine split dose")

adjust_imputed_binary_perarm <- function(d, delta_map) {
  d <- set_ref(d)
  d_obs <- d[!d$was_missing_epds12, , drop = FALSE]
  fit   <- suppressWarnings(glm(form_pred, data = d_obs, family = binomial()))
  eta0  <- suppressWarnings(predict(fit, newdata = d, type = "link"))
  p0    <- plogis(eta0); p_del <- p0
  for (arm in names(delta_map)) {
    idx <- d$was_missing_epds12 & d[[trt_col]] == arm
    if (any(idx)) p_del[idx] <- plogis(qlogis(p0[idx]) + delta_map[[arm]])
  }
  p_del <- pmin(pmax(p_del, 1e-8), 1 - 1e-8)
  idx_all <- which(d$was_missing_epds12 & d[[trt_col]] %in% names(delta_map) & !is.na(p_del))
  if (length(idx_all)) d$EPDS_12[idx_all] <- rbinom(length(idx_all), 1, p_del[idx_all])
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
  ests <- purrr::map_dfr(seq_len(m), function(i) {
    d_i <- complete(imp, i)
    d_i$was_missing_epds12 <- dat$was_missing_epds12
    d_i <- adjust_imputed_binary_perarm(d_i, delta_map = setNames(delta, shift_arm))
    fit <- suppressWarnings(glm(form_bin_unadj, data = set_ref(d_i), family = binomial()))
    out <- extract_pair_contrast(fit); out$.imp <- i; out
  }) %>% dplyr::filter(is.finite(est_vec), is.finite(se_vec))
  if (!nrow(ests)) return(tibble())
  ests %>% rubin_pool() %>%
    dplyr::transmute(contrast = paste(pair_levels, collapse = " vs "),
                     shifted_arm = shift_arm,
                     m_eff, estimate = Qbar, W, B, Tvar, se, lwr, upr, pval, delta)
}

delta_grid_pair <- seq(-2, 2, by = 0.25)
tp_pair_results_full  <- purrr::map_dfr(delta_grid_pair, ~ run_one_delta_pair(.x, shift_arm = "Mefloquine full dose"))
tp_pair_results_split <- purrr::map_dfr(delta_grid_pair, ~ run_one_delta_pair(.x, shift_arm = "Mefloquine split dose"))

# --- 5) Build one long DF with 4 panels --------------------------------------
df_full_vs_sp  <- tp_bin_results %>%
  dplyr::filter(contrast == "Mefloquine full dose") %>%
  dplyr::mutate(panel = "Full vs SP (δ on MQ imputed)")

df_split_vs_sp <- tp_bin_results %>%
  dplyr::filter(contrast == "Mefloquine split dose") %>%
  dplyr::mutate(panel = "Split vs SP (δ on MQ imputed)")

df_pair_full_shift  <- tp_pair_results_full  %>% dplyr::mutate(panel = "Full vs Split (shift Full)")
df_pair_split_shift <- tp_pair_results_split %>% dplyr::mutate(panel = "Full vs Split (shift Split)")

df_all <- dplyr::bind_rows(
  df_full_vs_sp  %>% dplyr::select(delta, estimate, lwr, upr, pval, panel),
  df_split_vs_sp %>% dplyr::select(delta, estimate, lwr, upr, pval, panel),
  df_pair_full_shift  %>% dplyr::select(delta, estimate, lwr, upr, pval, panel),
  df_pair_split_shift %>% dplyr::select(delta, estimate, lwr, upr, pval, panel)
) %>%
  dplyr::mutate(
    sig   = pval < 0.05,
    panel = factor(panel, levels = c(
      "Full vs SP (δ on MQ imputed)",
      "Split vs SP (δ on MQ imputed)",
      "Full vs Split (shift Full)",
      "Full vs Split (shift Split)"
    ))
  )

# --- 6) Tipping points per panel (first δ where 95% CI excludes 0) -----------
tip_one <- function(d) {
  nm <- names(d)
  delta_col <- intersect(c("delta","delta_shift","d","delta_v","shift","logodds_shift"), nm)[1]
  lo_col    <- intersect(c("lcl","ci_low","lo","lower","ymin","lwr","low","cil","ci_l"), nm)[1]
  hi_col    <- intersect(c("ucl","ci_high","hi","upper","ymax","upr","high","ciu","ci_u"), nm)[1]
  stopifnot(!is.na(delta_col), !is.na(lo_col), !is.na(hi_col))
  d <- d %>%
    dplyr::mutate(..delta = as.numeric(.data[[delta_col]]),
                  ..lo    = as.numeric(.data[[lo_col]]),
                  ..hi    = as.numeric(.data[[hi_col]])) %>%
    dplyr::arrange(..delta)
  excl0 <- (d$..lo > 0) | (d$..hi < 0)
  pos <- d %>% dplyr::filter(..delta >= 0, excl0) %>%
    dplyr::summarize(delta_pos = if (dplyr::n() > 0) min(..delta) else NA_real_) %>% dplyr::pull(delta_pos)
  neg <- d %>% dplyr::filter(..delta <= 0, excl0) %>%
    dplyr::summarize(delta_neg = if (dplyr::n() > 0) max(..delta) else NA_real_) %>% dplyr::pull(delta_neg)
  tibble::tibble(delta_pos = pos, delta_neg = neg)
}

tips_4 <- df_all %>%
  dplyr::group_by(panel) %>%
  dplyr::group_modify(~ tip_one(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    first_tip_delta = dplyr::case_when(
      !is.na(delta_pos) & !is.na(delta_neg) ~ ifelse(abs(delta_pos) <= abs(delta_neg), delta_pos, delta_neg),
      !is.na(delta_pos) ~ delta_pos,
      !is.na(delta_neg) ~ delta_neg,
      TRUE ~ NA_real_
    ),
    first_tip_dir = dplyr::case_when(
      is.na(first_tip_delta) ~ "No tipping in range",
      first_tip_delta > 0    ~ "CI > 0",
      TRUE                   ~ "CI < 0"
    )
  )

print(tips_4 %>% dplyr::mutate(dplyr::across(c(delta_pos, delta_neg, first_tip_delta), ~ round(., 2))))

tp_vlines4 <- tips_4 %>%
  tidyr::pivot_longer(c(delta_pos, delta_neg), names_to = "type", values_to = "delta_v") %>%
  dplyr::filter(!is.na(delta_v))

tips_4_long <- tips_4 %>%
  dplyr::select(panel, delta_pos, delta_neg) %>%
  tidyr::pivot_longer(c(delta_pos, delta_neg),
                      names_to = "direction_raw", values_to = "delta_star") %>%
  dplyr::filter(!is.na(delta_star)) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      direction_raw == "delta_pos" ~ "Positive δ (higher odds of depression)",
      direction_raw == "delta_neg" ~ "Negative δ (lower odds of depression)"
    )
  ) %>%
  dplyr::select(panel, direction, delta_star)

# --- 7) Single ggplot with 4 panels ------------------------------------------
y_tops <- df_all %>%
  dplyr::group_by(panel) %>%
  dplyr::summarise(y_top = max(upr, na.rm = TRUE), .groups = "drop")

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
    title    = "Tipping-point sensitivity to missing outcomes — EPDS ≥ 12 (binary)",
    subtitle = "δ represents a log-odds (≈ SD-unit) shift applied to imputed outcomes; models are unadjusted.",
    x = expression(delta~"(log-odds shift ≈ SD-unit)"),
    y = "Pooled log-odds ratio",
    shape = "p < 0.05",
    caption = "Dashed vertical lines mark the first δ (positive and negative) where the 95% CI excludes 0."
  )

print(p4)
ggplot2::ggsave("tp_epds12_four_panel_faceted.png", p4, width = 11, height = 9, dpi = 300)

# --- 8) Pooled results at δ = 0 (primary, unadjusted) -------------------------
mi_unadj <- purrr::map_dfr(seq_len(m), function(i) {
  d_i <- complete(imp, i) |> set_ref()
  fit <- suppressWarnings(glm(form_bin_unadj, data = d_i, family = binomial()))
  extract_contrasts(fit) |> dplyr::mutate(.imp = i)
}) |>
  dplyr::group_by(contrast) |>
  rubin_pool() |>
  dplyr::transmute(
    contrast,
    estimate = Qbar, se, lwr, upr, pval,
    model = "Unadjusted (δ = 0)"
  )

# View (log-odds) and as ORs
mi_unadj %>% dplyr::mutate(dplyr::across(c(estimate, se, lwr, upr, pval), ~ round(., 3))) %>% print(n = Inf)

mi_unadj_or <- mi_unadj %>%
  dplyr::transmute(
    contrast, 
    OR = exp(estimate),
    CI_low = exp(lwr),
    CI_high = exp(upr),
    pval
  )

mi_unadj_or %>%
  dplyr::mutate(dplyr::across(c(OR, CI_low, CI_high, pval), ~ round(., 3))) %>%
  print(n = Inf)



# interpretation of the delta shifts (percentage points added or subtracted)

# --- Make a 'long' tips table with both directions ----------------------------
tips_4_long <- tips_4 %>%
  dplyr::select(panel, delta_pos, delta_neg) %>%
  tidyr::pivot_longer(c(delta_pos, delta_neg),
                      names_to = "direction_raw", values_to = "delta_star") %>%
  dplyr::filter(!is.na(delta_star)) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      direction_raw == "delta_pos" ~ "Positive δ (higher odds of depression)",
      direction_raw == "delta_neg" ~ "Negative δ (lower odds of depression)"
    )
  ) %>%
  dplyr::select(panel, direction, delta_star)

get_p0_EPDS12_from_imp <- function(imp_obj, form_pred, trt_col, ref_label) {
  # NB: despite the name, this works for EPDS_12 too; it reads the outcome from form_pred
  outcome_var <- all.vars(form_pred)[1]
  m <- imp_obj$m
  n <- nrow(mice::complete(imp_obj, 1))
  pmat <- matrix(NA_real_, nrow = n, ncol = m)
  for (i in seq_len(m)) {
    d_i <- mice::complete(imp_obj, i)
    d_i[[trt_col]] <- stats::relevel(factor(d_i[[trt_col]]), ref = ref_label)
    obs_idx <- which(!is.na(d_i[[outcome_var]]))
    fit_i <- suppressWarnings(glm(form_pred, data = d_i[obs_idx, , drop = FALSE], family = binomial()))
    eta_i <- suppressWarnings(predict(fit_i, newdata = d_i, type = "link"))
    pmat[, i] <- plogis(eta_i)
  }
  list(p0 = rowMeans(pmat, na.rm = TRUE))
}

summarize_missing_arm_all <- function(d, p0_list, arm, delta_star) {
  # uses 'was_missing_epds12' and 'EPDS_12' present in your EPDS_12 dataset
  d2 <- set_ref(d)
  p0 <- p0_list$p0
  idx_miss <- d2$was_missing_epds12 & d2[[trt_col]] == arm
  M <- sum(idx_miss)
  if (M == 0) {
    return(dplyr::tibble(
      arm = arm, M_missing = 0,
      p0_mean = NA_real_, p_tip_mean = NA_real_,
      pct_needed = NA_real_, extra_cases = NA_real_
    ))
  }
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

tips_4_long <- tips_4 %>%
  dplyr::select(panel, delta_pos, delta_neg) %>%
  tidyr::pivot_longer(c(delta_pos, delta_neg),
                      names_to = "direction_raw", values_to = "delta_star") %>%
  dplyr::filter(!is.na(delta_star)) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      direction_raw == "delta_pos" ~ "Positive δ (higher odds of depression)",
      direction_raw == "delta_neg" ~ "Negative δ (lower odds of depression)"
    )
  ) %>%
  dplyr::select(panel, direction, delta_star)


# --- Compute interpretability rows for BOTH directions per panel --------------
calc_tip_requirements_all <- function(original_dat, tips_tbl_long) {
  d0 <- set_ref(original_dat)
  # Generic p0 getter (reads outcome name from 'form_pred')
  p0_list <- get_p0_EPDS12_from_imp(imp, form_pred, trt_col, ref_label)
  
  levs <- levels(d0[[trt_col]])
  nonref_arms <- setdiff(levs, ref_label)
  full_arm    <- "Mefloquine full dose"
  split_arm   <- "Mefloquine split dose"
  
  purrr::map_dfr(seq_len(nrow(tips_tbl_long)), function(i) {
    panel      <- tips_tbl_long$panel[i]
    direction  <- tips_tbl_long$direction[i]
    delta_star <- tips_tbl_long$delta_star[i]
    if (is.na(delta_star)) return(NULL)
    
    arms_to_shift <-
      if (panel %in% c("Full vs SP (δ on MQ imputed)", "Split vs SP (δ on MQ imputed)")) {
        nonref_arms
      } else if (panel == "Full vs Split (shift Full)") {
        full_arm
      } else if (panel == "Full vs Split (shift Split)") {
        split_arm
      } else character(0)
    
    purrr::map_dfr(arms_to_shift, ~ summarize_missing_arm_all(d0, p0_list, .x, delta_star)) %>%
      dplyr::mutate(
        panel = panel,
        direction = direction,
        delta_star = delta_star,
        odds_multiplier = exp(delta_star)
      ) %>%
      dplyr::select(panel, direction, delta_star, odds_multiplier, arm, M_missing,
                    p0_mean, p_tip_mean, pct_needed, extra_cases)
  })
}


tip_requirements_12 <- calc_tip_requirements_all(dat, tips_4_long)

# Pretty, two-direction table (like the others)
tip_requirements_epds12 <- tip_requirements_12 %>%
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
    `Mean p(δ*)`  = p_tip_mean,
    `% Depressed (δ*)` = pct_needed,
    `Δ Cases (δ* - MAR)` = extra_cases,
    `N Missing` = M_missing
  )

print(tip_requirements_epds12, n = Inf)

library(openxlsx)
write.xlsx(
  tip_requirements_epds12,
  "epds12_tipping_point_results_table.xlsx",
  rowNames = FALSE
)



# --- Pretty narrative lines for EPDS ≥ 12 tipping points ----------------------

# (Optional) nicer panel labels so the lines read well in text
panel_nice <- c(
  "Full vs SP (δ on MQ imputed)"  = "Full vs SP (δ on MQ imputed: Full & Split shifted)",
  "Split vs SP (δ on MQ imputed)" = "Split vs SP (δ on MQ imputed: Full & Split shifted)",
  "Full vs Split (shift Full)"    = "Full vs Split (shift Full)",
  "Full vs Split (shift Split)"   = "Full vs Split (shift Split)"
)

lines_tbl_epds12 <- tip_requirements_12 %>%
  dplyr::mutate(
    panel_lbl = dplyr::recode(panel, !!!panel_nice),
    # round for readability in text, keep underlying table unrounded
    delta_star_r = round(delta_star, 2),
    odds_mult_r  = round(odds_multiplier, 2),
    p0_pct_r     = round(100 * p0_mean, 1),
    p_tip_pct_r  = round(pct_needed, 1),
    extra_cases_r= round(extra_cases, 1)
  ) %>%
  dplyr::arrange(panel_lbl, direction, arm)

# Print one-liners to console
lines_epds12 <- sprintf(
  "%s – %s (%s): at the tip (δ = %.2f; odds × %.2f), ~%.1f%% of the missing would be depressed (vs %.1f%% MAR); ~%+.1f cases among %d missing.",
  lines_tbl_epds12$panel_lbl,
  lines_tbl_epds12$arm,
  lines_tbl_epds12$direction,
  lines_tbl_epds12$delta_star_r,
  lines_tbl_epds12$odds_mult_r,
  lines_tbl_epds12$p_tip_pct_r,
  lines_tbl_epds12$p0_pct_r,
  lines_tbl_epds12$extra_cases_r,
  lines_tbl_epds12$M_missing
)

cat(lines_epds12, sep = "\n")



