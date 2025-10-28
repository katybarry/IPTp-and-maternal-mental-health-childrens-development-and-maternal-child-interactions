# --- Libraries -----------------------------------------------------------------
library(mice)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

table(tipping_maternal_db$medgroup)

# --- Adjusted MI for HOME (TOV_hometot, continuous) ---------------------------

dat <- tipping_maternal_db
trt_col   <- "medgroup"
ref_label <- "Sulfadoxine-pryimethamine dose"  # exact spelling

# Keep ONLY variables used in the MI and analysis
dat <- dat %>%
  dplyr::select(
    TOV_hometot, !!trt_col,
    v01bas_age_2, v01bas_canread_4, v01bas_canwrite_5,
    v01bas_livingchd_cat, v01bas_weight_14, v01bas_gestage_17
  )

# Types
dat[[trt_col]]            <- as.factor(dat[[trt_col]])
dat$v01bas_canread_4      <- as.factor(dat$v01bas_canread_4)
dat$v01bas_canwrite_5     <- as.factor(dat$v01bas_canwrite_5)
dat$v01bas_livingchd_cat  <- ordered(dat$v01bas_livingchd_cat)  # ordered -> polr

# Analysis frame alias 
dat_anl <- dat %>%
  dplyr::mutate(!!trt_col := stats::relevel(factor(.data[[trt_col]]), ref = ref_label))

# Flag who was missing the outcome BEFORE imputation (used for δ-shift later)
was_missing <- is.na(dat_anl$TOV_hometot)

# --- Multiple imputation (MAR for TOV_hometot + impute missing predictors) ----
m <- 20
set.seed(817)

tmp  <- mice(dat, m = 1, maxit = 0)
meth <- tmp$method
pred <- tmp$predictorMatrix

meth[] <- ""

# Methods
meth["TOV_hometot"]          <- "pmm"   # continuous outcome
meth["v01bas_age_2"]         <- "pmm"   # continuous predictors (only used if NA)
meth["v01bas_weight_14"]     <- "pmm"
meth["v01bas_gestage_17"]    <- "pmm"
meth["v01bas_canread_4"]     <- "logreg"  # binary factors
meth["v01bas_canwrite_5"]    <- "logreg"
meth["v01bas_livingchd_cat"] <- "polr"    # ordered factor
meth[trt_col]                <- ""        # do NOT impute treatment

# Predictor matrix
diag(pred) <- 0                  # standard
pred[trt_col, ] <- 0             # don't impute treatment

# Predictors used for TOV_hometot 
pred["TOV_hometot", ] <- 0
pred["TOV_hometot", c(
  trt_col, "v01bas_age_2","v01bas_canread_4","v01bas_canwrite_5",
  "v01bas_livingchd_cat","v01bas_weight_14","v01bas_gestage_17"
)] <- 1

# Predictors used for v01bas_livingchd_cat 
pred["v01bas_livingchd_cat", ] <- 0
pred["v01bas_livingchd_cat", c(
  trt_col, "v01bas_age_2","v01bas_canread_4","v01bas_canwrite_5",
  "v01bas_weight_14","v01bas_gestage_17","TOV_hometot"
)] <- 1

# Restrict imputation to cells with NA only
where <- matrix(FALSE, nrow = nrow(dat), ncol = ncol(dat),
                dimnames = list(NULL, names(dat)))
for (nm in names(dat)) where[, nm] <- is.na(dat[[nm]])

# Run mice
imp <- mice(dat, m = m, method = meth, predictorMatrix = pred, where = where,
            maxit = 20, seed = 817)

# --- 2) SD scaling and δ-grid --------------------------------------------------
ref_label <- "Sulfadoxine-pryimethamine dose"
dat_anl$medgroup <- stats::relevel(factor(dat_anl$medgroup), ref = ref_label)

sd_nonref_obs <- dat_anl %>%
  dplyr::filter(medgroup != ref_label) %>%
  dplyr::summarise(sd = sd(TOV_hometot, na.rm = TRUE)) %>%
  dplyr::pull(sd)
if (!is.finite(sd_nonref_obs) || length(sd_nonref_obs) == 0) {
  sd_nonref_obs <- sd(dat_anl$TOV_hometot, na.rm = TRUE)
}
cat(sprintf("Scaling SD for TOV_hometot (non-ref observed): %.3f\n", sd_nonref_obs))

# Choose grid using both SD and a max-SE heuristic
max_se_home <- 2.39
range_by_se <- 4 * max_se_home
range_by_sd <- 3 * sd_nonref_obs
K_raw <- max(range_by_se, range_by_sd)
# --- Use exactly ±2 SD for the δ grid ---------------------------------------
# (Drop the SE-based heuristic)
K <- 2 * sd_nonref_obs    

# Step ~ SD/12, clamped to [0.5, 2] points for a smooth curve/runtime balance
step <- sd_nonref_obs / 12
step <- min(max(step, 0.5), 2)

delta_grid <- seq(-K, K, by = step)

cat(sprintf("δ grid set to ±2 SD: ±%.2f points (step=%.2f)  [SD≈%.2f]\n",
            K, step, sd_nonref_obs))


# --- 3) Contrast + Pooling helpers --------------------------------------------
contrast_stat <- function(fit, A, B, var = "medgroup") {
  levs <- levels(fit$model[[var]])
  ref  <- levs[1]
  cf   <- coef(fit); V <- vcov(fit); rn <- names(cf)
  cn   <- function(level) paste0(var, level)
  L <- setNames(rep(0, length(cf)), rn)
  if (A == ref && B != ref) {
    nb <- cn(B); if (!nb %in% rn) return(NULL); L[nb] <- -1
  } else if (B == ref && A != ref) {
    na <- cn(A); if (!na %in% rn) return(NULL); L[na] <- 1
  } else if (A != ref && B != ref) {
    na <- cn(A); nb <- cn(B); if (!na %in% rn || !nb %in% rn) return(NULL); L[na] <- 1; L[nb] <- -1
  } else return(NULL)
  est <- as.numeric(crossprod(L, cf))
  var <- as.numeric(t(L) %*% V %*% L)
  data.frame(contrast = sprintf("%s vs %s", A, B), est = est, var = var)
}
pool_scalar <- function(ests, vars, dfcom) {
  Qbar <- mean(ests); Ubar <- mean(vars); B <- var(ests); m <- length(ests)
  Tvar <- Ubar + (1 + 1/m) * B
  r    <- (1 + 1/m) * B / Ubar
  df   <- (m - 1) * (1 + 1/r)^2; if (!is.finite(df)) df <- dfcom
  se   <- sqrt(Tvar); tcr <- qt(0.975, df)
  lwr  <- Qbar - tcr * se; upr <- Qbar + tcr * se
  pval <- 2 * pt(-abs(Qbar / se), df)
  data.frame(estimate = Qbar, lwr = lwr, upr = upr, pval = pval, df = df)
}

# --- 4) δ-adjusted MI functions -----------------------------------------------
run_tipping_mi <- function(delta) {
  rows <- lapply(seq_len(m), function(i) {
    d <- mice::complete(imp, i)
    d$medgroup <- stats::relevel(factor(d$medgroup), ref = ref_label)
    d$TOV_hometot_delta <- d$TOV_hometot
    idx <- was_missing & d$medgroup != ref_label
    if (any(idx)) d$TOV_hometot_delta[idx] <- d$TOV_hometot_delta[idx] + delta
    fit <- lm(TOV_hometot_delta ~ medgroup, data = d)
    levs <- levels(d$medgroup)
    do.call(rbind, lapply(setdiff(levs, ref_label),
                          function(B) contrast_stat(fit, A = B, B = ref_label)))
  })
  long <- dplyr::bind_rows(rows)
  long %>%
    dplyr::group_by(contrast) %>%
    dplyr::summarise(
      estimate = pool_scalar(est, var, dfcom = nrow(dat_anl) - nlevels(dat_anl$medgroup))$estimate,
      lwr      = pool_scalar(est, var, dfcom = nrow(dat_anl) - nlevels(dat_anl$medgroup))$lwr,
      upr      = pool_scalar(est, var, dfcom = nrow(dat_anl) - nlevels(dat_anl$medgroup))$upr,
      pval     = pool_scalar(est, var, dfcom = nrow(dat_anl) - nlevels(dat_anl$medgroup))$pval,
      .groups = "drop"
    ) %>% dplyr::mutate(delta = delta)
}

tp_cont_sp <- dplyr::bind_rows(lapply(delta_grid, run_tipping_mi)) %>%
  dplyr::arrange(contrast, delta)

# Pairwise (Full vs Split): shift one arm at a time
pair_levels <- c("Mefloquine full dose", "Mefloquine split dose")
adjust_imputed_cont_perarm <- function(d, delta, shift_arm) {
  d$medgroup <- stats::relevel(factor(d$medgroup), ref = ref_label)
  d$TOV_hometot_delta <- d$TOV_hometot
  idx <- was_missing & d$medgroup == shift_arm
  if (any(idx)) d$TOV_hometot_delta[idx] <- d$TOV_hometot_delta[idx] + delta
  d
}
extract_pair_contrast_lm <- function(fit) {
  contrast_stat(fit, A = pair_levels[1], B = pair_levels[2], var = "medgroup")
}
run_one_delta_pair_cont <- function(delta, shift_arm) {
  rows <- lapply(seq_len(m), function(i) {
    d_i <- mice::complete(imp, i)
    d_i <- adjust_imputed_cont_perarm(d_i, delta, shift_arm)
    fit <- lm(TOV_hometot_delta ~ medgroup, data = d_i)
    extract_pair_contrast_lm(fit)
  })
  long <- dplyr::bind_rows(rows)
  pooled <- pool_scalar(long$est, long$var, dfcom = nrow(dat_anl) - nlevels(dat_anl$medgroup))
  dplyr::tibble(
    contrast = paste(pair_levels, collapse = " vs "),
    shifted_arm = shift_arm,
    estimate = pooled$estimate, lwr = pooled$lwr, upr = pooled$upr, pval = pooled$pval,
    delta = delta
  )
}
tp_pair_full  <- dplyr::bind_rows(lapply(delta_grid, run_one_delta_pair_cont, shift_arm = pair_levels[1]))
tp_pair_split <- dplyr::bind_rows(lapply(delta_grid, run_one_delta_pair_cont, shift_arm = pair_levels[2]))

# --- 5) Build 4 panels ---------------------------------------------------------

# --- Build the 4 panels into one long data frame ------------------------------
# Top row: vs SP (from tp_cont_sp)
df_full_vs_sp <- tp_cont_sp %>%
  dplyr::filter(contrast == "Mefloquine full dose vs Sulfadoxine-pryimethamine dose") %>%
  dplyr::mutate(panel = "Full vs SP (δ on MQ imputed)")

df_split_vs_sp <- tp_cont_sp %>%
  dplyr::filter(contrast == "Mefloquine split dose vs Sulfadoxine-pryimethamine dose") %>%
  dplyr::mutate(panel = "Split vs SP (δ on MQ imputed)")

# Bottom row: pairwise (from tp_pair_full / tp_pair_split)
df_pair_full_shift  <- tp_pair_full  %>% dplyr::mutate(panel = "Full vs Split (shift Full)")
df_pair_split_shift <- tp_pair_split %>% dplyr::mutate(panel = "Full vs Split (shift Split)")

# Bind and keep the columns the plotting code expects
df_all <- dplyr::bind_rows(
  df_full_vs_sp  %>% dplyr::select(delta, estimate, lwr, upr, pval, panel),
  df_split_vs_sp %>% dplyr::select(delta, estimate, lwr, upr, pval, panel),
  df_pair_full_shift  %>% dplyr::select(delta, estimate, lwr, upr, pval, panel),
  df_pair_split_shift %>% dplyr::select(delta, estimate, lwr, upr, pval, panel)
) %>%
  dplyr::mutate(
    sig = pval < 0.05,
    panel = factor(panel, levels = c(
      "Full vs SP (δ on MQ imputed)",
      "Split vs SP (δ on MQ imputed)",
      "Full vs Split (shift Full)",
      "Full vs Split (shift Split)"
    ))
  )

# --- 6) Tipping points ------------------------------
tip_two_dirs <- function(d) {
  o <- d[order(d$delta), , drop = FALSE]
  excl0 <- (o$lwr > 0) | (o$upr < 0)  # CI fully excludes 0
  
  # first δ ≥ 0 where CI excludes 0
  pos <- o$delta[o$delta >= 0 & excl0]
  delta_pos <- if (length(pos)) min(pos) else NA_real_
  
  # last δ ≤ 0 where CI excludes 0
  neg <- o$delta[o$delta <= 0 & excl0]
  delta_neg <- if (length(neg)) max(neg) else NA_real_
  
  tibble::tibble(delta_pos = delta_pos, delta_neg = delta_neg)
}

tips_4 <- df_all %>%
  dplyr::group_by(panel) %>%
  dplyr::group_modify(~ tip_two_dirs(.x)) %>%
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

# Make a long table with one row per panel × direction
tips_4_long <- tips_4 %>%
  dplyr::select(panel, delta_pos, delta_neg) %>%
  tidyr::pivot_longer(c(delta_pos, delta_neg),
                      names_to = "direction_raw", values_to = "delta_star") %>%
  dplyr::filter(!is.na(delta_star)) %>%
  dplyr::mutate(
    direction = dplyr::case_when(
      direction_raw == "delta_pos" ~ "positive (CI > 0)",
      direction_raw == "delta_neg" ~ "negative (CI < 0)"
    )
  ) %>%
  dplyr::select(panel, direction, delta_star)



# --- 7) Build data for dashed lines + labels ----------------------------------
y_tops <- df_all %>%
  dplyr::group_by(panel) %>%
  dplyr::summarise(y_top = max(upr, na.rm = TRUE), .groups = "drop")

tp_vlines4 <- tips_4 %>%
  tidyr::pivot_longer(c(delta_pos, delta_neg), names_to = "type", values_to = "delta_v") %>%
  dplyr::filter(!is.na(delta_v)) %>%
  dplyr::left_join(y_tops, by = "panel") %>%
  dplyr::mutate(lbl = paste0("δ=", round(delta_v, 2)))

# --- 8) Axes + plot -----------------------------------------------------------
ylim_all <- range(c(df_all$lwr, df_all$upr), na.rm = TRUE)
K        <- max(abs(range(df_all$delta)))
sdx      <- sd_nonref_obs
sd_lim   <- floor((K / sdx) * 2) / 2
sd_breaks <- seq(-sd_lim, sd_lim, by = 0.5)

p4 <- ggplot(df_all, aes(x = delta, y = estimate)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15) +
  geom_line() +
  geom_point(aes(shape = sig), size = 2) +
  geom_vline(data = tp_vlines4, aes(xintercept = delta_v),
             linetype = "dashed", linewidth = 0.4, alpha = 0.7) +
  geom_text(data = tp_vlines4,
            aes(x = delta_v, y = y_top, label = lbl),
            vjust = -0.35, size = 3) +
  facet_wrap(~ panel, ncol = 2) +
  # ---- Zoom y-axis here (adjust to taste) ----
coord_cartesian(ylim = c(-3, 5)) +
  # --------------------------------------------
scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  scale_x_continuous(
    limits = c(-K, K),
    sec.axis = sec_axis(~ . / sdx,
                        breaks = sd_breaks,
                        labels = function(x) sprintf("%.1f", x),
                        name = expression(delta~"(in SD units)"))
  ) +
  theme_bw() +
  labs(
    title    = "Tipping-point sensitivity to missing outcomes — HOME total (continuous)",
    subtitle = "δ is a points-shift applied only to imputed values in non-reference arms. Unadjusted analysis model.",
    x = expression(delta~"(points)"),
    y = "Mean difference (MQ minus comparator) with 95% CI",
    shape = "p < 0.05",
    caption = "Dashed lines indicate δ where the 95% CI first excludes 0 (both directions)."
  )

print(p4)
ggsave("tp_home_four_panel.png", p4, width = 11, height = 9, dpi = 300)

# --- 8) Primary (δ=0) pooled estimates table ----------------------------------
panel_lvls <- c(
  "Full vs SP (δ on MQ imputed)",
  "Split vs SP (δ on MQ imputed)",
  "Full vs Split (shift Full)",
  "Full vs Split (shift Split)"
)
sp0 <- run_tipping_mi(0)
primary_sp <- sp0 %>%
  dplyr::transmute(
    panel = dplyr::case_when(
      contrast == "Mefloquine full dose vs Sulfadoxine-pryimethamine dose"  ~ panel_lvls[1],
      contrast == "Mefloquine split dose vs Sulfadoxine-pryimethamine dose" ~ panel_lvls[2],
      TRUE ~ NA_character_
    ),
    estimate, lwr, upr, pval
  ) %>% dplyr::filter(!is.na(panel))
pair0_full  <- run_one_delta_pair_cont(0, shift_arm = "Mefloquine full dose")
pair0_split <- run_one_delta_pair_cont(0, shift_arm = "Mefloquine split dose")
primary_pair <- dplyr::bind_rows(
  pair0_full  %>% dplyr::mutate(panel = panel_lvls[3]),
  pair0_split %>% dplyr::mutate(panel = panel_lvls[4])
) %>% dplyr::select(panel, estimate, lwr, upr, pval)
primary_all <- dplyr::bind_rows(primary_sp, primary_pair) %>%
  dplyr::mutate(panel = factor(panel, levels = panel_lvls)) %>%
  dplyr::arrange(panel) %>%
  dplyr::mutate(dplyr::across(c(estimate, lwr, upr, pval), ~ round(., 3)))
print(primary_all, n = Inf)

# === 9) INTERPRETABILITY for CONTINUOUS outcome ===============================
# Average score shift among the missing (points & SD) implied at the tipping δ*

set_ref <- function(d) { d[[trt_col]] <- stats::relevel(factor(d[[trt_col]]), ref = ref_label); d }

get_mu0_HOME_from_imp <- function(imp_obj) {
  m <- imp_obj$m
  n <- nrow(mice::complete(imp_obj, 1))
  mu_mat <- matrix(NA_real_, nrow = n, ncol = m)
  
  for (i in seq_len(m)) {
    d_i <- mice::complete(imp_obj, i)
    d_i <- set_ref(d_i)
    obs_i <- !is.na(d_i$TOV_hometot)
    
    # Arm means among observed outcomes, then assign each row its arm's mean
    m_by_arm <- tapply(d_i$TOV_hometot[obs_i], d_i$medgroup[obs_i], mean, na.rm = TRUE)
    mu_mat[, i] <- unname(m_by_arm[as.character(d_i$medgroup)])
    
  }
  list(mu0 = rowMeans(mu_mat, na.rm = TRUE))
}

summarize_missing_arm_cont <- function(d, mu_list, arm, delta_star, sd_units) {
  d2 <- set_ref(d)
  mu0 <- mu_list$mu0
  idx_miss <- was_missing & d2[[trt_col]] == arm
  M <- sum(idx_miss)
  if (M == 0) {
    return(dplyr::tibble(
      arm = arm, M_missing = 0,
      mu0_mean = NA_real_, mu_tip_mean = NA_real_,
      shift_points = NA_real_, shift_d = NA_real_
    ))
  }
  mu0_sub  <- mu0[idx_miss]
  mu_tip   <- mu0_sub + delta_star
  dplyr::tibble(
    arm          = arm,
    M_missing    = M,
    mu0_mean     = mean(mu0_sub, na.rm = TRUE),
    mu_tip_mean  = mean(mu_tip, na.rm = TRUE),
    shift_points = mean(mu_tip - mu0_sub, na.rm = TRUE),
    shift_d      = (mean(mu_tip - mu0_sub, na.rm = TRUE)) / sd_units
  )
}

calc_tip_requirements_home <- function(original_dat_anl, tips_tbl_long, sd_units) {
  d0 <- set_ref(original_dat_anl)
  mu_list <- get_mu0_HOME_from_imp(imp)
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
    
    purrr::map_dfr(arms_to_shift, ~ summarize_missing_arm_cont(d0, mu_list, .x, delta_star, sd_units)) %>%
      dplyr::mutate(
        panel = panel,
        direction = direction,
        delta_star_points = delta_star,
        delta_star_SD     = delta_star / sd_units
      ) %>%
      dplyr::select(panel, direction, arm, M_missing,
                    mu0_mean, mu_tip_mean,
                    shift_points, shift_d,
                    delta_star_points, delta_star_SD)
  })
}


home_requirements <- calc_tip_requirements_home(dat_anl, tips_4_long, sd_units = sd_nonref_obs)

panel_nice <- c(
  "Full vs SP (δ on MQ imputed)"  = "Full vs SP (δ on MQ imputed: Full & Split shifted)",
  "Split vs SP (δ on MQ imputed)" = "Split vs SP (δ on MQ imputed: Full & Split shifted)",
  "Full vs Split (shift Full)"    = "Full vs Split (shift Full)",
  "Full vs Split (shift Split)"   = "Full vs Split (shift Split)"
)

home_requirements %>%
  dplyr::mutate(
    panel_lbl = panel_nice[panel],
    line = sprintf(
      "%s – %s – %s dir: at the tip (δ = %.2f points; %.2f SD), the missing would average %.2f points (vs %.2f MAR), i.e., a shift of %.2f points (%.2f SD) across %d missing.",
      panel_lbl, arm, direction, delta_star_points, delta_star_SD,
      mu_tip_mean, mu0_mean, shift_points, shift_d, M_missing
    )
  ) %>%
  dplyr::pull(line) %>%
  cat(sep = "\n")

home_requirements %>%
  dplyr::mutate(dplyr::across(
    c(mu0_mean, mu_tip_mean, shift_points, shift_d, delta_star_points, delta_star_SD),
    ~ round(., 3)
  )) %>%
  dplyr::arrange(panel, direction, arm) %>%
  print(n = Inf)


library(openxlsx)
write.xlsx(home_requirements, "home_tipping_point_results_table.xlsx", rowNames = FALSE)
