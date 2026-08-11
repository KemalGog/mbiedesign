# Summary outcomes and sample size by trial design: the main user-facing
# function, layering the sample-size and relative-efficiency calculation on top
# of the natural-history projection.

#' Summary outcomes and sample size by trial design
#'
#' The main function of the package. Given a calibrated natural-history
#' `rate_matrix` and a screening schedule, it projects the end-of-trial
#' late-stage outcomes and computes the required per-arm sample size and
#' relative efficiency for the traditional and intended-effect (IE) designs,
#' returning one row per design (a summary table). For the full
#' time-dependent early/late/overall outcomes, use
#' [get_trial_results_by_design()].
#'
#' @param rate_matrix Natural-history transition rate matrix (e.g. from
#'   [load_fitted_model()] or [fit_natural_history()]).
#' @param start_age Age at trial entry.
#' @param numscreens Number of screening examinations.
#' @param screen_int Interval between screens (years).
#' @param num_followup_intervals Follow-up after the final screen.
#' @param sens_e Early-stage test sensitivity.
#' @param sens_l Late-stage test sensitivity.
#' @param specificity Test specificity.
#' @param design Which design(s) to return: `"traditional"`, `"IE"`, or both.
#'   Defaults to both.
#' @param allocation_ratio Screen-arm to control-arm allocation ratio
#'   (n_screen / n_control); 1 is balanced (the default). Event rates are
#'   per-participant, so only the allocation ratio affects the required sample
#'   size — absolute arm sizes are not needed to estimate it.
#' @param alpha One-sided significance level.
#' @param power Target power (= 1 - beta).
#'
#' @return A data frame with one row per requested `design` and columns:
#'   \describe{
#'     \item{design}{`"traditional"` or `"IE"`.}
#'     \item{p_c, p_s}{Control- and screen-arm late-stage event probabilities
#'       (per randomized participant for traditional; per ever-positive for IE).}
#'     \item{rr, rd}{Risk ratio (`p_s / p_c`) and risk difference (`p_c - p_s`).}
#'     \item{fraction}{The subgroup fraction the design conditions on: `NA` for
#'       traditional, the ever-positive fraction (pi) for IE.}
#'     \item{n_per_arm}{Required per-arm sample size.}
#'     \item{relative_efficiency}{N_traditional / N_design (1 for traditional).}
#'   }
#' @export
get_summary_outcomes_by_design <- function(rate_matrix,
                         start_age, numscreens, screen_int,
                         num_followup_intervals,
                         sens_e, sens_l, specificity,
                         design = c("traditional", "IE"),
                         allocation_ratio = 1,
                         alpha = 0.025, power = 0.90) {

  design <- match.arg(design, c("traditional", "IE"), several.ok = TRUE)

  # Event rates are per-participant (independent of arm sizes), so unit arm
  # sizes are used here; the allocation ratio enters only the sample-size
  # formula below.
  res <- get_trial_results_by_design(
    start_age = start_age, numscreens = numscreens, screen_int = screen_int,
    num_followup_intervals = num_followup_intervals, rate_matrix = rate_matrix,
    sens_e = sens_e, sens_l = sens_l, specificity = specificity,
    n_control = 1, n_screen = 1
  )
  tl <- function(x) as.numeric(tail(x, 1))

  # Traditional
  p_c_std <- tl(res$control$standard$late$incidence$cumulative$control_arm_detected)
  p_s_std <- tl(res$screen$standard$late$incidence$cumulative$screen_arm_detected)
  # Intended-effect
  p_c_ie  <- tl(res$control$IE$late$incidence$cumulative$control_arm_detected)
  p_s_ie  <- tl(res$screen$IE$late$incidence$cumulative$screen_arm_detected)
  pi_pos  <- get_positivity_rate(
    sens_e = sens_e, sens_l = sens_l, specificity = specificity,
    screen_times = res$params$screen_times, start_age = start_age,
    rate_matrix = rate_matrix)

  R     <- allocation_ratio
  n_std <- n_required_by_design(power, p_c_std, p_s_std, "traditional", alpha, R)
  n_ie  <- n_required_by_design(power, p_c_ie,  p_s_ie,  "IE", alpha, R, pi_pos = pi_pos)

  out <- data.frame(
    design              = c("traditional", "IE"),
    p_c                 = c(p_c_std, p_c_ie),
    p_s                 = c(p_s_std, p_s_ie),
    rr                  = c(p_s_std / p_c_std, p_s_ie / p_c_ie),
    rd                  = c(p_c_std - p_s_std, p_c_ie - p_s_ie),
    fraction            = c(NA_real_, pi_pos),
    n_per_arm           = c(n_std, n_ie) / (1 + R),
    relative_efficiency = n_std / c(n_std, n_ie),
    stringsAsFactors    = FALSE
  )
  out[out$design %in% design, , drop = FALSE]
}

#' Relative efficiency of the intended-effect vs traditional design
#'
#' Ratio of the sample size required by the traditional design to that required
#' by the intended-effect (IE) design for the same event probabilities. Values
#' greater than 1 mean the IE design needs fewer participants.
#'
#' @param p_c_std,p_s_std Traditional control- and screen-arm event probabilities.
#' @param p_c_ie,p_s_ie IE control- and screen-arm event probabilities.
#' @param pi_pos Ever-positive fraction.
#' @param alpha One-sided significance level.
#' @param power Target power.
#' @param R Allocation ratio, screen / control.
#'
#' @return The relative efficiency (a single numeric value).
#' @export
relative_efficiency <- function(p_c_std, p_s_std, p_c_ie, p_s_ie, pi_pos,
                                alpha = 0.025, power = 0.90, R = 1) {
  n_std <- n_required_by_design(power, p_c_std, p_s_std, "traditional", alpha, R)
  n_ie  <- n_required_by_design(power, p_c_ie, p_s_ie, "IE", alpha, R, pi_pos = pi_pos)
  n_std / n_ie
}

# Look up the matched overall mean sojourn time (OMST) for a target early-stage
# mean sojourn time (EMST) at a given late-stage mean sojourn time (LMST), from
# a matched-pairs table (see load_fitted_model()$matched_pairs). Internal.
match_omst <- function(matched_pairs, target_emst, lmst, tol = 0.01) {
  d <- abs(matched_pairs$target_EMST - target_emst) +
       abs(matched_pairs$LMST - lmst)
  ok <- abs(matched_pairs$target_EMST - target_emst) <= tol &
        abs(matched_pairs$LMST - lmst) <= tol
  if (!any(ok)) {
    stop("No matched pair found for target EMST = ", target_emst,
         " and LMST = ", lmst, ".")
  }
  matched_pairs$matched_OMST[which.min(ifelse(ok, d, Inf))]
}
