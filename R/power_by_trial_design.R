# Power and sample-size calculations by trial design:
# traditional and intended-effect (IE).

#' Power for a given sample size, by design
#'
#' Computes the power to detect a difference in event probabilities for the
#' traditional or intended-effect (`"IE"`) design. The IE design inflates the
#' variance by `1 / pi_pos`, reflecting that inference is restricted to the
#' ever-positive subgroup.
#'
#' @param n_total Total sample size (both arms).
#' @param p_c,p_s Control- and screen-arm event probabilities (use the
#'   probabilities appropriate to the chosen design).
#' @param design One of `"traditional"` or `"IE"`.
#' @param alpha One-sided significance level.
#' @param R Allocation ratio, screen-arm size / control-arm size.
#' @param pi_pos Ever-positive fraction; required for `design = "IE"`.
#'
#' @return The power (a single numeric value).
#' @seealso [n_required_by_design()]
#' @export
power_by_design <- function(n_total,
                            p_c, p_s,
                            design = c("traditional", "IE"),
                            alpha = 0.025,
                            R = 1,
                            pi_pos = NULL) {

  design <- match.arg(design)

  n_c <- n_total / (1 + R)
  n_s <- n_total - n_c

  # pooled probability
  p <- (p_c + R * p_s) / (1 + R)

  z_alpha <- qnorm(1 - alpha)

  # variance scaling
  var_scale <- switch(design,
                      "traditional" = 1,
                      "IE" = {
                        if (is.null(pi_pos)) stop("pi_pos required for IE design")
                        1 / pi_pos
                      })

  se <- sqrt(var_scale * p * (1 - p) * (1 / n_c + 1 / n_s))
  z_effect <- (p_c - p_s) / se

  1 - pnorm(z_alpha - z_effect)
}

#' Required sample size for a target power, by design
#'
#' Computes the total sample size required to achieve a target power for the
#' traditional or intended-effect (`"IE"`) design. As in [power_by_design()],
#' the IE design inflates the variance by `1 / pi_pos`.
#'
#' @inheritParams power_by_design
#' @param power Target power.
#'
#' @return The required total sample size (a single numeric value).
#' @seealso [power_by_design()], [get_achieved_power()]
#' @export
n_required_by_design <- function(power = 0.80,
                                 p_c, p_s,
                                 design = c("traditional", "IE"),
                                 alpha = 0.025,
                                 R = 1,
                                 pi_pos = NULL) {

  design <- match.arg(design)

  # pooled probability
  p <- (p_c + R * p_s) / (1 + R)

  z_alpha <- qnorm(1 - alpha)
  z_beta  <- qnorm(power)

  # variance scaling
  var_scale <- switch(design,
                      "traditional" = 1,
                      "IE" = {
                        if (is.null(pi_pos)) stop("pi_pos required for IE design")
                        1 / pi_pos
                      })

  var_scale *
    ((1 + R)^2 * (z_alpha + z_beta)^2 * p * (1 - p)) /
    (R * (p_c - p_s)^2)
}

#' Achieved power for a given sample size
#'
#' Inverts a one-argument required-sample-size function to find the power at
#' which the required sample size equals a current sample size. Useful for
#' reporting the power a fixed trial achieves.
#'
#' @param n_required_fun A function of a single argument (power) returning the
#'   required sample size, e.g. [n_required_by_design()] with its other
#'   arguments held fixed.
#' @param N_current The available (fixed) sample size.
#' @param interval Search interval for the power root.
#'
#' @return The achieved power rounded to 3 digits, or the strings `">0.999"` /
#'   `"<0.001"` when the power lies outside the search interval.
#' @export
get_achieved_power <- function(n_required_fun, N_current, interval = c(0.001, 0.999)) {

  f <- function(p) n_required_fun(p) - N_current

  f_low  <- f(interval[1])
  f_high <- f(interval[2])

  if (f_low < 0 && f_high < 0) {
    return(">0.999")
  } else if (f_low > 0 && f_high > 0) {
    return("<0.001")
  } else {
    return(round(uniroot(f, interval = interval)$root, 3))
  }
}
