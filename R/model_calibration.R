# Calibration of the continuous-time multistate natural-history model
# (Lange et al. 2024) to population incidence data by maximum likelihood.
#
# Lange JM, Gogebakan KC, Gulati R, Etzioni R. Projecting the impact of
# multi-cancer early detection on late-stage incidence using multi-state
# disease modeling. Cancer Epidemiol Biomarkers Prev. 2024;33(6):830-7.
#
# All helpers below are internal; the exported entry point is
# fit_natural_history().

# Construct the CTMC rate matrix consistent with a target overall mean sojourn
# time (mst) and late-stage mean sojourn time (dmst). `params` are the
# log-transition rates for the early preclinical states.
get_rate_mat <- function(params, k, mst, dmst) {

  rate_mat <- matrix(0, nrow = k, ncol = k)

  # early preclinical transitions
  for (i in 1:(k - 3)) {
    rate_mat[i, i + 1] <- exp(params[i])
  }

  # late-stage preclinical rate structure (Lange constraints)
  c <- 1 / dmst
  b <- rate_mat[k - 3, k - 2]
  rate_mat[k - 3, k - 1] <- (c + b) / (c * mst) - b
  a <- rate_mat[k - 3, k - 1]
  rate_mat[k - 2, k] <- c

  # finalize diagonal
  diag(rate_mat) <- -apply(rate_mat, 1, sum)
  rate_mat
}

# Cause-specific hazards of early- (h_L) and late-stage (h_D) clinical diagnosis
# at a given age, from the CTMC forward solution MatrixExp(age).
cause_spec_hazards <- function(age, rate_mat, k) {

  prob_mat <- MatrixExp(t = age, mat = rate_mat)

  prob_l <- prob_mat[1, k - 1]    # early-stage clinical
  prob_d <- prob_mat[1, k]        # late-stage clinical
  denom  <- 1 - (prob_l + prob_d)

  num_L <- prob_mat[1, k - 3] * rate_mat[k - 3, k - 1]
  h_L   <- num_L / denom

  num_D <- prob_mat[1, k - 2] * rate_mat[k - 2, k]
  h_D   <- num_D / denom

  list(h_L = h_L, h_D = h_D,
       dens_L = num_L, dens_D = num_D,
       prob_L = prob_l, prob_D = prob_d)
}

# Poisson log-likelihood contribution for one age bin:
# observed_count ~ Poisson(PY * hazard_stage).
indiv_likelihood <- function(age, observed, PY, stage, rate_mat, k) {

  haz_out <- cause_spec_hazards(age = age, rate_mat = rate_mat, k = k)
  thehaz <- if (stage %in% c("local", "Early")) haz_out$h_L else haz_out$h_D

  mean_val <- PY * thehaz
  observed * log(mean_val) - mean_val
}

# Total log-likelihood over all age/stage bins, optimized over `params`.
likelihood_fun <- function(params, mst, dmst, k,
                           allage, allobserved, allPY, allstage) {

  rate_mat <- get_rate_mat(params, k, mst, dmst)
  loglike  <- mapply(allage, allobserved, allPY, allstage,
                     FUN = "indiv_likelihood",
                     MoreArgs = list(rate_mat = rate_mat, k = k))
  sum(unlist(loglike))
}

# Maximize the log-likelihood over `num_seeds` random starts, returning the
# best optim() fit. Uses L-BFGS with a box constraint on the last rate.
get_max_LL <- function(the_data, mst, num_seeds, k, dmst, mean1) {

  c <- 1 / dmst
  M <- mst
  if (c < 1 / M) {
    stop("dmst must be less than mst")
  }

  out_list <- list()
  for (i in 1:num_seeds) {
    out_list[[i]] <- optim(
      par    = -1 * abs(rnorm(n = (k - 3), mean = mean1)),
      fn     = likelihood_fun,
      method = "L-BFGS",
      lower  = rep(-Inf, (k - 3)),
      upper  = c(rep(Inf, (k - 4)), log(c / (c * M - 1))),
      mst = mst, dmst = dmst, k = k,
      allage = the_data$midage,
      allobserved = the_data$Count,
      allPY = the_data$PY,
      allstage = the_data$stage,
      control = list(fnscale = -1, trace = 2, maxit = 20000)
    )
  }

  maxLL <- which.max(unlist(lapply(out_list, "[[", "value")))
  out_list[[maxLL]]
}

# Final sojourn-time summaries for a fitted model.
summarize_fit <- function(the_data, fit, mst, dmst, k) {

  b <- exp(fit$par[k - 3])
  c <- 1 / dmst
  a <- (c + b) / (c * mst) - b

  sojourn_time  <- b / (c * (a + b)) + 1 / (a + b)
  EMST          <- 1 / (a + b)

  list(a = a, b = b,
       sojourn_time = sojourn_time,
       OMST = mst,
       LMST = dmst,
       EMST = EMST)
}

# Optional diagnostic: predicted vs observed incidence by age and stage.
# Returns NULL unless both ggplot2 and reshape2 are available.
plot_observed_expected <- function(the_data, fit, mst, dmst, k) {

  b <- exp(fit$par[k - 3])
  c <- 1 / dmst
  a <- (c + b) / (c * mst) - b

  sojourn_time  <- b / (c * (a + b)) + 1 / (a + b)
  earlypreclin  <- 1 / (a + b)

  rate_mat <- get_rate_mat(fit$par, k = k, mst = mst, dmst = dmst)

  pred_haz <- mapply(
    the_data$midage[1:17],
    FUN = "cause_spec_hazards",
    MoreArgs = list(rate_mat = rate_mat, k = k)
  )

  the_data$pred <- c(unlist(pred_haz[2, ]) * 1e5,
                     unlist(pred_haz[1, ]) * 1e5)

  if (!requireNamespace("reshape2", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    return(list(p1 = NULL, the_data = the_data))
  }

  the_data_m <- reshape2::melt(the_data, measure.vars = c("pred", "rate"))

  p1 <- ggplot2::ggplot(the_data_m, ggplot2::aes(x = .data$midage, y = .data$value)) +
    ggplot2::geom_point(ggplot2::aes(color = .data$variable)) +
    ggplot2::geom_line(ggplot2::aes(color = .data$variable, group = .data$variable)) +
    ggplot2::facet_grid(~stage) +
    ggplot2::xlab("Age") + ggplot2::ylab("Rate per 100,000") +
    ggplot2::ggtitle(paste0(
      "Mean sojourn time = ", round(sojourn_time, 2),
      ", Early preclinical = ", round(earlypreclin, 2)
    ))

  list(p1 = p1, the_data = the_data_m)
}

# Load and prepare provided incidence data for a cancer site (internal). Applies
# the risk multiplier and the standard person-years / rate / age-filter prep
# that fit_natural_history() expects. Made explicit in the calibration vignette.
load_incidence <- function(cancer = "lung", risk_multiplier = 1) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is needed to load provided incidence data; ",
         "install it, or pass your own data via `the_data`.", call. = FALSE)
  }
  xlsx <- system.file("extdata", paste0(cancer, "_incidence.xlsx"),
                      package = "mbiedesign")
  if (!nzchar(xlsx)) {
    stop("No provided incidence data for cancer = '", cancer,
         "'. Pass your own data via `the_data`.", call. = FALSE)
  }
  inc <- openxlsx::read.xlsx(xlsx)
  inc$PY    <- inc$Pop * inc$years
  inc$rate  <- risk_multiplier * inc$Count / inc$PY * 1e5
  inc$Count <- inc$Count * risk_multiplier
  inc[inc$midage > 0 & inc$midage < 80, , drop = FALSE]
}

#' Calibrate the natural-history model to incidence data
#'
#' Fits the continuous-time multistate natural-history model of Lange et al.
#' (2024) to age- and stage-specific incidence data by maximum likelihood,
#' targeting a specified overall mean sojourn time (`OMST`) and late-stage mean
#' sojourn time (`LMST`). By default it loads and calibrates to the provided
#' incidence data for `cancer`; supply `the_data` to fit your own. The returned
#' object supplies the `rate.matrix` that the projection functions (e.g.
#' [get_trial_results_by_design()]) consume. See the calibration vignette
#' (`vignette("calibration", package = "mbiedesign")`) for the data-prep
#' steps done internally.
#'
#' @param OMST Target overall mean sojourn time, in years.
#' @param LMST Target late-stage mean sojourn time, in years.
#' @param cancer Cancer site whose incidence data to fit (currently
#'   `"lung"`). Ignored when `the_data` is supplied.
#' @param k Number of states in the CTMC.
#' @param num_seeds Number of random optimizer restarts; the best fit is kept.
#' @param mean1 Mean of the normal distribution used to draw starting values.
#' @param risk_multiplier Incidence inflation factor applied to the provided
#'   data, e.g. `3.1` to reflect an NLST-eligible high-risk population. Ignored
#'   when `the_data` is supplied.
#' @param the_data Optional incidence data frame (columns `midage`, `Count`,
#'   `Pop`, `years`, `stage`) to fit instead of the provided `cancer` data. Used
#'   as-is, with no `risk_multiplier` applied.
#'
#' @return A list with the raw `optim` fit (`the_fit`), an optional diagnostic
#'   plot (`outplot`, `NULL` if ggplot2/reshape2 are unavailable), the fitted
#'   `rate.matrix`, and a `summary_out` list of sojourn-time quantities
#'   (`OMST`, `LMST`, and the derived `EMST`, the early-stage mean sojourn time).
#'
#' @examples
#' \dontrun{
#' # Fit the lung base case: loads provided incidence, applies the 3.1 multiplier
#' fit <- fit_natural_history(OMST = 4.0, LMST = 1.35)
#' fit$summary_out$EMST  # early-stage mean sojourn time (derived)
#' }
#' @export
fit_natural_history <- function(OMST, LMST, cancer = "lung",
                                k = 16, num_seeds = 1, mean1 = 3,
                                risk_multiplier = 3.1, the_data = NULL) {

  if (is.null(the_data)) {
    the_data <- load_incidence(cancer = cancer, risk_multiplier = risk_multiplier)
  }
  mst  <- OMST
  dmst <- LMST

  the_fit <- tryCatch(
    get_max_LL(the_data = the_data, num_seeds = num_seeds,
               mst = mst, dmst = dmst, k = k, mean1 = mean1),
    error = function(e) e
  )

  while (any(class(the_fit) == "error")) {
    the_fit <- tryCatch(
      get_max_LL(the_data = the_data, num_seeds = num_seeds,
                 mst = mst, dmst = dmst, k = k, mean1 = mean1),
      error = function(e) e
    )
  }

  outplot <- plot_observed_expected(
    the_data = the_data, fit = the_fit, mst = mst, dmst = dmst, k = k
  )

  rate.matrix <- get_rate_mat(params = the_fit$par, k = k, mst = mst, dmst = dmst)
  summary_out <- summarize_fit(the_data, the_fit, mst, dmst, k)

  list(
    the_fit     = the_fit,
    outplot     = outplot$p1,
    fitted_data = outplot$the_data,
    rate.matrix = rate.matrix,
    summary_out = summary_out
  )
}
