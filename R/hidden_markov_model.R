# Hidden Markov model internals for the natural-history / screening process:
# initial-state distribution, emission matrices, interval transition
# probabilities, forward algorithm, and HMM likelihood. Most are unexported;
# get_init(), get_emission_test(), and forward_prob() are exported but marked
# @keywords internal for advanced use (e.g. custom pathway decompositions).

#' Initial state distribution (internal HMM helper)
#'
#' Exposed for advanced use (e.g. custom pathway decompositions); most users
#' will not call this directly.
#'
#' @param rate_matrix Natural-history transition rate matrix.
#' @param start_age Age at trial entry.
#' @return A vector of initial state-occupancy probabilities.
#' @keywords internal
#' @export
get_init <- function(rate_matrix, start_age) {
  # --------------------------------------------------
  # Description:
  # Computes the initial state distribution at the starting age, conditional on
  # not yet being clinically diagnosed. The distribution is restricted to
  # preclinical states, including healthy and preclinical disease states.
  #
  # Inputs:
  # - rate_matrix: transition rate matrix
  # - start_age: starting age
  #
  # Outputs:
  # - Initial state distribution vector over non-clinical states
  # --------------------------------------------------
  k <- nrow(rate_matrix)
  init <- vector()
  init[k] <- 0
  init[k - 1] <- 0
  
  prob_mat_a1 <- MatrixExp(t = start_age, mat = rate_matrix)
  denom <- 1 - prob_mat_a1[1, (k - 1)] - prob_mat_a1[1, k]
  
  for (j in 1:(k - 2)) {
    init[j] <- prob_mat_a1[1, j] / denom
  }
  
  init
}

get_emission_confirmation <- function(k, sens_e, sens_l) {
  # --------------------------------------------------
  # Description:
  # Builds the confirmation-level emission probability matrix linking latent
  # disease states to observed screen-detected and clinically detected states.
  # This matrix is used after the screening episode includes diagnostic
  # confirmation.
  #
  # Appendix mapping:
  # - This is emission matrix E in Equation \ref{eq:emissionE}.
  # - It is used for S_2, S_3, J_4, J_5, and Q_{j,m,k}.
  # - Observed states are:
  #   1 = no observed cancer
  #   2 = screen-detected early-stage cancer
  #   3 = screen-detected late-stage cancer
  #   4 = early-stage clinical diagnosis
  #   5 = late-stage clinical diagnosis
  #
  # Inputs:
  # - k: number of latent states
  # - sens_e: early-stage sensitivity
  # - sens_l: late-stage sensitivity
  #
  # Outputs:
  # - Confirmation-level emission probability matrix E
  # --------------------------------------------------
  emit <- matrix(0, nrow = k, ncol = 5)
  
  emit[1:(k - 4), 1] <- 1
  
  emit[(k - 3), 2] <- sens_e
  emit[(k - 3), 1] <- 1 - sens_e
  
  emit[(k - 2), 3] <- sens_l
  emit[(k - 2), 1] <- 1 - sens_l
  
  emit[(k - 1), 4] <- 1
  emit[k, 5] <- 1
  
  emit
}

#' Primary-test emission matrix (internal HMM helper)
#'
#' Exposed for advanced use; most users will not call this directly.
#'
#' @param k Number of states in the model.
#' @param sens_e Early-stage test sensitivity.
#' @param sens_l Late-stage test sensitivity.
#' @param specificity Test specificity.
#' @return The primary-test emission probability matrix.
#' @keywords internal
#' @export
get_emission_test <- function(k, sens_e, sens_l, specificity) {
  # --------------------------------------------------
  # Description:
  # Builds the primary-test emission probability matrix for IE positivity.
  # This matrix models whether the screening test itself is positive or negative
  # before diagnostic confirmation.
  #
  # Appendix mapping:
  # - This is emission matrix E* in Equation \ref{eq:emissionFP}.
  # - It is used for pi(T), S_test(a_l), EP_{m,k}, and alpha_m.
  # - We use a 5-column version of E* so that clinical states have the same
  #   labels as E:
  #   1 = test negative
  #   2 = test positive
  #   3 = unused/reserved
  #   4 = early-stage clinical diagnosis
  #   5 = late-stage clinical diagnosis
  # - Therefore, in the appendix notation, O*_m = 2 represents first test
  #   positivity. The condition X_m in no-cancer latent states makes it a
  #   false-positive/unconfirmed-positive pathway.
  #
  # Inputs:
  # - k: number of latent states
  # - sens_e: early-stage sensitivity
  # - sens_l: late-stage sensitivity
  # - specificity: test specificity
  #
  # Outputs:
  # - Primary-test emission probability matrix E*
  # --------------------------------------------------
  emit <- matrix(0, nrow = k, ncol = 5)
  
  # No-cancer latent states: positive test is false-positive.
  emit[1:(k - 4), 1] <- specificity
  emit[1:(k - 4), 2] <- 1 - specificity
  
  # Early preclinical disease: positive test is true-positive.
  emit[(k - 3), 1] <- 1 - sens_e
  emit[(k - 3), 2] <- sens_e
  
  # Late preclinical disease: positive test is true-positive.
  emit[(k - 2), 1] <- 1 - sens_l
  emit[(k - 2), 2] <- sens_l
  
  # Clinical disease states.
  emit[(k - 1), 4] <- 1
  emit[k, 5] <- 1
  
  emit
}

get_times <- function(screen_times, preclin = FALSE, post = FALSE, extended_followup = NULL) {
  # --------------------------------------------------
  # Description:
  # Creates times for screen detection, interval detection, or postscreen
  # detection under the standard HMM likelihood calculation.
  #
  # Inputs:
  # - screen_times: vector of screen times
  # - preclin: whether the event is screen-detected
  # - post: whether the event is postscreen clinical detection
  # - extended_followup: additional follow-up time after the last screen
  #
  # Outputs:
  # - List of observation-time vectors
  # --------------------------------------------------
  if (post) {
    if (is.null(extended_followup)) stop("Need extended_followup when post = TRUE")
    return(list(c(screen_times, max(screen_times) + extended_followup)))
  }
  
  n <- length(screen_times)
  the_times <- list()
  
  if (n == 1) {
    if (preclin) {
      the_times[[1]] <- c(0)
    } else {
      if (is.null(extended_followup)) {
        stop("Need extended_followup after single exam when preclin = FALSE")
      } else {
        the_times[[1]] <- c(0, extended_followup)
      }
    }
    return(the_times)
  }
  
  if (preclin) {
    for (i in seq_len(n)) {
      the_times[[i]] <- screen_times[1:i]
    }
  } else {
    for (i in seq_len(n - 1)) {
      the_times[[i]] <- screen_times[1:(i + 1)]
    }
    if (!is.null(extended_followup)) {
      the_times[[length(the_times) + 1]] <- c(screen_times, tail(screen_times, 1) + extended_followup)
    }
  }
  
  the_times
}

get_data <- function(screen_times,
                     stage = c("early", "late"),
                     preclin = FALSE,
                     post = FALSE,
                     extended_followup = NULL) {
  # --------------------------------------------------
  # Description:
  # Creates observed states for screen-detected, interval-detected, or
  # postscreen-detected outcomes under the standard HMM likelihood calculation.
  #
  # Inputs:
  # - screen_times: vector of screen times
  # - stage: "early" or "late"
  # - preclin: whether the event is screen-detected
  # - post: whether the event is postscreen clinical detection
  # - extended_followup: additional follow-up time after the last screen
  #
  # Outputs:
  # - List of observed-state vectors
  # --------------------------------------------------
  stage <- match.arg(stage)
  
  preclin_code <- if (stage == "early") 2 else 3
  clinical_code <- if (stage == "early") 4 else 5
  
  if (post) {
    if (is.null(extended_followup)) stop("Need extended_followup when post = TRUE")
    return(list(c(rep(1, length(screen_times)), clinical_code)))
  }
  
  n <- length(screen_times)
  the_data <- list()
  
  if (n == 1) {
    if (preclin) {
      the_data[[1]] <- c(preclin_code)
    } else {
      if (is.null(extended_followup)) {
        stop("Need extended_followup after single exam when preclin = FALSE")
      } else {
        the_data[[1]] <- c(1, clinical_code)
      }
    }
    return(the_data)
  }
  
  if (preclin) {
    the_data[[1]] <- c(preclin_code)
    for (i in seq_len(n - 1)) {
      the_data[[i + 1]] <- c(rep(1, i), preclin_code)
    }
  } else {
    for (i in seq_len(n - 1)) {
      the_data[[i]] <- c(rep(1, i), clinical_code)
    }
    if (!is.null(extended_followup)) {
      the_data[[length(the_data) + 1]] <- c(rep(1, n), clinical_code)
    }
  }
  
  the_data
}

transition_prob_all <- function(time_intervals_list, rate_matrix_list, exact_time_rank = NULL) {
  # --------------------------------------------------
  # author JL 2/7/2011
  # Description:
  # Computes lists of transition probability matrices for several observation
  # intervals under a time-homogeneous continuous-time Markov chain:
  # P(t2 - t1), P(t3 - t2), ..., P(tn - t_{n-1}).
  #
  # If an interval corresponds to an exact event time, that interval is evaluated
  # using a transition density rather than a transition probability. This is
  # implemented by multiplying the transition probability matrix by the
  # off-diagonal part of the rate matrix.
  #
  # Inputs:
  # - time_intervals_list: list of interval vectors
  # - rate_matrix_list: list of transition rate matrices
  # - exact_time_rank: optional index of the interval to evaluate as a density
  #   (i.e., exact event time)
  #
  # Outputs:
  # - List of transition probability-matrix lists, one for each observation pattern
  # --------------------------------------------------
  mapply(
    FUN = function(time.intervals, rate.matrix) {
      probs.list <- lapply(time.intervals, FUN = "MatrixExp", mat = rate.matrix)
      
      if (!is.null(exact_time_rank)) {
        if (exact_time_rank < 1 || exact_time_rank > length(probs.list)) {
          stop("exact_time_rank must be between 1 and length(time.intervals)")
        }
        
        rate_nodiag <- rate.matrix
        diag(rate_nodiag) <- 0
        probs.list[[exact_time_rank]] <- probs.list[[exact_time_rank]] %*% rate_nodiag
      }
      
      probs.list
    },
    time.intervals = time_intervals_list,
    rate.matrix = rate_matrix_list,
    SIMPLIFY = FALSE
  )
}

#' Forward algorithm (internal HMM helper)
#'
#' Exposed for advanced use (e.g. custom pathway decompositions); most users
#' will not call this directly.
#'
#' @param obs Observation sequence.
#' @param times Observation times.
#' @param delta Initial state distribution.
#' @param rate_matrix Natural-history transition rate matrix.
#' @param emission_matrix Emission probability matrix.
#' @param return_state Logical; if `TRUE`, return the state-probability vector.
#' @param start_time Optional start time for the recursion.
#' @return The forward probability (or state vector if `return_state = TRUE`).
#' @keywords internal
#' @export
forward_prob <- function(obs, times, delta, rate_matrix, emission_matrix,
                         return_state = FALSE, start_time = NULL) {
  # --------------------------------------------------
  # Description:
  # Computes a forward probability for an observed sequence without the backward
  # recursion. This is needed for the IE screen-arm clinical pathway because the
  # calculation uses two emission matrices: first E* to identify first test
  # positivity, then E to propagate confirmed screening/clinical outcomes.
  #
  # New function.
  #
  # Appendix mapping:
  # - With E* and return_state = TRUE, this returns the raw forward vector after
  #   O*_0=1,...,O*_{m-1}=1,O*_m=2. After retaining no-cancer latent states, this
  #   corresponds to alpha_m(1_j) B_FP.
  # - With E, the second call computes the Q_{j,m,k} pathway: no observed cancer
  #   at later screens followed by clinical diagnosis at the endpoint.
  # - start_time is required for the second call because the starting vector is
  #   already located at a_m.
  #
  # Inputs:
  # - obs: observed state sequence
  # - times: observation times corresponding to obs
  # - delta: starting state vector; can be normalized or unnormalized
  # - rate_matrix: transition rate matrix
  # - emission_matrix: emission matrix for the observed process
  # - return_state: if TRUE, return the raw state vector; otherwise return its sum
  # - start_time: optional time corresponding to the starting vector delta
  #
  # Outputs:
  # - Total probability of the observed sequence or the raw state vector
  # --------------------------------------------------
  if (length(obs) != length(times)) {
    stop("obs and times must have the same length.")
  }
  
  phi <- matrix(delta, nrow = 1)
  
  for (i in seq_along(obs)) {
    if (i == 1 && !is.null(start_time)) {
      if (times[i] < start_time) stop("times must be >= start_time.")
      P <- MatrixExp(t = times[i] - start_time, mat = rate_matrix)
      phi <- phi %*% P
    }
    
    if (i > 1) {
      P <- MatrixExp(t = times[i] - times[i - 1], mat = rate_matrix)
      phi <- phi %*% P
    }
    
    phi <- phi * matrix(emission_matrix[, obs[i]], nrow = 1)
  }
  
  if (return_state) {
    return(as.numeric(phi))
  }
  
  as.numeric(sum(phi))
}


likelihood <- function(rates_list, init_list, emission_list,
                       obs_data_list, obs_times_list,
                       exact_time_rank = NULL) {
  # --------------------------------------------------
  # Description:
  # Computes the likelihood across one or more observation patterns under the
  # hidden Markov model.
  #
  # Inputs:
  # - rates_list: list of transition rate matrices
  # - init_list: list of initial state distributions
  # - emission_list: list of emission matrices
  # - obs_data_list: list of observed state sequences
  # - obs_times_list: list of observation times
  # - exact_time_rank: optional index of interval with exact-event adjustment
  #
  # Outputs:
  # - Total likelihood across the given observations
  # --------------------------------------------------
  time_diffs_list <- lapply(obs_times_list, diff)
  
  transition_probabilities_list <- transition_prob_all(
    time_intervals_list = time_diffs_list,
    rate_matrix_list = rates_list,
    exact_time_rank = exact_time_rank
  )
  
  likelihood_forward_backward_list <- mapply(
    obs_data_list,
    transition_probabilities_list,
    init_list,
    emission_list,
    FUN = forwardback,
    SIMPLIFY = FALSE
  )
  
  LL <- vapply(likelihood_forward_backward_list, function(z) as.numeric(z$LL), numeric(1))
  sum(exp(LL))
}


forwardback <- function(x, Pi, delta, emission_matrix) {
  # --------------------------------------------------
  # Author: JL, 2/4/2011
  # Description:
  # Computes forward and backward probabilities and the observed-data log
  # likelihood for a hidden Markov model using a single emission matrix.
  # This function is retained for the standard HMM likelihood calculations.
  #
  # Inputs:
  # - x: observed data sequence
  # - Pi: list of transition probability matrices from t1-t2, t2-t3, ..., t_{n-1}-t_n
  # - delta: initial distribution of hidden states
  # - emission_matrix: matrix of emission probabilities; row i corresponds to
  #   latent state i and column k to observed state k
  #
  # Outputs:
  # - List containing logalpha, logbeta, and LL, the observed-data log likelihood
  # --------------------------------------------------
  if (length(x) == 1) {
    LL <- as.numeric(log(sum(emission_matrix[, x] * delta)))
    return(list(LL = LL))
  }
  
  m <- nrow(Pi[[1]])
  n <- length(x)
  
  phi <- matrix(delta, nrow = 1)
  logalpha <- matrix(rep(NA_real_, m * n), nrow = n)
  lscale <- 0
  
  for (i in seq_len(n)) {
    if (i > 1) {
      phi <- phi %*% Pi[[i - 1]]
    }
    phi <- phi %*% diag(emission_matrix[, x[i]])
    sumphi <- as.numeric(sum(phi))
    phi <- phi / sumphi
    lscale <- lscale + log(sumphi)
    logalpha[i, ] <- as.numeric(log(phi) + lscale)
  }
  
  LL1 <- as.numeric(lscale)
  
  logbeta <- matrix(rep(NA_real_, m * n), nrow = n)
  logbeta[n, ] <- 0
  phi <- matrix(rep(1 / m, m), ncol = 1)
  lscale <- log(m)
  
  for (i in seq(n - 1, 1, -1)) {
    phi <- Pi[[i]] %*% diag(emission_matrix[, x[i + 1]]) %*% phi
    logbeta[i, ] <- as.numeric(log(phi) + lscale)
    sumphi <- as.numeric(sum(phi))
    phi <- phi / sumphi
    lscale <- lscale + log(sumphi)
  }
  
  list(
    logalpha = logalpha,
    logbeta = logbeta,
    LL = LL1
  )
}
