# Trial-outcome projection by design: traditional (standard) and
# intended-effect (IE). 

#####################################################################
# 2. STANDARD TRIAL DESIGN INCIDENCE PROJECTION FUNCTIONS
#####################################################################

get_control_incidence_standard <- function(rate_matrix, age1, age2, stage = c("early", "late")) {
  # --------------------------------------------------
  # Description:
  # Computes control-arm incidence in an age interval under the standard design.
  #
  # Inputs:
  # - rate_matrix: transition rate matrix
  # - age1: start of interval
  # - age2: end of interval
  # - stage: "early" or "late"
  #
  # Outputs:
  # - Interval incidence in the control arm
  # --------------------------------------------------
  stage <- match.arg(stage)
  
  prob_mat_a1 <- MatrixExp(t = age1, mat = rate_matrix)
  prob_mat_a2 <- MatrixExp(t = age2, mat = rate_matrix)
  
  k <- nrow(rate_matrix)
  clinical_idx <- if (stage == "early") (k - 1) else k
  
  (prob_mat_a2[1, clinical_idx] - prob_mat_a1[1, clinical_idx]) /
    (1 - prob_mat_a1[1, k] - prob_mat_a1[1, k - 1])
}

get_screen_incidence_standard <- function(
    sens_e, sens_l,
    screen_times, age1, rate_matrix,
    stage = c("early", "late"),
    preclin = FALSE,
    post = FALSE,
    extended_followup = NULL
) {
  # --------------------------------------------------
  # Description:
  # Computes screen-arm incidence under the standard design (screen-detected cases,
  # interval-detected cases, or cases detected in follow-up period after last screen).
  #
  # Inputs:
  # - sens_e: early-stage sensitivity
  # - sens_l: late-stage sensitivity
  # - screen_times: vector of screen times
  # - age1: starting age
  # - rate_matrix: transition rate matrix
  # - stage: "early" or "late"
  # - preclin: whether the event is screen-detected
  # - post: whether the event is postscreen clinical detection
  # - extended_followup: additional follow-up after the last screen
  #
  # Outputs:
  # - Cumulative incidence for the specified screen-arm detection type
  # --------------------------------------------------
  stage <- match.arg(stage)
  
  k <- nrow(rate_matrix)
  init_dist <- get_init(rate_matrix, age1)
  emission_dist <- get_emission_confirmation(k, sens_e = sens_e, sens_l = sens_l) # Appendix E
  
  obs_data_list <- get_data(
    screen_times = screen_times,
    stage = stage,
    preclin = preclin,
    post = post,
    extended_followup = extended_followup
  )
  
  obs_times_list <- get_times(
    screen_times = screen_times,
    preclin = preclin,
    post = post,
    extended_followup = extended_followup
  )
  
  init_list     <- rep(list(init_dist), length(obs_data_list))
  emission_list <- rep(list(emission_dist), length(obs_data_list))
  rates_list    <- rep(list(rate_matrix), length(obs_data_list))
  
  likelihood(
    rates_list = rates_list,
    init_list = init_list,
    emission_list = emission_list,
    obs_data_list = obs_data_list,
    obs_times_list = obs_times_list
  )
}

#####################################################################
# 3. INTENDED EFFECT DESIGN INCIDENCE PROJECTION FUNCTIONS
#####################################################################

#' Ever-positive fraction under the intended-effect design
#'
#' Computes the intended-effect (IE) ever-positive fraction: the probability of
#' at least one positive primary screening test across the scheduled screening
#' examinations, including both true-positive and false-positive primary-test
#' results.
#'
#' @param sens_e Early-stage test sensitivity.
#' @param sens_l Late-stage test sensitivity.
#' @param specificity Test specificity.
#' @param screen_times Numeric vector of scheduled screen times.
#' @param start_age Age at trial entry.
#' @param rate_matrix Natural-history transition rate matrix (e.g. from a fitted
#'   model via [load_fitted_model()]).
#'
#' @return A single numeric value: the ever-positive probability.
#' @export
get_positivity_rate <- function(sens_e, sens_l, specificity, screen_times, start_age, rate_matrix) {
  # --------------------------------------------------
  # Description:
  # Computes the IE ever-positive fraction: the probability of at least one
  # positive primary screening test across scheduled screening examinations.
  # This includes both true-positive and false-positive primary test results.
  #
  # Appendix mapping:
  # - Computes pi(T) = sum_l S_test(a_l).
  # - Uses E* from get_emission_test().
  # - S_test(a_m) = Pr(O*_0=1,...,O*_{m-1}=1,O*_m=2).
  #
  # Inputs:
  # - sens_e: early-stage sensitivity
  # - sens_l: late-stage sensitivity
  # - specificity: test specificity
  # - screen_times: vector of screen times
  # - start_age: starting age
  # - rate_matrix: transition rate matrix
  #
  # Outputs:
  # - Ever-positive probability pi(T)
  # --------------------------------------------------
  k <- nrow(rate_matrix)
  init_dist <- get_init(rate_matrix, start_age)
  emission_test <- get_emission_test(k, sens_e = sens_e, sens_l = sens_l, specificity = specificity) # Appendix E*
  
  sum(sapply(seq_along(screen_times), function(m) {
    forward_prob(
      obs = c(rep(1, m - 1), 2),
      times = screen_times[seq_len(m)],
      delta = init_dist,
      rate_matrix = rate_matrix,
      emission_matrix = emission_test # Appendix S_test(a_m) using E*
    )
  }))
}

get_control_incidence_ie <- function(
    sens_e, sens_l, specificity,
    screen_times, start_age, age1, age2,
    rate_matrix,
    stage = c("early", "late")
) {
  # --------------------------------------------------
  # Description:
  # Computes control-arm interval incidence under the intended-effect design,
  # conditional on ever-positive status.
  #
  # Appendix mapping:
  # - Computes L_4(a_l,a_{l+1}] or L_5(a_l,a_{l+1}] divided by pi(T).
  # - Uses E* because the control-arm IE population is defined by hypothetical
  #   primary test positivity.
  # - EP_{m,k} = {O*_0=1,...,O*_{m-1}=1,O*_m=2,Y_stage <= a_k}.
  # - The line k_time <= screen_times[m] returns 0, so the interval immediately
  #   after first positivity, (a_m,a_{m+1}], is included correctly.
  #
  # Inputs:
  # - sens_e: early-stage sensitivity
  # - sens_l: late-stage sensitivity
  # - specificity: test specificity
  # - screen_times: full vector of screen times
  # - start_age: starting age
  # - age1: start of interval
  # - age2: end of interval
  # - rate_matrix: transition rate matrix
  # - stage: "early" or "late"
  #
  # Outputs:
  # - Interval incidence in the control arm under the intended-effect design
  # --------------------------------------------------
  stage <- match.arg(stage)
  clinical_code <- if (stage == "early") 4 else 5
  
  positivity_rate <- get_positivity_rate( # Appendix pi(T)
    sens_e = sens_e,
    sens_l = sens_l,
    specificity = specificity,
    screen_times = screen_times,
    start_age = start_age,
    rate_matrix = rate_matrix
  )
  
  if (positivity_rate <= 0) return(0)
  
  k <- nrow(rate_matrix)
  init_dist <- get_init(rate_matrix, start_age)
  emission_test <- get_emission_test(k, sens_e = sens_e, sens_l = sens_l, specificity = specificity) # Appendix E*
  
  t1 <- age1 - start_age
  t2 <- age2 - start_age
  
  if (t2 <= t1) {
    stop("age2 must be greater than age1.")
  }
  
  prob_first_positive_at_m_and_clinical_by_k <- function(m, k_time) {
    if (m < 1 || m > length(screen_times)) {
      stop("m is out of range.")
    }
    
    # Appendix EP_{m,m}: simultaneous first positivity and clinical diagnosis
    # at a_m has probability zero.
    if (k_time <= screen_times[m]) {
      return(0)
    }
    
    forward_prob(
      obs = c(rep(1, m - 1), 2, clinical_code),
      times = c(screen_times[seq_len(m)], k_time),
      delta = init_dist,
      rate_matrix = rate_matrix,
      emission_matrix = emission_test # Appendix EP_{m,k} using E*
    )
  }
  
  valid_screens <- which(screen_times <= t1)
  
  if (length(valid_screens) == 0) {
    return(0)
  }
  
  interval_probability <- sum(
    sapply(valid_screens, function(m) {
      prob_first_positive_at_m_and_clinical_by_k(m, t2) -
        prob_first_positive_at_m_and_clinical_by_k(m, t1)
    })
  ) / positivity_rate
  
  interval_probability
}

get_screen_incidence_ie <- function(
    sens_e, sens_l, specificity,
    screen_times, start_age,
    rate_matrix,
    stage = c("early", "late"),
    preclin = FALSE,
    post = FALSE,
    extended_followup = NULL,
    positivity_rate = NULL
) {
  # --------------------------------------------------
  # Description:
  # Computes screen-arm incidence under the intended-effect design for
  # screen-detected, interval-detected, or postscreen outcomes.
  #
  # Appendix mapping:
  # - If preclin = TRUE: computes S_2 or S_3 divided by pi(T), using E.
  # - If preclin = FALSE: computes clinical diagnosis after entering the IE
  #   population through a positive primary test without confirmed cancer.
  # - alpha line: E* gives the raw state vector after first test positivity.
  #   Retaining no-cancer states gives alpha_m(1_j) B_FP.
  # - Q line: E is then used for later screening episodes and clinical diagnosis.
  # - p(k) = min(k-1,s-1) is implemented by later_screens, which contains only
  #   actual screening times after the first positive test and before k_time.
  # - R line: cumulative EFP_{m,k} contributions are computed at endpoints;
  #   interval R values are obtained by differencing these cumulative values over time.
  #
  # Inputs:
  # - sens_e: early-stage sensitivity
  # - sens_l: late-stage sensitivity
  # - specificity: test specificity
  # - screen_times: vector of screen times
  # - start_age: starting age
  # - rate_matrix: transition rate matrix
  # - stage: "early" or "late"
  # - preclin: whether the event is screen-detected
  # - post: whether the event is postscreen clinical detection
  # - extended_followup: additional follow-up after the last screen
  # - positivity_rate: optional fixed ever-positive probability pi(T)
  #
  # Outputs:
  # - Cumulative incidence for the specified IE screen-arm detection type
  # --------------------------------------------------
  stage <- match.arg(stage)
  clinical_code <- if (stage == "early") 4 else 5
  
  if (is.null(positivity_rate)) {
    positivity_rate <- get_positivity_rate( # Appendix pi(T)
      sens_e = sens_e,
      sens_l = sens_l,
      specificity = specificity,
      screen_times = screen_times,
      start_age = start_age,
      rate_matrix = rate_matrix
    )
  }
  
  if (positivity_rate <= 0) return(0)
  
  # Appendix S_2/S_3 line:
  # Screen-detected IE incidence uses E. Prior unconfirmed positive tests are
  # observed as O=1 under E, so later screen-detected cancers are already included.
  if (preclin) {
    return(
      get_screen_incidence_standard(
        sens_e = sens_e,
        sens_l = sens_l,
        screen_times = screen_times,
        age1 = start_age,
        rate_matrix = rate_matrix,
        stage = stage,
        preclin = TRUE
      ) / positivity_rate
    )
  }
  
  if (is.null(extended_followup)) {
    stop("Need extended_followup when preclin = FALSE")
  }
  
  k <- nrow(rate_matrix)
  init_dist <- get_init(rate_matrix, start_age)
  emission_test <- get_emission_test(k, sens_e = sens_e, sens_l = sens_l, specificity = specificity) # Appendix E*
  emission_confirmation <- get_emission_confirmation(k, sens_e = sens_e, sens_l = sens_l) # Appendix E
  
  # Endpoint times mimic the standard interval/postscreen structure.
  if (post) {
    endpoint_times <- max(screen_times) + extended_followup
  } else if (length(screen_times) == 1) {
    endpoint_times <- max(screen_times) + extended_followup
  } else {
    endpoint_times <- screen_times[-1]
    endpoint_times <- c(endpoint_times, max(screen_times) + extended_followup)
  }
  
  # Appendix R line:
  # For each endpoint k_time, the inner sum computes cumulative EFP_{m,k}
  # contributions. Interval R values are obtained by differencing.
  clinical_probability <- sum(sapply(endpoint_times, function(k_time) {
    valid_screens <- which(screen_times < k_time)
    
    if (length(valid_screens) == 0) {
      return(0)
    }
    
    sum(sapply(valid_screens, function(m) {
      # Appendix alpha line:
      # alpha = Pr(O*_0=1,...,O*_{m-1}=1,O*_m=2,X_m=j), computed with E*.
      alpha <- forward_prob(
        obs = c(rep(1, m - 1), 2),
        times = screen_times[seq_len(m)],
        delta = init_dist,
        rate_matrix = rate_matrix,
        emission_matrix = emission_test, # Appendix E*
        return_state = TRUE
      )
      
      # Appendix alpha_m(1_j) B_FP line:
      # Retain only no-cancer latent states. This isolates first positive tests
      # that did not correspond to confirmed cancer at the same screening episode.
      alpha_unconfirmed <- rep(0, k)
      alpha_unconfirmed[1:(k - 4)] <- alpha[1:(k - 4)]
      
      if (sum(alpha_unconfirmed) <= 0) {
        return(0)
      }
      
      # Appendix p(k) line:
      # p(k)=min(k-1,s-1) is implemented by including only actual later
      # screening times after a_m and before k_time.
      later_screens <- screen_times[
        screen_times > screen_times[m] &
          screen_times < k_time
      ]
      
      # Appendix Q line:
      # Q_{j,m,k} is computed with E: no observed cancer at later screening
      # episodes, followed by clinical diagnosis at k_time.
      forward_prob(
        obs = c(rep(1, length(later_screens)), clinical_code),
        times = c(later_screens, k_time),
        delta = alpha_unconfirmed,
        rate_matrix = rate_matrix,
        emission_matrix = emission_confirmation, # Appendix E
        start_time = screen_times[m]
      )
    }))
  }))
  
  clinical_probability / positivity_rate
}

#####################################################################
# 4. STAGE-SPECIFIC TIME-DEPENDENT INCIDENCE PROJECTION FUNCTIONS
#    BY TRIAL DESIGN AND ARM
#####################################################################

get_screen_arm_trial_results <- function(stage, design, params) {
  # --------------------------------------------------
  # Description:
  # Computes screen-arm incidence and counts for a given stage and trial design.
  # The same structure is used for standard and IE designs: screen-detected,
  # interval-detected, and postscreen-detected components are computed and then
  # combined into total screen-arm detected incidence.
  #
  # Appendix mapping:
  # - For standard design, this assembles p_s(T) from S_2/S_3 and J_4/J_5.
  # - For IE design, this assembles p_{s,IE}(T) among ever-positive participants.
  # - IE screen-detected calls use E through get_screen_incidence_ie(preclin=TRUE).
  # - IE clinical calls use E* followed by E through get_screen_incidence_ie(preclin=FALSE).
  #
  # Inputs:
  # - stage: "early" or "late"
  # - design: "standard" or "IE"
  # - params: list of trial design, test, natural history, and sample size inputs
  #
  # Outputs:
  # - List with yearly and cumulative incidence/counts for screen-detected,
  #   interval-detected, and total screen-arm detected outcomes
  # --------------------------------------------------
  stage <- match.arg(stage, c("early", "late"))
  design <- match.arg(design, c("standard", "IE"))
  
  n_intervals <- length(params$interval_ends)
  
  
  positivity_rate <- if (design == "IE") {
    get_positivity_rate( # Appendix pi(T)
      sens_e = params$sens_e,
      sens_l = params$sens_l,
      specificity = params$specificity,
      screen_times = params$screen_times,
      start_age = params$start_age,
      rate_matrix = params$rate_matrix
    )
  } else {
    1
  }
  
  n_screen_effective <- params$n_screen * positivity_rate
  
  get_inc <- function(times, preclin = TRUE, post = FALSE, followup = NULL) {
    switch(
      design,
      standard = get_screen_incidence_standard(
        sens_e = params$sens_e,
        sens_l = params$sens_l,
        screen_times = times,
        age1 = params$start_age,
        rate_matrix = params$rate_matrix,
        stage = stage,
        preclin = preclin,
        post = post,
        extended_followup = followup
      ),
      IE = get_screen_incidence_ie(
        sens_e = params$sens_e,
        sens_l = params$sens_l,
        specificity = params$specificity,
        screen_times = times,
        start_age = params$start_age,
        rate_matrix = params$rate_matrix,
        stage = stage,
        preclin = preclin,
        post = post,
        extended_followup = followup,
        positivity_rate = positivity_rate
      )
    )
  }
  
  screen_detected_screen <- sapply(seq_len(params$numscreens), function(i) {
    get_inc(params$screen_times[seq_len(i)], preclin = TRUE)
  })
  
  cum_screen_detected <- numeric(n_intervals)
  
  for (i in seq_len(n_intervals)) {
    t_end <- params$interval_ends[i]
    screen_index <- sum(params$screen_times < t_end)
    
    if (screen_index == 0) {
      cum_screen_detected[i] <- 0
    } else {
      screen_index <- min(screen_index, length(screen_detected_screen))
      cum_screen_detected[i] <- screen_detected_screen[screen_index]
    }
  }
  
  interval_screen <- sapply(seq_len(params$numscreens), function(i) {
    get_inc(
      params$screen_times[seq_len(i)],
      preclin = FALSE,
      post = FALSE,
      followup = params$screen_int
    )
  })
  
  cum_interval_detected <- numeric(n_intervals)
  
  for (i in seq_len(n_intervals)) {
    t_end <- params$interval_ends[i]
    screen_index <- sum(params$screen_times < t_end)
    
    if (screen_index == 0) {
      cum_interval_detected[i] <- 0
    } else if (screen_index < params$numscreens) {
      screen_index <- min(screen_index, length(interval_screen))
      cum_interval_detected[i] <- interval_screen[screen_index]
    } else {
      base_before_post <- if (params$numscreens > 1) {
        interval_screen[params$numscreens - 1]
      } else {
        0
      }
      
      cum_interval_detected[i] <- base_before_post +
        get_inc(
          params$screen_times,
          preclin = FALSE,
          post = TRUE,
          followup = t_end - max(params$screen_times)
        )
    }
  }
  
  cum_screen_arm_detected <- cum_screen_detected + cum_interval_detected
  
  yearly_screen_detected <- diff(c(0, cum_screen_detected))
  yearly_interval_detected <- diff(c(0, cum_interval_detected))
  yearly_screen_arm_detected <- yearly_screen_detected + yearly_interval_detected
  
  list(
    incidence = list(
      yearly = list(
        screen_detected = yearly_screen_detected,
        interval_detected = yearly_interval_detected,
        screen_arm_detected = yearly_screen_arm_detected
      ),
      cumulative = list(
        screen_detected = cum_screen_detected,
        interval_detected = cum_interval_detected,
        screen_arm_detected = cum_screen_arm_detected
      )
    ),
    count = list(
      yearly = list(
        screen_detected = yearly_screen_detected * n_screen_effective,
        interval_detected = yearly_interval_detected * n_screen_effective,
        screen_arm_detected = yearly_screen_arm_detected * n_screen_effective
      ),
      cumulative = list(
        screen_detected = cum_screen_detected * n_screen_effective,
        interval_detected = cum_interval_detected * n_screen_effective,
        screen_arm_detected = cum_screen_arm_detected * n_screen_effective
      )
    )
  )
}

get_control_arm_trial_results <- function(stage, design, params) {
  # --------------------------------------------------
  # Description:
  # Computes control-arm incidence and counts for a given stage and trial design.
  #
  # Appendix mapping:
  # - For standard design, this computes p_c(T) from clinical incidence intervals.
  # - For IE design, this computes p_{c,IE}(T), conditioning on hypothetical
  #   ever-positive status using E*.
  #
  # Inputs:
  # - stage: "early" or "late"
  # - design: "standard" or "IE"
  # - params: list of trial design, test, natural history, and sample size inputs
  #
  # Outputs:
  # - List with yearly and cumulative incidence/counts for control-arm detected outcomes
  # --------------------------------------------------
  stage <- match.arg(stage, c("early", "late"))
  design <- match.arg(design, c("standard", "IE"))
  
  n_intervals <- length(params$interval_ends)
  
  
  positivity_rate <- if (design == "IE") {
    get_positivity_rate( # Appendix pi(T)
      sens_e = params$sens_e,
      sens_l = params$sens_l,
      specificity = params$specificity,
      screen_times = params$screen_times,
      start_age = params$start_age,
      rate_matrix = params$rate_matrix
    )
  } else {
    1
  }
  
  n_control_effective <- params$n_control * positivity_rate
  
  yearly_control_arm_detected <- mapply(
    function(t1, t2) {
      switch(
        design,
        standard = get_control_incidence_standard(
          rate_matrix = params$rate_matrix,
          age1 = params$start_age + t1,
          age2 = params$start_age + t2,
          stage = stage
        ),
        IE = get_control_incidence_ie(
          sens_e = params$sens_e,
          sens_l = params$sens_l,
          specificity = params$specificity,
          screen_times = params$screen_times,
          start_age = params$start_age,
          age1 = params$start_age + t1,
          age2 = params$start_age + t2,
          rate_matrix = params$rate_matrix,
          stage = stage
        )
      )
    },
    t1 = params$interval_starts,
    t2 = params$interval_ends
  )
  
  cum_control_arm_detected <- cumsum(yearly_control_arm_detected)
  
  list(
    incidence = list(
      yearly = list(control_arm_detected = yearly_control_arm_detected),
      cumulative = list(control_arm_detected = cum_control_arm_detected)
    ),
    count = list(
      yearly = list(control_arm_detected = yearly_control_arm_detected * n_control_effective),
      cumulative = list(control_arm_detected = cum_control_arm_detected * n_control_effective)
    )
  )
}

#####################################################################
# 5. MAIN FUNCTIONS
#####################################################################

#' Project trial outcomes for every design
#'
#' Main projection function. Given a calibrated natural-history model and a
#' screening schedule, projects control- and screen-arm cancer outcomes under
#' the traditional (`"standard"`) and intended-effect (`"IE"`) designs, for
#' early-stage, late-stage, and overall cancer.
#'
#' @param start_age Age at trial entry.
#' @param numscreens Number of screening examinations.
#' @param screen_int Interval between screens (years).
#' @param num_followup_intervals Length of follow-up after the final screen.
#' @param rate_matrix Natural-history transition rate matrix (e.g. from a fitted
#'   model via [load_fitted_model()]).
#' @param sens_e Early-stage test sensitivity.
#' @param sens_l Late-stage test sensitivity.
#' @param specificity Test specificity.
#' @param n_control Control-arm sample size.
#' @param n_screen Screen-arm sample size.
#'
#' @return A list with the projection `params`, the trial timeline, and
#'   control-arm and screen-arm outcomes plus stage-shift summaries for the
#'   standard and IE designs.
#' @seealso [stage_shift_by_design()] to extract a formatted table for one design.
#' @export
get_trial_results_by_design <- function(
    start_age,
    numscreens,
    screen_int,
    num_followup_intervals,
    rate_matrix,
    sens_e,
    sens_l,
    specificity,
    n_control,
    n_screen
) {
  # --------------------------------------------------
  # Description:
  # Main projection function. Generates trial outcomes for the standard, IE,
  # in both the control and screen arms. Results are computed
  # for early-stage, late-stage, and overall cancer outcomes.
  #
  # Inputs:
  # - start_age: baseline age at trial entry
  # - numscreens: number of screening examinations
  # - screen_int: screening interval
  # - num_followup_intervals: amount of follow-up after the final screen
  # - rate_matrix: transition rate matrix
  # - sens_e: early-stage sensitivity
  # - sens_l: late-stage sensitivity
  # - specificity: test specificity
  # - n_control: control-arm sample size
  # - n_screen: screen-arm sample size
  #
  # Outputs:
  # - List containing params, timeline, control-arm outcomes, screen-arm outcomes,
  #   and stage-shift summaries for standard and IE designs
  # --------------------------------------------------
  screen_times <- seq(0, (numscreens - 1) * screen_int, by = screen_int)
  
  final_time <- max(screen_times) + num_followup_intervals
  
  trial_times <- seq(0, final_time, by = screen_int)
  
  if (tail(trial_times, 1) < final_time) {
    trial_times <- c(trial_times, final_time)
  }
  
  trial_times <- sort(unique(trial_times))
  
  interval_starts <- head(trial_times, -1)
  interval_ends <- tail(trial_times, -1)
  
  params <- list(
    start_age = start_age,
    numscreens = numscreens,
    screen_int = screen_int,
    num_followup_intervals = num_followup_intervals,
    screen_times = screen_times,
    trial_times = trial_times,
    interval_starts = interval_starts,
    interval_ends = interval_ends,
    final_time = final_time,
    rate_matrix = rate_matrix,
    sens_e = sens_e,
    sens_l = sens_l,
    specificity = specificity,
    n_control = n_control,
    n_screen = n_screen
  )
  
  timeline <- data.frame(
    time_since_randomization = interval_ends,
    age1 = start_age + interval_starts,
    age2 = start_age + interval_ends,
    period = ifelse(interval_starts < max(screen_times), "screen", "postscreen")
  )
  
  designs <- c("standard", "IE")
  
  control <- setNames(vector("list", length(designs)), designs)
  screen <- setNames(vector("list", length(designs)), designs)
  
  for (d in designs) {
    control[[d]] <- list(
      early = get_control_arm_trial_results("early", d, params),
      late = get_control_arm_trial_results("late", d, params)
    )
    
    screen[[d]] <- list(
      early = get_screen_arm_trial_results("early", d, params),
      late = get_screen_arm_trial_results("late", d, params)
    )
  }
  
  add_control_overall <- function(x) {
    x$overall <- list(
      incidence = list(
        yearly = list(
          control_arm_detected =
            x$early$incidence$yearly$control_arm_detected +
            x$late$incidence$yearly$control_arm_detected
        ),
        cumulative = list(
          control_arm_detected =
            x$early$incidence$cumulative$control_arm_detected +
            x$late$incidence$cumulative$control_arm_detected
        )
      ),
      count = list(
        yearly = list(
          control_arm_detected =
            x$early$count$yearly$control_arm_detected +
            x$late$count$yearly$control_arm_detected
        ),
        cumulative = list(
          control_arm_detected =
            x$early$count$cumulative$control_arm_detected +
            x$late$count$cumulative$control_arm_detected
        )
      )
    )
    x
  }
  
  add_screen_overall <- function(x) {
    x$overall <- list(
      incidence = list(
        yearly = list(
          screen_detected =
            x$early$incidence$yearly$screen_detected +
            x$late$incidence$yearly$screen_detected,
          interval_detected =
            x$early$incidence$yearly$interval_detected +
            x$late$incidence$yearly$interval_detected,
          screen_arm_detected =
            x$early$incidence$yearly$screen_arm_detected +
            x$late$incidence$yearly$screen_arm_detected
        ),
        cumulative = list(
          screen_detected =
            x$early$incidence$cumulative$screen_detected +
            x$late$incidence$cumulative$screen_detected,
          interval_detected =
            x$early$incidence$cumulative$interval_detected +
            x$late$incidence$cumulative$interval_detected,
          screen_arm_detected =
            x$early$incidence$cumulative$screen_arm_detected +
            x$late$incidence$cumulative$screen_arm_detected
        )
      ),
      count = list(
        yearly = list(
          screen_detected =
            x$early$count$yearly$screen_detected +
            x$late$count$yearly$screen_detected,
          interval_detected =
            x$early$count$yearly$interval_detected +
            x$late$count$yearly$interval_detected,
          screen_arm_detected =
            x$early$count$yearly$screen_arm_detected +
            x$late$count$yearly$screen_arm_detected
        ),
        cumulative = list(
          screen_detected =
            x$early$count$cumulative$screen_detected +
            x$late$count$cumulative$screen_detected,
          interval_detected =
            x$early$count$cumulative$interval_detected +
            x$late$count$cumulative$interval_detected,
          screen_arm_detected =
            x$early$count$cumulative$screen_arm_detected +
            x$late$count$cumulative$screen_arm_detected
        )
      )
    )
    x
  }
  
  for (d in designs) {
    control[[d]] <- add_control_overall(control[[d]])
    screen[[d]] <- add_screen_overall(screen[[d]])
  }
  
  get_stage_shift <- function(design) {
    control_late_yearly <- control[[design]]$late$incidence$yearly$control_arm_detected
    screen_late_yearly <- screen[[design]]$late$incidence$yearly$screen_arm_detected
    
    control_late_cumulative <- control[[design]]$late$incidence$cumulative$control_arm_detected
    screen_late_cumulative <- screen[[design]]$late$incidence$cumulative$screen_arm_detected
    
    list(
      yearly = ifelse(
        control_late_yearly > 0,
        100 * (control_late_yearly - screen_late_yearly) / control_late_yearly,
        NA_real_
      ),
      cumulative = ifelse(
        control_late_cumulative > 0,
        100 * (control_late_cumulative - screen_late_cumulative) / control_late_cumulative,
        NA_real_
      )
    )
  }
  
  stage_shift <- setNames(vector("list", length(designs)), designs)
  
  for (d in designs) {
    stage_shift[[d]] <- get_stage_shift(d)
  }
  
  round_nested <- function(x, digits) {
    if (is.numeric(x)) {
      round(x, digits)
    } else if (is.list(x)) {
      lapply(x, round_nested, digits = digits)
    } else {
      x
    }
  }
  
  round_count_only <- function(x) {
    if (!is.list(x)) return(x)
    
    for (nm in names(x)) {
      if (nm == "count") {
        x[[nm]] <- round_nested(x[[nm]], 0)
      } else if (is.list(x[[nm]])) {
        x[[nm]] <- round_count_only(x[[nm]])
      }
    }
    
    x
  }
  
  control <- round_count_only(control)
  screen <- round_count_only(screen)
  stage_shift <- round_nested(stage_shift, 0)
  
  list(
    params = params,
    timeline = timeline,
    control = control,
    screen = screen,
    stage_shift = stage_shift
  )
}

#' Extract a stage-shift table for one design
#'
#' Converts the results object from [get_trial_results_by_design()] into
#' a table for a selected design, time scale, and outcome scale, including
#' control- and screen-arm outcomes and the implied stage shift.
#'
#' @param results Results object returned by [get_trial_results_by_design()].
#' @param design One of `"standard"` or `"IE"`.
#' @param yearly_or_cumulative Report `"yearly"` or `"cumulative"` values.
#' @param incidence_or_count Report `"incidence"` (rates) or `"count"`.
#' @param early_stage Logical; include early-stage outcomes.
#' @param overall Logical; include overall (any-stage) outcomes.
#'
#' @return A data frame with time, age interval, control/screen outcomes, and
#'   the stage shift.
#' @export
stage_shift_by_design <- function(
    results,
    design = c("standard", "IE"),
    yearly_or_cumulative = c("yearly", "cumulative"),
    incidence_or_count = c("incidence", "count"),
    early_stage = FALSE,
    overall = FALSE
) {
  # --------------------------------------------------
  # Description:
  # Converts the nested results object from get_trial_results_by_design() into
  # a formatted table for a selected design, time scale, and outcome scale.
  #
  # Inputs:
  # - results: results object returned by get_trial_results_by_design()
  # - design: "standard" or "IE"
  # - yearly_or_cumulative: "yearly" or "cumulative"
  # - incidence_or_count: "incidence" or "count"
  # - early_stage: whether to include early-stage outcomes
  # - overall: whether to include overall outcomes
  #
  # Outputs:
  # - Data frame with time, age interval, control/screen outcomes, and stage shift
  # --------------------------------------------------
  design <- match.arg(design)
  yearly_or_cumulative <- match.arg(yearly_or_cumulative)
  incidence_or_count <- match.arg(incidence_or_count)
  
  get_stage_cols <- function(stage_name, suffix = stage_name) {
    data.frame(
      setNames(
        list(results$control[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$control_arm_detected),
        paste0("control_arm_detected_", suffix)
      ),
      setNames(
        list(results$screen[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$screen_detected),
        paste0("screen_detected_", suffix)
      ),
      setNames(
        list(results$screen[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$interval_detected),
        paste0("interval_detected_", suffix)
      ),
      setNames(
        list(results$screen[[design]][[stage_name]][[incidence_or_count]][[yearly_or_cumulative]]$screen_arm_detected),
        paste0("screen_arm_detected_", suffix)
      )
    )
  }
  
  out <- data.frame(
    time_since_randomization = results$timeline$time_since_randomization,
    age1 = results$timeline$age1,
    age2 = results$timeline$age2,
    period = results$timeline$period
  )
  
  out <- cbind(out, get_stage_cols("late", "late"))
  
  if (isTRUE(early_stage)) {
    out <- cbind(out, get_stage_cols("early", "early"))
  }
  
  if (isTRUE(overall)) {
    out <- cbind(out, get_stage_cols("overall", "overall"))
  }
  
  out$stage_shift <- results$stage_shift[[design]][[yearly_or_cumulative]]
  
  out
}
