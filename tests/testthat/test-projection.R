test_that("load_fitted_model returns the expected structure", {
  fit <- load_fitted_model("lung")
  expect_type(fit, "list")
  expect_true(all(c("base_fit", "base_emst", "metadata", "fits",
                    "matched_pairs") %in% names(fit)))
  expect_equal(dim(fit$base_fit$rate.matrix), c(16, 16))
  expect_true(is.finite(fit$base_emst) && fit$base_emst > 0)
})

test_that("get_trial_results_by_design produces both designs and both arms", {
  fit <- load_fitted_model("lung")
  res <- get_trial_results_by_design(
    start_age = 62, numscreens = 3, screen_int = 1, num_followup_intervals = 4.5,
    rate_matrix = fit$base_fit$rate.matrix,
    sens_e = 0.35, sens_l = 0.82, specificity = 0.855,
    n_control = 26722, n_screen = 26722)
  expect_true(all(c("standard", "IE") %in% names(res$screen)))
  expect_true(all(c("standard", "IE") %in% names(res$control)))
})

test_that("base-case get_summary_outcomes_by_design reproduces the paper's headline results", {
  rm  <- load_fitted_model("lung", OMST = 4.0, LMST = 1.35)$rate.matrix
  res <- get_summary_outcomes_by_design(
    rate_matrix = rm,
    start_age = 62, numscreens = 3, screen_int = 1, num_followup_intervals = 4.5,
    sens_e = 0.35, sens_l = 0.82, specificity = 0.855,
    allocation_ratio = 1, alpha = 0.025, power = 0.90)

  # one row per design
  expect_setequal(res$design, c("traditional", "IE"))
  trad <- res[res$design == "traditional", ]
  ie   <- res[res$design == "IE", ]

  # ever-positive fraction calibrated to ~0.38 (the IE subgroup fraction)
  expect_equal(ie$fraction, 0.38, tolerance = 0.03)
  # traditional is the reference (RE == 1); IE is more efficient
  expect_equal(trad$relative_efficiency, 1)
  expect_gt(ie$relative_efficiency, 1)
  expect_lt(ie$n_per_arm, trad$n_per_arm)
  # screen arm has fewer late-stage cancers than control in each design
  expect_true(all(res$p_s < res$p_c))

  # design argument filters rows
  expect_equal(nrow(get_summary_outcomes_by_design(
    rate_matrix = rm, start_age = 62, numscreens = 3, screen_int = 1,
    num_followup_intervals = 4.5, sens_e = 0.35, sens_l = 0.82,
    specificity = 0.855, design = "IE")), 1)
})

test_that("power and sample-size calculators are mutually consistent", {
  # sample size sized for 90% power should return ~90% power
  n_total <- n_required_by_design(power = 0.90, p_c = 0.015, p_s = 0.011,
                                  design = "traditional", alpha = 0.025, R = 1)
  pw <- power_by_design(n_total = n_total, p_c = 0.015, p_s = 0.011,
                        design = "traditional", alpha = 0.025, R = 1)
  expect_equal(pw, 0.90, tolerance = 1e-3)

  # IE design requires its ever-positive fraction
  expect_error(
    n_required_by_design(power = 0.9, p_c = 0.02, p_s = 0.015, design = "IE"),
    "pi_pos"
  )
})
