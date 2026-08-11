# Summary outcomes and sample size by trial design

The main function of the package. Given a calibrated natural-history
`rate_matrix` and a screening schedule, it projects the end-of-trial
late-stage outcomes and computes the required per-arm sample size and
relative efficiency for the traditional and intended-effect (IE)
designs, returning one row per design (a summary table). For the full
time-dependent early/late/overall outcomes, use
[`get_trial_results_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_trial_results_by_design.md).

## Usage

``` r
get_summary_outcomes_by_design(
  rate_matrix,
  start_age,
  numscreens,
  screen_int,
  num_followup_intervals,
  sens_e,
  sens_l,
  specificity,
  design = c("traditional", "IE"),
  allocation_ratio = 1,
  alpha = 0.025,
  power = 0.9
)
```

## Arguments

- rate_matrix:

  Natural-history transition rate matrix (e.g. from
  [`load_fitted_model()`](https://kemalgog.github.io/mbiedesign/reference/load_fitted_model.md)
  or
  [`fit_natural_history()`](https://kemalgog.github.io/mbiedesign/reference/fit_natural_history.md)).

- start_age:

  Age at trial entry.

- numscreens:

  Number of screening examinations.

- screen_int:

  Interval between screens (years).

- num_followup_intervals:

  Follow-up after the final screen.

- sens_e:

  Early-stage test sensitivity.

- sens_l:

  Late-stage test sensitivity.

- specificity:

  Test specificity.

- design:

  Which design(s) to return: `"traditional"`, `"IE"`, or both. Defaults
  to both.

- allocation_ratio:

  Screen-arm to control-arm allocation ratio (n_screen / n_control); 1
  is balanced (the default). Event rates are per-participant, so only
  the allocation ratio affects the required sample size — absolute arm
  sizes are not needed to estimate it.

- alpha:

  One-sided significance level.

- power:

  Target power (= 1 - beta).

## Value

A data frame with one row per requested `design` and columns:

- design:

  `"traditional"` or `"IE"`.

- p_c, p_s:

  Control- and screen-arm late-stage event probabilities (per randomized
  participant for traditional; per ever-positive for IE).

- rr, rd:

  Risk ratio (`p_s / p_c`) and risk difference (`p_c - p_s`).

- fraction:

  The subgroup fraction the design conditions on: `NA` for traditional,
  the ever-positive fraction (pi) for IE.

- n_per_arm:

  Required per-arm sample size.

- relative_efficiency:

  N_traditional / N_design (1 for traditional).
