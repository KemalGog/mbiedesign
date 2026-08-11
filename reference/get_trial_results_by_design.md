# Project trial outcomes for every design

Main projection function. Given a calibrated natural-history model and a
screening schedule, projects control- and screen-arm cancer outcomes
under the traditional (`"standard"`) and intended-effect (`"IE"`)
designs, for early-stage, late-stage, and overall cancer.

## Usage

``` r
get_trial_results_by_design(
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
)
```

## Arguments

- start_age:

  Age at trial entry.

- numscreens:

  Number of screening examinations.

- screen_int:

  Interval between screens (years).

- num_followup_intervals:

  Length of follow-up after the final screen.

- rate_matrix:

  Natural-history transition rate matrix (e.g. from a fitted model via
  [`load_fitted_model()`](https://kemalgog.github.io/mbiedesign/reference/load_fitted_model.md)).

- sens_e:

  Early-stage test sensitivity.

- sens_l:

  Late-stage test sensitivity.

- specificity:

  Test specificity.

- n_control:

  Control-arm sample size.

- n_screen:

  Screen-arm sample size.

## Value

A list with the projection `params`, the trial timeline, and control-arm
and screen-arm outcomes plus stage-shift summaries for the standard and
IE designs.

## See also

[`stage_shift_by_design()`](https://kemalgog.github.io/mbiedesign/reference/stage_shift_by_design.md)
to extract a formatted table for one design.
