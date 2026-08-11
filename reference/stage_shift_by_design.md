# Extract a stage-shift table for one design

Converts the results object from
[`get_trial_results_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_trial_results_by_design.md)
into a table for a selected design, time scale, and outcome scale,
including control- and screen-arm outcomes and the implied stage shift.

## Usage

``` r
stage_shift_by_design(
  results,
  design = c("standard", "IE"),
  yearly_or_cumulative = c("yearly", "cumulative"),
  incidence_or_count = c("incidence", "count"),
  early_stage = FALSE,
  overall = FALSE
)
```

## Arguments

- results:

  Results object returned by
  [`get_trial_results_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_trial_results_by_design.md).

- design:

  One of `"standard"` or `"IE"`.

- yearly_or_cumulative:

  Report `"yearly"` or `"cumulative"` values.

- incidence_or_count:

  Report `"incidence"` (rates) or `"count"`.

- early_stage:

  Logical; include early-stage outcomes.

- overall:

  Logical; include overall (any-stage) outcomes.

## Value

A data frame with time, age interval, control/screen outcomes, and the
stage shift.
