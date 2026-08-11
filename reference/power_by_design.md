# Power for a given sample size, by design

Computes the power to detect a difference in event probabilities for the
traditional or intended-effect (`"IE"`) design. The IE design inflates
the variance by `1 / pi_pos`, reflecting that inference is restricted to
the ever-positive subgroup.

## Usage

``` r
power_by_design(
  n_total,
  p_c,
  p_s,
  design = c("traditional", "IE"),
  alpha = 0.025,
  R = 1,
  pi_pos = NULL
)
```

## Arguments

- n_total:

  Total sample size (both arms).

- p_c, p_s:

  Control- and screen-arm event probabilities (use the probabilities
  appropriate to the chosen design).

- design:

  One of `"traditional"` or `"IE"`.

- alpha:

  One-sided significance level.

- R:

  Allocation ratio, screen-arm size / control-arm size.

- pi_pos:

  Ever-positive fraction; required for `design = "IE"`.

## Value

The power (a single numeric value).

## See also

[`n_required_by_design()`](https://kemalgog.github.io/mbiedesign/reference/n_required_by_design.md)
