# Required sample size for a target power, by design

Computes the total sample size required to achieve a target power for
the traditional or intended-effect (`"IE"`) design. As in
[`power_by_design()`](https://kemalgog.github.io/mbiedesign/reference/power_by_design.md),
the IE design inflates the variance by `1 / pi_pos`.

## Usage

``` r
n_required_by_design(
  power = 0.8,
  p_c,
  p_s,
  design = c("traditional", "IE"),
  alpha = 0.025,
  R = 1,
  pi_pos = NULL
)
```

## Arguments

- power:

  Target power.

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

The required total sample size (a single numeric value).

## See also

[`power_by_design()`](https://kemalgog.github.io/mbiedesign/reference/power_by_design.md),
[`get_achieved_power()`](https://kemalgog.github.io/mbiedesign/reference/get_achieved_power.md)
