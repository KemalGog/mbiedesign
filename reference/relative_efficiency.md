# Relative efficiency of the intended-effect vs traditional design

Ratio of the sample size required by the traditional design to that
required by the intended-effect (IE) design for the same event
probabilities. Values greater than 1 mean the IE design needs fewer
participants.

## Usage

``` r
relative_efficiency(
  p_c_std,
  p_s_std,
  p_c_ie,
  p_s_ie,
  pi_pos,
  alpha = 0.025,
  power = 0.9,
  R = 1
)
```

## Arguments

- p_c_std, p_s_std:

  Traditional control- and screen-arm event probabilities.

- p_c_ie, p_s_ie:

  IE control- and screen-arm event probabilities.

- pi_pos:

  Ever-positive fraction.

- alpha:

  One-sided significance level.

- power:

  Target power.

- R:

  Allocation ratio, screen / control.

## Value

The relative efficiency (a single numeric value).
