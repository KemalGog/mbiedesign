# Ever-positive fraction under the intended-effect design

Computes the intended-effect (IE) ever-positive fraction: the
probability of at least one positive primary screening test across the
scheduled screening examinations, including both true-positive and
false-positive primary-test results.

## Usage

``` r
get_positivity_rate(
  sens_e,
  sens_l,
  specificity,
  screen_times,
  start_age,
  rate_matrix
)
```

## Arguments

- sens_e:

  Early-stage test sensitivity.

- sens_l:

  Late-stage test sensitivity.

- specificity:

  Test specificity.

- screen_times:

  Numeric vector of scheduled screen times.

- start_age:

  Age at trial entry.

- rate_matrix:

  Natural-history transition rate matrix (e.g. from a fitted model via
  [`load_fitted_model()`](https://kemalgog.github.io/mbiedesign/reference/load_fitted_model.md)).

## Value

A single numeric value: the ever-positive probability.
