# Achieved power for a given sample size

Inverts a one-argument required-sample-size function to find the power
at which the required sample size equals a current sample size. Useful
for reporting the power a fixed trial achieves.

## Usage

``` r
get_achieved_power(n_required_fun, N_current, interval = c(0.001, 0.999))
```

## Arguments

- n_required_fun:

  A function of a single argument (power) returning the required sample
  size, e.g.
  [`n_required_by_design()`](https://kemalgog.github.io/mbiedesign/reference/n_required_by_design.md)
  with its other arguments held fixed.

- N_current:

  The available (fixed) sample size.

- interval:

  Search interval for the power root.

## Value

The achieved power rounded to 3 digits, or the strings `">0.999"` /
`"<0.001"` when the power lies outside the search interval.
