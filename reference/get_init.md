# Initial state distribution (internal HMM helper)

Exposed for advanced use (e.g. custom pathway decompositions); most
users will not call this directly.

## Usage

``` r
get_init(rate_matrix, start_age)
```

## Arguments

- rate_matrix:

  Natural-history transition rate matrix.

- start_age:

  Age at trial entry.

## Value

A vector of initial state-occupancy probabilities.
