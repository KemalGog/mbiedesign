# Forward algorithm (internal HMM helper)

Exposed for advanced use (e.g. custom pathway decompositions); most
users will not call this directly.

## Usage

``` r
forward_prob(
  obs,
  times,
  delta,
  rate_matrix,
  emission_matrix,
  return_state = FALSE,
  start_time = NULL
)
```

## Arguments

- obs:

  Observation sequence.

- times:

  Observation times.

- delta:

  Initial state distribution.

- rate_matrix:

  Natural-history transition rate matrix.

- emission_matrix:

  Emission probability matrix.

- return_state:

  Logical; if `TRUE`, return the state-probability vector.

- start_time:

  Optional start time for the recursion.

## Value

The forward probability (or state vector if `return_state = TRUE`).
