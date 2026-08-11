# Primary-test emission matrix (internal HMM helper)

Exposed for advanced use; most users will not call this directly.

## Usage

``` r
get_emission_test(k, sens_e, sens_l, specificity)
```

## Arguments

- k:

  Number of states in the model.

- sens_e:

  Early-stage test sensitivity.

- sens_l:

  Late-stage test sensitivity.

- specificity:

  Test specificity.

## Value

The primary-test emission probability matrix.
