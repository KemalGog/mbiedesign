# Calibrate the natural-history model to incidence data

Fits the continuous-time multistate natural-history model of Lange et
al. (2024) to age- and stage-specific incidence data by maximum
likelihood, targeting a specified overall mean sojourn time (`OMST`) and
late-stage mean sojourn time (`LMST`). By default it loads and
calibrates to the provided incidence data for `cancer`; supply
`the_data` to fit your own. The returned object supplies the
`rate.matrix` that the projection functions (e.g.
[`get_trial_results_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_trial_results_by_design.md))
consume. See the calibration vignette
([`vignette("calibration", package = "mbiedesign")`](https://kemalgog.github.io/mbiedesign/articles/calibration.md))
for the data-prep steps done internally.

## Usage

``` r
fit_natural_history(
  OMST,
  LMST,
  cancer = "lung",
  k = 16,
  num_seeds = 1,
  mean1 = 3,
  risk_multiplier = 3.1,
  the_data = NULL
)
```

## Arguments

- OMST:

  Target overall mean sojourn time, in years.

- LMST:

  Target late-stage mean sojourn time, in years.

- cancer:

  Cancer site whose incidence data to fit (currently `"lung"`). Ignored
  when `the_data` is supplied.

- k:

  Number of states in the CTMC.

- num_seeds:

  Number of random optimizer restarts; the best fit is kept.

- mean1:

  Mean of the normal distribution used to draw starting values.

- risk_multiplier:

  Incidence inflation factor applied to the provided data, e.g. `3.1` to
  reflect an NLST-eligible high-risk population. Ignored when `the_data`
  is supplied.

- the_data:

  Optional incidence data frame (columns `midage`, `Count`, `Pop`,
  `years`, `stage`) to fit instead of the provided `cancer` data. Used
  as-is, with no `risk_multiplier` applied.

## Value

A list with the raw `optim` fit (`the_fit`), an optional diagnostic plot
(`outplot`, `NULL` if ggplot2/reshape2 are unavailable), the fitted
`rate.matrix`, and a `summary_out` list of sojourn-time quantities
(`OMST`, `LMST`, and the derived `EMST`, the early-stage mean sojourn
time).

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit the lung base case: loads provided incidence, applies the 3.1 multiplier
fit <- fit_natural_history(OMST = 4.0, LMST = 1.35)
fit$summary_out$EMST  # early-stage mean sojourn time (derived)
} # }
```
