# mbiedesign

**Model-based power and sample-size evaluation of cancer screening trial
designs with late-stage incidence as an endpoint.**

`mbiedesign` projects the late-stage outcomes and estimates per-arm
sample size required to detect a reduction in late-stage incidence under
two trial designs:

- **Traditional** — compares late-stage incidence across all randomized
  participants.
- **Intended-effect (IE)** — restricts the analysis to *ever-positive*
  participants.

The function,
[`get_summary_outcomes_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_summary_outcomes_by_design.md),
returns the projected outcomes, the required per-arm sample size, and
the relative efficiency of the IE versus the traditional design.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("KemalGog/mbiedesign", build_vignettes = TRUE)
```

`build_vignettes = TRUE` installs the worked-example tutorials listed
below.

## Start here

After installing, the vignettes walk through everything from a quick
application of the model in single trial setting to reproducing the
results in the manuscript:

``` r

library(mbiedesign)

vignette("mbiedesign")           # getting started: choose a model
vignette("calibration")          # fit the natural-history model to your own data
vignette("manuscript-analysis")  # reproduce the analyses in the manuscript
```

## A minimal example

``` r

library(mbiedesign)

# A precalibrated lung model (pick the overall and late-stage mean sojourn times;
# the early-stage sojourn time is derived).
model <- load_fitted_model("lung", OMST = 4.0, LMST = 1.35)

# Trial specifications: 3 annual LDCT screens, 4.5 years of follow-up.
get_summary_outcomes_by_design(
  rate_matrix = model$rate.matrix,
  start_age = 62, numscreens = 3, screen_int = 1, num_followup_intervals = 4.5,
  sens_e = 0.35, sens_l = 0.82, specificity = 0.855,   # test performance
  power = 0.90, alpha = 0.025                           # 90% power, one-sided level
)
#>        design    p_c    p_s   rr   fraction  n_per_arm  relative_efficiency
#>   traditional 0.0150 0.0111 0.74      NA        17296                 1.00
#>            IE 0.0252 0.0151 0.60   0.38 (pi)     10796                 1.60
```

Each row carries the control- and screen-arm late-stage event rates, the
risk ratio and difference, the subgroup fraction the design conditions
on, the required per-arm sample size, and the relative efficiency
(N_traditional / N_design). See
[`vignette("mbiedesign")`](https://kemalgog.github.io/mbiedesign/articles/mbiedesign.md)
for a full walk-through.

## Package structure

| Path | Contents |
|:---|:---|
| `R/` | package functions — calibration, projection, and sample-size / power |
| `inst/extdata/` | precalibrated lung model and SEER incidence data |
| `data-raw/` | reproducible calibration driver (`fit_lung.R`) |
| `vignettes/` | getting-started, calibration, and manuscript-reproduction tutorials |
| `tests/` | unit tests |

## Reference

> Lange JM, Gogebakan KC, Katki H, Etzioni R. Power evaluation of cancer
> screening trial designs with incidence-based endpoints using a
> multi-state disease model
