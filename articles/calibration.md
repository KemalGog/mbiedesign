# Calibrating the natural-history model

Every projection in `mbiedesign` starts from a **calibrated
natural-history model** — a continuous-time multistate model whose
transition rates are fit by maximum likelihood to age- and
stage-specific cancer incidence. This vignette shows the full path from
raw incidence data to a fitted model.

For convenience,
[`fit_natural_history()`](https://kemalgog.github.io/mbiedesign/reference/fit_natural_history.md)
and
[`load_fitted_model()`](https://kemalgog.github.io/mbiedesign/reference/load_fitted_model.md)
do the data handling for you. Here we open that up so you can see
exactly what happens — and adapt it to your own data.

``` r

library(mbiedesign)
```

## 1. The incidence data

The package includes SEER-style lung cancer incidence, one row per age
band and stage:

``` r

xlsx <- system.file("extdata", "lung_incidence.xlsx", package = "mbiedesign")
inc  <- openxlsx::read.xlsx(xlsx)
head(inc)
#>           age Rate Count      Pop    stage years midage       PY         rate
#> 1    00 years  0.0     0  5732775 Advanced     1    0.5  5732775 0.000000e+00
#> 2 01-04 years  0.0     0 23147084 Advanced     1    3.0 23147084 0.000000e+00
#> 3 05-09 years  0.0     0 28781431 Advanced     1    7.5 28781431 0.000000e+00
#> 4 10-14 years  0.0     0 29521104 Advanced     1   12.5 29521104 0.000000e+00
#> 5 15-19 years  0.0     7 31007881 Advanced     1   17.5 31007881 2.257491e-07
#> 6 20-24 years  0.1    30 30673778 Advanced     1   22.5 30673778 9.780341e-07
```

SEER\*Stat gives incidence by 5-year age band, with case counts,
populations, and rates. The model needs five columns (one row per age
band and stage):

- **`midage`** — midpoint of the age band, in years (e.g. 62.5 for
  60–64),
- **`Count`** — number of incident cancers in that age band and stage,
- **`Pop`** — population (denominator),
- **`years`** — person-time multiplier used to form the denominator
  (person-years = `Pop * years`),
- **`stage`** — stage at diagnosis, in two levels: **`"Early"`** for
  early-stage rows and **`"Advanced"`** for late-stage rows.

**The early vs advanced split is your modeling choice.** Label
early-stage rows `"Early"` and late-stage rows `"Advanced"`; *which*
stages count as early is up to you. In this paper we used AJCC7 stage
I–II as early and III–IV as advanced, but you might instead use SEER
Summary Stage, AJCC6 or AJCC9, or a different cut point (for example
only stage I as early). How you handle unstaged or missing-stage records
is likewise your choice — in this paper we distributed them across early
and advanced in proportion to the observed (known-stage) distribution.
Any extra columns (here `age`, `Rate`) are ignored. To calibrate a
**different cancer site**, export the same age-by-stage incidence from
SEER\*Stat and supply it via `the_data` (see the last section).

## 2. Preparing the data

Two derived quantities and one filter turn the raw counts into what the
likelihood uses:

``` r

risk_multiplier <- 3.1   # incidence inflation factor (see below)

inc$PY    <- inc$Pop * inc$years                      # person-years
inc$rate  <- risk_multiplier * inc$Count / inc$PY * 1e5   # rate per 100,000
inc$Count <- inc$Count * risk_multiplier              # inflate the counts too
inc <- inc[inc$midage > 0 & inc$midage < 80, ]        # restrict to ages 0-80

head(inc[c("midage", "stage", "Count", "PY", "rate")])
#>   midage    stage Count       PY       rate
#> 1    0.5 Advanced   0.0  5732775 0.00000000
#> 2    3.0 Advanced   0.0 23147084 0.00000000
#> 3    7.5 Advanced   0.0 28781431 0.00000000
#> 4   12.5 Advanced   0.0 29521104 0.00000000
#> 5   17.5 Advanced  21.7 31007881 0.06998221
#> 6   22.5 Advanced  93.0 30673778 0.30319056
```

The **risk multiplier** is the one modeling choice worth dwelling on.
General- population SEER incidence understates the risk in a
*screening-trial* population, which is enriched for risk factors (in the
NLST, heavy smoking history). The multiplier of `3.1` inflates incidence
to reflect that NLST-eligible high-risk cohort. It is a parameter you
can change — a different trial population, or a different cancer, will
warrant a different value.

This is exactly what
[`fit_natural_history()`](https://kemalgog.github.io/mbiedesign/reference/fit_natural_history.md)
does internally via its `cancer` and `risk_multiplier` arguments; the
block above simply makes it explicit.

## 3. Fitting the model

With the prepared data, calibration is a single call. You specify the
two target sojourn times — the overall mean sojourn time (`OMST`) and
the late-stage mean sojourn time (`LMST`) — and the fitter finds
transition rates consistent with them that best reproduce the observed
incidence:

``` r

# Equivalent one-liner (loads + prepares the provided data for you):
fit <- fit_natural_history(OMST = 4.0, LMST = 1.35)

# ...or fit the data we prepared above explicitly:
fit <- fit_natural_history(OMST = 4.0, LMST = 1.35, the_data = inc)
```

The optimizer arguments are optional with sensible defaults:

- **`k`** (default 16) — number of states in the model,
- **`num_seeds`** (default 1) — random optimizer restarts; the best fit
  is kept,
- **`mean1`** (default 3) — mean of the distribution used to draw
  starting values,
- **`risk_multiplier`** (default 3.1) — the inflation factor from
  Section 2,
- **`the_data`** — supply your own incidence frame to override the
  provided `cancer`.

Fitting runs a maximum-likelihood optimization, so it takes a little
time; the chunk above is not evaluated when this vignette builds. The
package provides the result of exactly this calibration.

## 4. Inspecting the fitted model

Load the provided fit and look at what calibration produced:

``` r

fit <- load_fitted_model("lung")$base_fit
#> Warning: namespace 'mbtrialdesign' is not available and has been replaced
#> by .GlobalEnv when processing object 'base_fit'

# The fitted transition-rate matrix (the object projections consume)
dim(fit$rate.matrix)
#> [1] 16 16

# Sojourn-time summaries, including the early-stage mean sojourn time (EMST)
str(fit$summary_out)
#> List of 6
#>  $ a           : num 0.0876
#>  $ b           : num 0.245
#>  $ sojourn_time: num 4
#>  $ OMST        : num 4
#>  $ LMST        : num 1.35
#>  $ EMST        : num 3.01
```

The `rate.matrix` is the hand-off point: pass it to
[`get_trial_results_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_trial_results_by_design.md)
or
[`get_summary_outcomes_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_summary_outcomes_by_design.md)
to project trial outcomes (see
[`vignette("mbiedesign")`](https://kemalgog.github.io/mbiedesign/articles/mbiedesign.md)).

## Fitting your own data or another cancer site

The provided model is lung, but the natural-history model calibrates to
**any cancer site**. Prepare an incidence data frame with the five
columns from Section 1 — `midage`, `Count`, `Pop`, `years`, `stage`
(`"Early"` vs `"Advanced"`) — typically a SEER\*Stat age-by-stage
incidence export for your site, then pass it directly:

``` r

# OMST and LMST are your choice; use sojourn times reasonable for your cancer,
# not the lung values below.
fit <- fit_natural_history(OMST = 4.0, LMST = 1.35, the_data = my_incidence)
```

The `OMST` and `LMST` targets are modeling inputs *you* choose —
reasonable values depend on the cancer and should come from prior
estimates or the literature for your site; the 4.0 / 1.35 used here are
lung-specific. No `risk_multiplier` is applied to data you supply
yourself — scale it as appropriate for your target population before
fitting.
