# Load a calibrated natural-history model

Loads a lung cancer natural-history model calibrated by maximum
likelihood to SEER age- and stage-specific incidence (scaled by a risk
multiplier of 3.1 to reflect the elevated risk of an NLST-eligible
population).

## Usage

``` r
load_fitted_model(cancer = "lung", OMST = NULL, LMST = NULL, EMST = NULL)
```

## Arguments

- cancer:

  Cancer site. Currently only `"lung"` is available.

- OMST, EMST:

  Overall (`OMST`) or early-stage (`EMST`) mean sojourn time used to
  select a single model. Give exactly one of them, together with `LMST`.

- LMST:

  Late-stage mean sojourn time; required to select a model.

## Value

If no sojourn times are given, a list with:

- base_fit:

  The base-case model; use `base_fit$rate.matrix`.

- base_omst, base_lmst, base_emst:

  Base-case OMST (4.0, chosen), LMST (1.35, chosen), and EMST (about
  3.01, derived by the fit).

- metadata:

  Metadata for the sensitivity-analysis models (`fitID`, OMST, DMST).

- fits:

  List of calibrated models over the sensitivity analyses.

- matched_pairs:

  Lookup mapping a target EMST and LMST to the matched OMST / `fitID`.

If a model is selected, that single model object (with `rate.matrix` and
`summary_out`), as returned by
[`fit_natural_history()`](https://kemalgog.github.io/mbiedesign/reference/fit_natural_history.md).

## Details

With no sojourn-time arguments, returns all the calibrated models (the
base-case model plus the sensitivity-analysis models and the lookup
table). Given `LMST` together with either `OMST` or `EMST`, returns a
**single** calibrated model directly — the same object as
[`fit_natural_history()`](https://kemalgog.github.io/mbiedesign/reference/fit_natural_history.md)
— so you can pass `model$rate.matrix` to the projection functions
without re-calibrating.

Two selection axes are supported because the models are organized that
way: the base case is defined by its chosen `OMST` and `LMST` (4.0,
1.35), while the sensitivity-analysis models are organized by a target
early-stage mean sojourn time (`EMST`) and `LMST`, with `OMST` matched
to reproduce that EMST. Select the base case by `OMST` + `LMST`, and the
sensitivity models by `EMST` + `LMST`. Use
[`fit_natural_history()`](https://kemalgog.github.io/mbiedesign/reference/fit_natural_history.md)
for other combinations.

## Examples

``` r
# All calibrated models
all_models <- load_fitted_model("lung")
#> Warning: namespace ‘mbtrialdesign’ is not available and has been replaced
#> by .GlobalEnv when processing object ‘base_fit’
c(OMST = all_models$base_omst, LMST = all_models$base_lmst,
  EMST = all_models$base_emst)
#>    OMST    LMST    EMST 
#> 4.00000 1.35000 3.00525 

# Select the base-case model directly (instant; no re-calibration)
model <- load_fitted_model("lung", OMST = 4.0, LMST = 1.35)
#> Warning: namespace ‘mbtrialdesign’ is not available and has been replaced
#> by .GlobalEnv when processing object ‘base_fit’
dim(model$rate.matrix)
#> [1] 16 16

# Select a sensitivity-analysis model by target EMST + LMST
m2 <- load_fitted_model("lung", EMST = 4.0, LMST = 1.35)
#> Warning: namespace ‘mbtrialdesign’ is not available and has been replaced
#> by .GlobalEnv when processing object ‘base_fit’
m2$summary_out$OMST   # the matched OMST
#> [1] 5.1
```
