# Package index

## Main function

Summary outcomes and per-arm sample size by trial design.

- [`get_summary_outcomes_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_summary_outcomes_by_design.md)
  : Summary outcomes and sample size by trial design

## Calibration

Fit the multi-state natural-history model to incidence data, or load a
calibrated model.

- [`fit_natural_history()`](https://kemalgog.github.io/mbiedesign/reference/fit_natural_history.md)
  : Calibrate the natural-history model to incidence data
- [`load_fitted_model()`](https://kemalgog.github.io/mbiedesign/reference/load_fitted_model.md)
  : Load a calibrated natural-history model

## Projection

Project control- and screen-arm outcomes by trial design.

- [`get_trial_results_by_design()`](https://kemalgog.github.io/mbiedesign/reference/get_trial_results_by_design.md)
  : Project trial outcomes for every design
- [`stage_shift_by_design()`](https://kemalgog.github.io/mbiedesign/reference/stage_shift_by_design.md)
  : Extract a stage-shift table for one design
- [`get_positivity_rate()`](https://kemalgog.github.io/mbiedesign/reference/get_positivity_rate.md)
  : Ever-positive fraction under the intended-effect design

## Sample size and efficiency

Size or power a single design, or compute relative efficiency.

- [`n_required_by_design()`](https://kemalgog.github.io/mbiedesign/reference/n_required_by_design.md)
  : Required sample size for a target power, by design
- [`power_by_design()`](https://kemalgog.github.io/mbiedesign/reference/power_by_design.md)
  : Power for a given sample size, by design
- [`get_achieved_power()`](https://kemalgog.github.io/mbiedesign/reference/get_achieved_power.md)
  : Achieved power for a given sample size
- [`relative_efficiency()`](https://kemalgog.github.io/mbiedesign/reference/relative_efficiency.md)
  : Relative efficiency of the intended-effect vs traditional design
