# Load and select the calibrated natural-history models provided with
# the package.

#' Load a calibrated natural-history model
#'
#' Loads a lung cancer natural-history model calibrated by maximum likelihood to
#' SEER age- and stage-specific incidence (scaled by a risk multiplier of 3.1 to
#' reflect the elevated risk of an NLST-eligible population).
#'
#' With no sojourn-time arguments, returns all the calibrated models (the
#' base-case model plus the sensitivity-analysis models and the lookup table).
#' Given `LMST` together with either `OMST` or `EMST`, returns a **single**
#' calibrated model directly — the same object as [fit_natural_history()] — so
#' you can pass `model$rate.matrix` to the projection functions without
#' re-calibrating.
#'
#' Two selection axes are supported because the models are organized that way:
#' the base case is defined by its chosen `OMST` and `LMST` (4.0, 1.35), while
#' the sensitivity-analysis models are organized by a target early-stage mean
#' sojourn time (`EMST`) and `LMST`, with `OMST` matched to reproduce that EMST.
#' Select the base case by `OMST` + `LMST`, and the sensitivity models by
#' `EMST` + `LMST`. Use [fit_natural_history()] for other combinations.
#'
#' @param cancer Cancer site. Currently only `"lung"` is available.
#' @param OMST,EMST Overall (`OMST`) or early-stage (`EMST`) mean
#'   sojourn time used to select a single model. Give exactly one of them,
#'   together with `LMST`.
#' @param LMST Late-stage mean sojourn time; required to select a model.
#'
#' @return If no sojourn times are given, a list with:
#'   \describe{
#'     \item{base_fit}{The base-case model; use `base_fit$rate.matrix`.}
#'     \item{base_omst, base_lmst, base_emst}{Base-case OMST (4.0, chosen), LMST
#'       (1.35, chosen), and EMST (about 3.01, derived by the fit).}
#'     \item{metadata}{Metadata for the sensitivity-analysis models (`fitID`,
#'       OMST, DMST).}
#'     \item{fits}{List of calibrated models over the sensitivity analyses.}
#'     \item{matched_pairs}{Lookup mapping a target EMST and LMST to the matched
#'       OMST / `fitID`.}
#'   }
#'   If a model is selected, that single model object (with `rate.matrix` and
#'   `summary_out`), as returned by [fit_natural_history()].
#'
#' @examples
#' # All calibrated models
#' all_models <- load_fitted_model("lung")
#' c(OMST = all_models$base_omst, LMST = all_models$base_lmst,
#'   EMST = all_models$base_emst)
#'
#' # Select the base-case model directly (instant; no re-calibration)
#' model <- load_fitted_model("lung", OMST = 4.0, LMST = 1.35)
#' dim(model$rate.matrix)
#'
#' # Select a sensitivity-analysis model by target EMST + LMST
#' m2 <- load_fitted_model("lung", EMST = 4.0, LMST = 1.35)
#' m2$summary_out$OMST   # the matched OMST
#' @export
load_fitted_model <- function(cancer = "lung", OMST = NULL, LMST = NULL, EMST = NULL) {
  cancer <- match.arg(cancer, choices = "lung")

  rdata <- system.file("extdata", paste0(cancer, "_outfile.Rdata"),
                       package = "mbiedesign")
  if (!nzchar(rdata)) {
    stop("Fitted model file for '", cancer, "' not found in the package.")
  }

  e <- new.env(parent = emptyenv())
  load(rdata, envir = e)

  mp_file <- system.file("extdata", "matched_pairs.csv", package = "mbiedesign")
  matched_pairs <- if (nzchar(mp_file)) read.csv(mp_file) else NULL

  models <- list(
    base_fit      = e$base_fit,
    base_omst     = e$base_fit$summary_out$OMST,
    base_lmst     = e$base_fit$summary_out$LMST,
    base_emst     = e$target_base_emst,
    metadata      = e$metadata_lung,
    fits          = e$lungout,
    matched_pairs = matched_pairs
  )

  if (is.null(OMST) && is.null(LMST) && is.null(EMST)) {
    return(models)
  }
  select_fitted_model(models, OMST = OMST, LMST = LMST, EMST = EMST)
}

# Select a single calibrated model by (OMST, LMST) or (EMST, LMST). Internal.
select_fitted_model <- function(models, OMST = NULL, LMST = NULL, EMST = NULL,
                                tol = 0.02) {
  if (is.null(LMST)) {
    stop("Give LMST together with OMST or EMST to select a model.", call. = FALSE)
  }
  if (!is.null(OMST) && !is.null(EMST)) {
    stop("Select by either OMST or EMST (with LMST), not both.", call. = FALSE)
  }

  if (!is.null(OMST)) {
    if (abs(OMST - models$base_omst) <= tol && abs(LMST - models$base_lmst) <= tol) {
      return(models$base_fit)
    }
    md  <- models$metadata
    hit <- which(abs(md$OMST - OMST) <= tol & abs(md$DMST - LMST) <= tol)
    if (length(hit) == 0) {
      stop("No calibrated model for OMST = ", OMST, ", LMST = ", LMST,
           ". Available (OMST, LMST): ",
           paste0("(", round(c(models$base_omst, md$OMST), 2), ", ",
                  round(c(models$base_lmst, md$DMST), 2), ")", collapse = ", "),
           ". Use fit_natural_history() for other values.", call. = FALSE)
    }
    return(models$fits[[md$fitID[hit[1]]]])
  }

  mp  <- models$matched_pairs
  hit <- which(abs(mp$target_EMST - EMST) <= tol & abs(mp$LMST - LMST) <= tol)
  if (length(hit) == 0) {
    stop("No calibrated model for EMST = ", EMST, ", LMST = ", LMST,
         ". Available (EMST, LMST): ",
         paste0("(", round(mp$target_EMST, 2), ", ", round(mp$LMST, 2), ")",
                collapse = ", "),
         ". Use fit_natural_history() for other values.", call. = FALSE)
  }
  fid <- mp$fitID[hit[1]]
  if (is.na(fid)) models$base_fit else models$fits[[fid]]
}
