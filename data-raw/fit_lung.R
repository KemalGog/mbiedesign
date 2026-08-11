## data-raw/fit_lung.R
##
## Regenerates the calibrated lung natural-history model provided with the
## package:
##   inst/extdata/lung_outfile.Rdata   (base_fit, target_base_emst,
##                                       metadata_lung, lungout)
##   inst/extdata/matched_pairs.csv    (target EMST/LMST -> matched OMST)
##
## Fits the base case (OMST 4.0, LMST 1.35) plus the sensitivity-analysis
## pairs used by the paper: two LMST-sensitivity and two EMST-sensitivity
## combinations. (The factorial off-diagonal fits from the original pipeline
## are intentionally omitted.) Uses SEER-style lung incidence with a risk
## multiplier of 3.1.
##
## Run from the package root, e.g.:  Rscript data-raw/fit_lung.R
##
## Provenance note: the fits provided in inst/extdata are the exact calibration
## used in the paper (the original unified fit, restricted to the sensitivity
## pairs below). This script re-implements that calibration from scratch as a
## transparent, self-contained recipe; because optim() uses random starting
## values, a fresh run differs only negligibly from the provided fit.

set.seed(1)
library(openxlsx)
if (requireNamespace("mbiedesign", quietly = TRUE)) {
  library(mbiedesign)
} else {
  pkgload::load_all(".", quiet = TRUE)   # development: use source tree
}

k <- 16; num_seeds <- 1; mean1 <- 3; risk_multiplier <- 3.1

## ---- 1. Incidence data (lung), risk multiplier 3.1 -------------------------
xlsx <- system.file("extdata", "lung_incidence.xlsx", package = "mbiedesign")
if (!nzchar(xlsx)) xlsx <- "inst/extdata/lung_incidence.xlsx"
inc <- read.xlsx(xlsx)
inc$PY    <- inc$Pop * inc$years
inc$rate  <- risk_multiplier * inc$Count / inc$PY * 1e5
inc$Count <- inc$Count * risk_multiplier
inc <- inc[inc$midage > 0 & inc$midage < 80, , drop = FALSE]

## ---- 2. Base case exact fit (OMST 4.0, LMST 1.35) --------------------------
base_fit <- fit_natural_history(OMST = 4.0, LMST = 1.35, the_data = inc,
                                num_seeds = num_seeds, k = k, mean1 = mean1)
target_base_emst <- as.numeric(base_fit$summary_out$EMST)

## ---- 3. Paper target (target_EMST, LMST) pairs (no factorial) --------------
target_combinations <- data.frame(
  target_EMST = c(target_base_emst, target_base_emst, 1.00, 4.00),
  LMST        = c(0.50,             1.75,             1.35, 1.35)
)

## ---- 4. Search grid of OMST around each center, then fit -------------------
grids <- lapply(seq_len(nrow(target_combinations)), function(i) {
  center <- round(target_combinations$target_EMST[i] + target_combinations$LMST[i], 2)
  data.frame(OMST = round(seq(center - 0.25, center + 0.25, by = 0.1), 2),
             LMST = target_combinations$LMST[i])
})
grid <- unique(do.call(rbind, grids))
grid <- grid[grid$OMST > grid$LMST, , drop = FALSE]

metadata_lung <- data.frame(
  fitID = seq_len(nrow(grid)), site = "lung",
  OMST = grid$OMST, DMST = grid$LMST, risk_multiplier = risk_multiplier
)

lungout <- mapply(
  function(omst, dmst) fit_natural_history(OMST = omst, LMST = dmst,
                                           the_data = inc, num_seeds = num_seeds,
                                           k = k, mean1 = mean1),
  metadata_lung$OMST, metadata_lung$DMST, SIMPLIFY = FALSE
)

save(base_fit, target_base_emst, metadata_lung, lungout,
     file = "inst/extdata/lung_outfile.Rdata")

## ---- 5. Match each target to the closest fitted OMST ----------------------
emst_of <- vapply(lungout, function(f) as.numeric(f$summary_out$EMST), numeric(1))
grid_emst <- cbind(metadata_lung, EMST = emst_of)

matched <- do.call(rbind, lapply(seq_len(nrow(target_combinations)), function(i) {
  te  <- target_combinations$target_EMST[i]
  lm  <- target_combinations$LMST[i]
  sub <- grid_emst[abs(grid_emst$DMST - lm) < 1e-8, ]
  best <- sub[which.min(abs(sub$EMST - te)), ]
  data.frame(target_EMST = te, LMST = lm, matched_OMST = best$OMST,
             matched_EMST = best$EMST, abs_emst_diff = abs(best$EMST - te),
             fitID = best$fitID, risk_multiplier = risk_multiplier)
}))

## base case: exact fit, no search
matched <- rbind(matched, data.frame(
  target_EMST = target_base_emst, LMST = 1.35, matched_OMST = 4.0,
  matched_EMST = target_base_emst, abs_emst_diff = 0,
  fitID = NA_integer_, risk_multiplier = risk_multiplier))
matched <- matched[order(matched$LMST, matched$target_EMST), ]

write.csv(matched, "inst/extdata/matched_pairs.csv", row.names = FALSE)

cat("\n=== fit_lung.R done ===\n")
cat("base EMST:", round(target_base_emst, 4), "| grid fits:", nrow(metadata_lung),
    "| matched pairs:", nrow(matched), "\n")
