############################################################
## config.R
##
## Central path configuration for all analysis scripts.
##
## Every script in scripts/ sources this file and then refers
## to PROJECT_ROOT instead of any absolute path.
##
## To reproduce the analysis on your own machine:
##   Option A (recommended) - set an environment variable:
##       export FFPE_PROJECT_ROOT=/path/to/your/data
##   Option B - edit the default below.
############################################################

PROJECT_ROOT <- Sys.getenv(
  "FFPE_PROJECT_ROOT",
  unset = normalizePath(".", mustWork = FALSE)
)

if (!dir.exists(PROJECT_ROOT)) {
  stop(
    "PROJECT_ROOT does not exist: ", PROJECT_ROOT,
    "\nSet the FFPE_PROJECT_ROOT environment variable or edit config.R."
  )
}

## Expected layout under PROJECT_ROOT:
##   PROJECT_ROOT/
##     raw/                                     raw sequences and taxonomy
##     data/clinical/                           clinical metadata (not distributed)
##     output/decontamination_8controls/        first decontamination run; read only
##                                              to reconstruct the cohort definition
##     output/decontamination_7controls_final/  decontamination used throughout the
##                                              manuscript (seven negative controls)
##     output/analysis/                         all analysis outputs (modules 01-08)

OUTPUT_ROOT <- file.path(PROJECT_ROOT, "output", "analysis")

## Global random seed. Scripts that perform permutation, bootstrap,
## rarefaction or DMM fitting set their own seed; this is the fallback.
GLOBAL_SEED <- 20260726
set.seed(GLOBAL_SEED)
