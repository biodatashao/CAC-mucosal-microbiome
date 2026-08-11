#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_06_Figure7H_RMI_Cliffs_delta_CI.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Cliff's delta with confidence interval for the Recurrence-associated
## Microbial Index.
############################################################



options(stringsAsFactors = FALSE)
options(width = 220)

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

SOURCE_FILE <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "04_Figure7_clinical_recurrence/",
  "02_DHI_original_method/",
  "Figure7HI_RMI_LOOCV/",
  "Figure7HI_original_LOOCV_RMI_source_7KB.tsv"
)

OUT_DIR <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "05_manuscript_audit"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(SOURCE_FILE)) {
  stop(
    paste0(
      "Source file not found:\n",
      SOURCE_FILE
    )
  )
}

# ==============================================================================
# 2. Parameters
# ==============================================================================

BOOT_SEED <- 20260726L
N_BOOT <- 5000L

# ==============================================================================
# 3. Read formal Figure 7H source
# ==============================================================================

dat <- fread(
  SOURCE_FILE,
  data.table = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "SampleID",
  "Recurrence_plot",
  "RMI_two_genus"
)

missing_cols <- setdiff(
  required_cols,
  colnames(dat)
)

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}

# ==============================================================================
# 4. Standardize recurrence status
# ==============================================================================

status_raw <- tolower(
  trimws(
    as.character(
      dat$Recurrence_plot
    )
  )
)

status <- rep(
  NA_character_,
  length(status_raw)
)

status[
  status_raw %in% c(
    "recurrence",
    "recurrent",
    "yes",
    "y",
    "1",
    "true"
  )
] <- "Recurrence"

status[
  status_raw %in% c(
    "no recurrence",
    "no_recurrence",
    "non-recurrent",
    "nonrecurrent",
    "no",
    "n",
    "0",
    "false"
  )
] <- "No recurrence"

rmi <- suppressWarnings(
  as.numeric(
    dat$RMI_two_genus
  )
)

audit <- data.frame(
  SampleID = dat$SampleID,
  Recurrence = status,
  RMI_two_genus = rmi,
  stringsAsFactors = FALSE
)

audit <- audit[
  !is.na(audit$Recurrence) &
    is.finite(audit$RMI_two_genus),
  ,
  drop = FALSE
]

# ==============================================================================
# 5. Split groups
# ==============================================================================

x <- audit$RMI_two_genus[
  audit$Recurrence == "Recurrence"
]

y <- audit$RMI_two_genus[
  audit$Recurrence == "No recurrence"
]

if (length(x) != 10 || length(y) != 13) {
  warning(
    paste0(
      "Expected recurrence n=10 and no recurrence n=13, found ",
      length(x),
      " and ",
      length(y),
      "."
    )
  )
}

# ==============================================================================
# 6. Wilcoxon test
# ==============================================================================

wilcox_res <- wilcox.test(
  x,
  y,
  alternative = "two.sided",
  exact = FALSE
)

# ==============================================================================
# 7. Cliff's delta
# ==============================================================================

calc_cliff <- function(a, b) {
  
  cmp <- outer(
    a,
    b,
    FUN = "-"
  )
  
  (
    sum(cmp > 0) -
      sum(cmp < 0)
  ) / (
    length(a) * length(b)
  )
}

delta <- calc_cliff(
  x,
  y
)

# ==============================================================================
# 8. Stratified bootstrap 95% CI
#
# Resample recurrence and no-recurrence groups independently,
# preserving group sizes.
# ==============================================================================

set.seed(
  BOOT_SEED
)

boot_delta <- numeric(
  N_BOOT
)

for (i in seq_len(N_BOOT)) {
  
  xb <- sample(
    x,
    size = length(x),
    replace = TRUE
  )
  
  yb <- sample(
    y,
    size = length(y),
    replace = TRUE
  )
  
  boot_delta[i] <- calc_cliff(
    xb,
    yb
  )
}

ci <- quantile(
  boot_delta,
  probs = c(
    0.025,
    0.975
  ),
  na.rm = TRUE,
  names = FALSE,
  type = 7
)

# ==============================================================================
# 9. Descriptive statistics
# ==============================================================================

desc <- data.frame(
  Group = c(
    "Recurrence",
    "No recurrence"
  ),
  N = c(
    length(x),
    length(y)
  ),
  Median = c(
    median(x),
    median(y)
  ),
  Q1 = c(
    unname(quantile(x, 0.25)),
    unname(quantile(y, 0.25))
  ),
  Q3 = c(
    unname(quantile(x, 0.75)),
    unname(quantile(y, 0.75))
  ),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 10. Pairwise comparison audit
# ==============================================================================

cmp <- outer(
  x,
  y,
  FUN = "-"
)

n_greater <- sum(
  cmp > 0
)

n_less <- sum(
  cmp < 0
)

n_ties <- sum(
  cmp == 0
)

n_pairs <- length(x) * length(y)

# ==============================================================================
# 11. Console output
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 7H FORMAL ORIGINAL-METHOD RMI AUDIT\n")
cat("============================================================\n")

cat(
  sprintf(
    "Recurrence N = %d\n",
    length(x)
  )
)

cat(
  sprintf(
    "No recurrence N = %d\n",
    length(y)
  )
)

cat(
  sprintf(
    "Wilcoxon rank-sum P = %.15f\n",
    wilcox_res$p.value
  )
)

cat(
  sprintf(
    "Cliff's delta = %.15f\n",
    delta
  )
)

cat(
  sprintf(
    "Cliff's delta, 3 decimals = %.3f\n",
    delta
  )
)

cat(
  sprintf(
    "Bootstrap replicates = %d\n",
    N_BOOT
  )
)

cat(
  sprintf(
    "Bootstrap seed = %d\n",
    BOOT_SEED
  )
)

cat(
  sprintf(
    "Bootstrap 95%% CI = %.15f to %.15f\n",
    ci[1],
    ci[2]
  )
)

cat(
  sprintf(
    "Bootstrap 95%% CI, 3 decimals = %.3f to %.3f\n",
    ci[1],
    ci[2]
  )
)

cat(
  sprintf(
    "Pairwise comparisons: greater=%d, less=%d, ties=%d, total=%d\n",
    n_greater,
    n_less,
    n_ties,
    n_pairs
  )
)

cat("\nDescriptive statistics:\n")

write.table(
  desc,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

# ==============================================================================
# 12. Save outputs
# ==============================================================================

summary_out <- data.frame(
  Metric = c(
    "Recurrence_N",
    "No_recurrence_N",
    "Wilcoxon_rank_sum_P",
    "Cliffs_delta",
    "Cliffs_delta_CI_lower_95",
    "Cliffs_delta_CI_upper_95",
    "Bootstrap_replicates",
    "Bootstrap_seed"
  ),
  Value = c(
    length(x),
    length(y),
    wilcox_res$p.value,
    delta,
    ci[1],
    ci[2],
    N_BOOT,
    BOOT_SEED
  ),
  stringsAsFactors = FALSE
)

fwrite(
  summary_out,
  file.path(
    OUT_DIR,
    "Figure7H_original_RMI_Cliffs_delta_CI_audit_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  desc,
  file.path(
    OUT_DIR,
    "Figure7H_original_RMI_descriptive_statistics_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  audit,
  file.path(
    OUT_DIR,
    "Figure7H_original_RMI_sample_values_audit_7KB.tsv"
  ),
  sep = "\t"
)

cat("\n")
cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")
cat(
  OUT_DIR,
  "\n"
)