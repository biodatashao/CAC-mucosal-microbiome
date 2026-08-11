#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_05_Figure7H_Cliffs_delta.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Cliff's delta effect size for the recurrence comparison shown in Figure 7H.
## Computed directly from all pairwise comparisons between groups.
############################################################



options(stringsAsFactors = FALSE)
options(width = 220)

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- file.path(PROJECT_ROOT, "output/analysis/04_Figure7_clinical_recurrence")

if (!dir.exists(ROOT)) {
  stop(
    paste0(
      "Figure 7 output directory not found:\n",
      ROOT
    )
  )
}


# ==============================================================================
# 1. Locate candidate Figure 7H source tables
# ==============================================================================

all_files <- list.files(
  ROOT,
  recursive = TRUE,
  full.names = TRUE
)

candidate_files <- all_files[
  grepl(
    "index|Figure7H|two.genus|two_genus|recurrence",
    basename(all_files),
    ignore.case = TRUE
  ) &
    grepl(
      "\\.(csv|tsv|txt)$",
      all_files,
      ignore.case = TRUE
    )
]

cat("\n")
cat("============================================================\n")
cat("CANDIDATE FIGURE 7H FILES\n")
cat("============================================================\n")

for (f in candidate_files) {
  cat(f, "\n")
}


# ==============================================================================
# 2. Find a table containing recurrence status and microbial index
# ==============================================================================

selected_file <- NA_character_
selected_dat <- NULL
selected_status_col <- NA_character_
selected_index_col <- NA_character_

status_patterns <- c(
  "recurrence",
  "recur",
  "recurrence_status",
  "Recurrence"
)

index_patterns <- c(
  "index",
  "microbial_index",
  "two_genus",
  "score",
  "combined"
)

for (f in candidate_files) {
  
  dat <- tryCatch(
    fread(
      f,
      data.table = FALSE,
      check.names = FALSE,
      showProgress = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat) || nrow(dat) == 0) {
    next
  }
  
  cn <- colnames(dat)
  
  status_hits <- cn[
    grepl(
      paste(status_patterns, collapse = "|"),
      cn,
      ignore.case = TRUE
    )
  ]
  
  index_hits <- cn[
    grepl(
      paste(index_patterns, collapse = "|"),
      cn,
      ignore.case = TRUE
    )
  ]
  
  numeric_index_hits <- index_hits[
    vapply(
      dat[index_hits],
      is.numeric,
      logical(1)
    )
  ]
  
  if (
    length(status_hits) >= 1 &&
    length(numeric_index_hits) >= 1
  ) {
    
    selected_file <- f
    selected_dat <- dat
    selected_status_col <- status_hits[1]
    selected_index_col <- numeric_index_hits[1]
    
    break
  }
}


# ==============================================================================
# 3. Stop if automatic detection fails
# ==============================================================================

if (is.null(selected_dat)) {
  
  cat("\n")
  cat("============================================================\n")
  cat("NO AUTOMATIC SOURCE TABLE FOUND\n")
  cat("============================================================\n")
  
  cat(
    paste0(
      "Could not automatically identify a table containing both ",
      "recurrence status and the Figure 7H microbial index.\n"
    )
  )
  
  cat(
    "Please paste the CANDIDATE FIGURE 7H FILES section back to me.\n"
  )
  
  quit(
    save = "no",
    status = 0
  )
}


# ==============================================================================
# 4. Report detected source
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("DETECTED FIGURE 7H SOURCE\n")
cat("============================================================\n")

cat(
  "File: ",
  selected_file,
  "\n",
  sep = ""
)

cat(
  "Recurrence column: ",
  selected_status_col,
  "\n",
  sep = ""
)

cat(
  "Index column: ",
  selected_index_col,
  "\n",
  sep = ""
)


# ==============================================================================
# 5. Standardize recurrence status
# ==============================================================================

status_raw <- trimws(
  as.character(
    selected_dat[, selected_status_col]
  )
)

status_key <- tolower(
  status_raw
)

recurrence <- rep(
  NA_character_,
  length(status_key)
)

recurrence[
  status_key %in% c(
    "yes",
    "y",
    "1",
    "true",
    "recurrent",
    "recurrence",
    "rec"
  )
] <- "Recurrence"

recurrence[
  status_key %in% c(
    "no",
    "n",
    "0",
    "false",
    "non-recurrent",
    "nonrecurrent",
    "no recurrence",
    "without recurrence",
    "nonrec"
  )
] <- "No recurrence"


# ==============================================================================
# 6. Extract index values
# ==============================================================================

index_value <- suppressWarnings(
  as.numeric(
    selected_dat[, selected_index_col]
  )
)

audit <- data.frame(
  Recurrence_raw = status_raw,
  Recurrence = recurrence,
  Index = index_value,
  stringsAsFactors = FALSE
)

audit <- audit[
  !is.na(audit$Recurrence) &
    is.finite(audit$Index),
  ,
  drop = FALSE
]


# ==============================================================================
# 7. Check expected CAC sample counts
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("GROUP COUNTS\n")
cat("============================================================\n")

group_counts <- table(
  audit$Recurrence
)

print(
  group_counts
)

if (nrow(audit) != 23) {
  warning(
    paste0(
      "Expected 23 evaluable CAC samples, found ",
      nrow(audit),
      "."
    )
  )
}


# ==============================================================================
# 8. Define direction
#
# Figure 7H index is expected to be higher in recurrent CAC:
# higher UCG-005 + lower Lactococcus.
#
# Cliff's delta > 0 therefore means:
# Recurrence group has higher index than No recurrence group.
# ==============================================================================

x <- audit$Index[
  audit$Recurrence == "Recurrence"
]

y <- audit$Index[
  audit$Recurrence == "No recurrence"
]

if (length(x) == 0 || length(y) == 0) {
  stop("One recurrence group has zero observations.")
}


# ==============================================================================
# 9. Wilcoxon rank-sum test
# ==============================================================================

wilcox_res <- wilcox.test(
  x,
  y,
  alternative = "two.sided",
  exact = FALSE
)


# ==============================================================================
# 10. Cliff's delta from all pairwise comparisons
#
# delta = P(X > Y) - P(X < Y)
#
# Ties contribute 0.
# ==============================================================================

pair_cmp <- outer(
  x,
  y,
  FUN = "-"
)

n_greater <- sum(
  pair_cmp > 0
)

n_less <- sum(
  pair_cmp < 0
)

n_ties <- sum(
  pair_cmp == 0
)

n_pairs <- length(x) * length(y)

cliff_delta <- (
  n_greater - n_less
) / n_pairs


# ==============================================================================
# 11. Bootstrap 95% CI
#
# Stratified bootstrap:
# resample recurrent and non-recurrent groups independently,
# preserving group sizes.
#
# Fixed seed for reproducibility.
# ==============================================================================

BOOT_SEED <- 20260726L
N_BOOT <- 10000L

set.seed(
  BOOT_SEED
)

boot_delta <- numeric(
  N_BOOT
)

for (b in seq_len(N_BOOT)) {
  
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
  
  cmp_b <- outer(
    xb,
    yb,
    FUN = "-"
  )
  
  boot_delta[b] <- (
    sum(cmp_b > 0) -
      sum(cmp_b < 0)
  ) / (
    length(xb) * length(yb)
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
# 12. Descriptive statistics
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
    unname(
      quantile(
        x,
        0.25
      )
    ),
    unname(
      quantile(
        y,
        0.25
      )
    )
  ),
  Q3 = c(
    unname(
      quantile(
        x,
        0.75
      )
    ),
    unname(
      quantile(
        y,
        0.75
      )
    )
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 13. Console output
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("FIGURE 7H AUDIT RESULTS\n")
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
    "Wilcoxon rank-sum P = %.12f\n",
    wilcox_res$p.value
  )
)

cat(
  sprintf(
    "Cliff's delta = %.12f\n",
    cliff_delta
  )
)

cat(
  sprintf(
    "Cliff's delta, 3 decimals = %.3f\n",
    cliff_delta
  )
)

cat(
  sprintf(
    "Bootstrap 95%% CI = %.12f to %.12f\n",
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
# 14. Save audit outputs
# ==============================================================================

OUT_DIR <- file.path(
  file.path(PROJECT_ROOT, "output/analysis"),
  "05_manuscript_audit"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

summary_out <- data.frame(
  Metric = c(
    "Recurrence_N",
    "No_recurrence_N",
    "Wilcoxon_rank_sum_P",
    "Cliffs_delta",
    "Cliffs_delta_CI_lower_95",
    "Cliffs_delta_CI_upper_95",
    "Bootstrap_iterations",
    "Bootstrap_seed"
  ),
  Value = c(
    length(x),
    length(y),
    wilcox_res$p.value,
    cliff_delta,
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
    "Figure7H_two_genus_index_Cliffs_delta_audit_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  desc,
  file.path(
    OUT_DIR,
    "Figure7H_two_genus_index_descriptive_statistics_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  audit,
  file.path(
    OUT_DIR,
    "Figure7H_two_genus_index_sample_values_audit_7KB.tsv"
  ),
  sep = "\t"
)

cat("\n")
cat("============================================================\n")
cat("DONE\n")
cat("============================================================\n")

cat(
  "Audit outputs saved to:\n",
  OUT_DIR,
  "\n",
  sep = ""
)