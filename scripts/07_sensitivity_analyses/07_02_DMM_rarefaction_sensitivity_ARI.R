#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_02_DMM_rarefaction_sensitivity_ARI.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Stability of the DMM community-state solution across rarefaction depths,
## quantified by the adjusted Rand index.
############################################################



options(stringsAsFactors = FALSE)
options(width = 220)

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required.")
}

if (!requireNamespace("DirichletMultinomial", quietly = TRUE)) {
  stop("Package 'DirichletMultinomial' is required.")
}


# ==============================================================================
# 1. Paths
# ==============================================================================

ROOT <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "03_Figure4_DMM"
)

PRIMARY_DIR <- file.path(
  ROOT,
  "DMM_progression127_primary_prevalence10_7KB"
)

COUNT_FILE <- file.path(
  PRIMARY_DIR,
  "tables",
  "genus_raw_count_progression127_prevalence10_frozen_7KB.tsv"
)

ASSIGN_FILE <- file.path(
  PRIMARY_DIR,
  "tables",
  "DMM_sample_assignments_optimalK3_native_7KB.tsv"
)

OUT_DIR <- file.path(
  ROOT,
  "DMM_rarefaction_sensitivity_K3_7KB"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(OUT_DIR, "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Parameters
# ==============================================================================

RAREFACTION_SEED <- 20260726L
K_FIXED <- 3L


# ==============================================================================
# 3. Read native DMM input
# ==============================================================================

count_df <- fread(
  COUNT_FILE,
  data.table = FALSE,
  check.names = FALSE
)

assign_df <- fread(
  ASSIGN_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (!"SampleID" %in% colnames(count_df)) {
  stop("Count table does not contain SampleID.")
}

if (!all(
  c(
    "SampleID",
    "DMM_cluster"
  ) %in% colnames(assign_df)
)) {
  stop(
    "Assignment table must contain SampleID and DMM_cluster."
  )
}


# ==============================================================================
# 4. Count matrix
# ==============================================================================

sample_id <- as.character(
  count_df$SampleID
)

count_mat <- as.matrix(
  count_df[
    ,
    setdiff(
      colnames(count_df),
      "SampleID"
    ),
    drop = FALSE
  ]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- sample_id


if (any(count_mat < 0, na.rm = TRUE)) {
  stop("Negative counts detected.")
}

if (any(!is.finite(count_mat))) {
  stop("Non-finite counts detected.")
}


cat("\n")
cat("============================================================\n")
cat("7KB DMM RAREFACTION SENSITIVITY\n")
cat("============================================================\n")

cat(
  "Samples = ",
  nrow(count_mat),
  "\n",
  sep = ""
)

cat(
  "Genera = ",
  ncol(count_mat),
  "\n",
  sep = ""
)


# ==============================================================================
# 5. Native sequencing depths
# ==============================================================================

native_depth <- rowSums(
  count_mat
)

if (any(native_depth <= 0)) {
  stop("At least one sample has zero genus-level library size.")
}

rare_depth <- min(
  native_depth
)

cat(
  "Minimum native genus-level depth = ",
  rare_depth,
  "\n",
  sep = ""
)

cat(
  "Median native genus-level depth = ",
  median(native_depth),
  "\n",
  sep = ""
)

cat(
  "Maximum native genus-level depth = ",
  max(native_depth),
  "\n",
  sep = ""
)


# ==============================================================================
# 6. Rarefy all 127 samples to uniform depth
# ==============================================================================

set.seed(
  RAREFACTION_SEED
)

rare_mat <- vegan::rrarefy(
  count_mat,
  sample = rare_depth
)

rare_depth_check <- rowSums(
  rare_mat
)

if (!all(rare_depth_check == rare_depth)) {
  stop("Rarefaction depth check failed.")
}

cat(
  "Rarefaction depth = ",
  rare_depth,
  "\n",
  sep = ""
)

cat(
  "Seed = ",
  RAREFACTION_SEED,
  "\n",
  sep = ""
)


# ==============================================================================
# 7. Fit DMM at fixed K = 3
# ==============================================================================

cat("\nFitting rarefied DMM with K = 3...\n")

fit_rare <- DirichletMultinomial::dmn(
  rare_mat,
  k = K_FIXED,
  verbose = FALSE
)


# ==============================================================================
# 8. Extract posterior probabilities and assignments
# ==============================================================================

post <- DirichletMultinomial::mixture(
  fit_rare
)

post <- as.matrix(
  post
)

if (nrow(post) != nrow(rare_mat)) {
  stop(
    paste0(
      "Posterior row count mismatch: posterior=",
      nrow(post),
      ", samples=",
      nrow(rare_mat)
    )
  )
}

rare_cluster_num <- max.col(
  post,
  ties.method = "first"
)

rare_cluster_raw <- paste0(
  "R",
  rare_cluster_num
)


# ==============================================================================
# 9. Native assignments
# ==============================================================================

native <- assign_df[
  match(
    rownames(rare_mat),
    assign_df$SampleID
  ),
  c(
    "SampleID",
    "DMM_cluster"
  ),
  drop = FALSE
]

if (any(is.na(native$SampleID))) {
  stop("Some rarefied samples were not found in native assignment table.")
}

native_cluster <- as.character(
  native$DMM_cluster
)

if (!all(native_cluster %in% c("C1", "C2", "C3"))) {
  stop("Unexpected native DMM cluster labels.")
}


# ==============================================================================
# 10. Adjusted Rand Index
#
# Implemented directly so no additional package is required.
# ==============================================================================

adjusted_rand_index <- function(x, y) {
  
  tab <- table(
    x,
    y
  )
  
  choose2 <- function(z) {
    z * (z - 1) / 2
  }
  
  n <- sum(tab)
  
  sum_cells <- sum(
    choose2(tab)
  )
  
  row_sum <- rowSums(tab)
  col_sum <- colSums(tab)
  
  sum_rows <- sum(
    choose2(row_sum)
  )
  
  sum_cols <- sum(
    choose2(col_sum)
  )
  
  total_pairs <- choose2(n)
  
  expected <- (
    sum_rows * sum_cols
  ) / total_pairs
  
  max_index <- 0.5 * (
    sum_rows + sum_cols
  )
  
  denominator <- max_index - expected
  
  if (denominator == 0) {
    return(NA_real_)
  }
  
  (
    sum_cells - expected
  ) / denominator
}


ari <- adjusted_rand_index(
  native_cluster,
  rare_cluster_raw
)


# ==============================================================================
# 11. Map rarefied components to C1 / C2 / C3
#
# ARI itself is label-invariant.
# Mapping below is only for exact-agreement reporting.
# For K=3, evaluate all six possible permutations.
# ==============================================================================

permutations <- rbind(
  c(1, 2, 3),
  c(1, 3, 2),
  c(2, 1, 3),
  c(2, 3, 1),
  c(3, 1, 2),
  c(3, 2, 1)
)

mapping_eval <- data.frame()

for (i in seq_len(nrow(permutations))) {
  
  perm <- permutations[i, ]
  
  mapped <- paste0(
    "C",
    perm[rare_cluster_num]
  )
  
  n_match <- sum(
    mapped == native_cluster
  )
  
  tmp <- data.frame(
    Permutation = paste(
      perm,
      collapse = "-"
    ),
    N_exact_match = n_match,
    Exact_agreement = n_match / length(native_cluster),
    stringsAsFactors = FALSE
  )
  
  mapping_eval <- rbind(
    mapping_eval,
    tmp
  )
}

best_idx <- which.max(
  mapping_eval$N_exact_match
)

best_perm <- permutations[
  best_idx,
  ,
  drop = TRUE
]

rare_cluster_mapped <- paste0(
  "C",
  best_perm[rare_cluster_num]
)

exact_matches <- sum(
  rare_cluster_mapped == native_cluster
)

exact_agreement <- exact_matches /
  length(native_cluster)


# ==============================================================================
# 12. Contingency table
# ==============================================================================

cont <- table(
  Native = native_cluster,
  Rarefied = rare_cluster_mapped
)


# ==============================================================================
# 13. Console summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("RAREFACTION ROBUSTNESS RESULTS\n")
cat("============================================================\n")

cat(
  sprintf(
    "Adjusted Rand index = %.12f\n",
    ari
  )
)

cat(
  sprintf(
    "Adjusted Rand index, 3 decimals = %.3f\n",
    ari
  )
)

cat(
  "Best rarefied component mapping: ",
  paste(
    paste0(
      "R",
      1:3,
      "->C",
      best_perm
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Exact agreement = %d/%d = %.1f%%\n",
    exact_matches,
    length(native_cluster),
    100 * exact_agreement
  )
)

cat("\nNative vs rarefied K3 contingency:\n")

for (native_i in c("C1", "C2", "C3")) {
  
  vals <- cont[
    native_i,
    c("C1", "C2", "C3")
  ]
  
  cat(
    sprintf(
      "%s: C1=%d  C2=%d  C3=%d\n",
      native_i,
      vals["C1"],
      vals["C2"],
      vals["C3"]
    )
  )
}


# ==============================================================================
# 14. Save assignment table
# ==============================================================================

assignment_out <- data.frame(
  SampleID = rownames(rare_mat),
  Native_cluster = native_cluster,
  Rarefied_component = rare_cluster_raw,
  Rarefied_cluster_mapped = rare_cluster_mapped,
  Exact_agreement = rare_cluster_mapped == native_cluster,
  Rarefied_maximum_posterior = apply(
    post,
    1,
    max
  ),
  Native_genus_level_depth = native_depth[
    rownames(rare_mat)
  ],
  Rarefied_depth = rare_depth_check,
  stringsAsFactors = FALSE
)

fwrite(
  assignment_out,
  file.path(
    OUT_DIR,
    "tables",
    "DMM_K3_rarefied_vs_native_assignments_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 15. Save contingency
# ==============================================================================

cont_df <- as.data.frame.matrix(
  cont
)

cont_df$Native_cluster <- rownames(
  cont_df
)

cont_df <- cont_df[
  ,
  c(
    "Native_cluster",
    "C1",
    "C2",
    "C3"
  )
]

rownames(cont_df) <- NULL

fwrite(
  cont_df,
  file.path(
    OUT_DIR,
    "tables",
    "DMM_K3_native_vs_rarefied_contingency_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 16. Save mapping audit
# ==============================================================================

fwrite(
  mapping_eval,
  file.path(
    OUT_DIR,
    "tables",
    "DMM_K3_rarefied_component_mapping_audit_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 17. Save summary
# ==============================================================================

summary_out <- data.frame(
  Parameter = c(
    "Number_of_samples",
    "Number_of_genera",
    "Fixed_K",
    "Rarefaction_seed",
    "Uniform_rarefaction_depth",
    "Adjusted_Rand_index",
    "Exact_matches",
    "Exact_agreement_rate"
  ),
  Value = c(
    nrow(rare_mat),
    ncol(rare_mat),
    K_FIXED,
    RAREFACTION_SEED,
    rare_depth,
    ari,
    exact_matches,
    exact_agreement
  ),
  stringsAsFactors = FALSE
)

fwrite(
  summary_out,
  file.path(
    OUT_DIR,
    "tables",
    "DMM_K3_rarefaction_robustness_summary_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 18. Save rarefied matrix
# ==============================================================================

rare_out <- data.frame(
  SampleID = rownames(rare_mat),
  rare_mat,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

fwrite(
  rare_out,
  file.path(
    OUT_DIR,
    "tables",
    "genus_raw_count_progression127_prevalence10_rarefied_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 19. Session info
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_DIR,
    "sessionInfo_DMM_rarefaction_sensitivity_7KB.txt"
  )
)


cat("\n")
cat("============================================================\n")
cat("OUTPUT DIRECTORY\n")
cat("============================================================\n")

cat(
  OUT_DIR,
  "\n"
)

cat("\nDone.\n")