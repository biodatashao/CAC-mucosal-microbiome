#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_03_DMM_rarefaction_K3_multirepeat_ARI.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Stability of the K = 3 DMM solution across repeated rarefactions,
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
  "DMM_rarefaction_sensitivity_K3_multirepeat_7KB"
)

TABLE_DIR <- file.path(
  OUT_DIR,
  "tables"
)

MODEL_DIR <- file.path(
  OUT_DIR,
  "models"
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  MODEL_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Parameters
# ==============================================================================

RAREFACTION_SEED <- 20260726L

K_FIXED <- 3L

N_REPEATS <- 20L

FIT_SEED_BASE <- 20360726L

FIT_SEEDS <- FIT_SEED_BASE + seq_len(N_REPEATS) - 1L


# ==============================================================================
# 3. Check input files
# ==============================================================================

if (!file.exists(COUNT_FILE)) {
  stop(
    paste0(
      "Count file not found:\n",
      COUNT_FILE
    )
  )
}

if (!file.exists(ASSIGN_FILE)) {
  stop(
    paste0(
      "Native assignment file not found:\n",
      ASSIGN_FILE
    )
  )
}


# ==============================================================================
# 4. Read frozen primary DMM count table and native assignments
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

required_assign_cols <- c(
  "SampleID",
  "DMM_cluster"
)

missing_assign_cols <- setdiff(
  required_assign_cols,
  colnames(assign_df)
)

if (length(missing_assign_cols) > 0) {
  stop(
    paste0(
      "Native assignment table missing column(s): ",
      paste(
        missing_assign_cols,
        collapse = ", "
      )
    )
  )
}


# ==============================================================================
# 5. Construct count matrix
# ==============================================================================

sample_ids <- as.character(
  count_df$SampleID
)

feature_cols <- setdiff(
  colnames(count_df),
  "SampleID"
)

count_mat <- as.matrix(
  count_df[
    ,
    feature_cols,
    drop = FALSE
  ]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- sample_ids

if (any(!is.finite(count_mat))) {
  stop("Non-finite values detected in count matrix.")
}

if (any(count_mat < 0)) {
  stop("Negative values detected in count matrix.")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop(
    "Count matrix contains non-integer values. ",
    "This analysis requires raw genus-level read counts."
  )
}

count_mat <- round(
  count_mat
)


# ==============================================================================
# 6. Basic audit
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

cat("\n")
cat("============================================================\n")
cat("7KB DMM RAREFACTION K3 MULTI-REPEAT SENSITIVITY\n")
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

cat(
  "Rarefaction seed = ",
  RAREFACTION_SEED,
  "\n",
  sep = ""
)

cat(
  "Fixed K = ",
  K_FIXED,
  "\n",
  sep = ""
)

cat(
  "DMM repeats = ",
  N_REPEATS,
  "\n",
  sep = ""
)


# ==============================================================================
# 7. Generate ONE fixed rarefied matrix
#
# This deliberately uses the same rarefaction seed and depth as script 25b.
# The matrix is then held fixed across all 20 DMM fits.
# ==============================================================================

set.seed(
  RAREFACTION_SEED
)

rare_mat <- suppressWarnings(
  vegan::rrarefy(
    count_mat,
    sample = rare_depth
  )
)

rare_mat <- as.matrix(
  rare_mat
)

storage.mode(rare_mat) <- "numeric"

rownames(rare_mat) <- rownames(
  count_mat
)

colnames(rare_mat) <- colnames(
  count_mat
)

rare_depth_check <- rowSums(
  rare_mat
)

if (!all(rare_depth_check == rare_depth)) {
  stop("Rarefaction depth check failed.")
}

cat(
  "Uniform rarefaction depth = ",
  rare_depth,
  "\n",
  sep = ""
)


# ==============================================================================
# 8. Save fixed rarefied matrix
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
    TABLE_DIR,
    "genus_raw_count_progression127_prevalence10_rarefied_fixed_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 9. Helper: extract scalar fit statistics safely
# ==============================================================================

get_laplace <- function(fit) {
  
  value <- tryCatch(
    DirichletMultinomial::laplace(fit),
    error = function(e) NA_real_
  )
  
  value <- as.numeric(
    value
  )
  
  if (length(value) == 0) {
    return(NA_real_)
  }
  
  value[1]
}


get_aic <- function(fit) {
  
  value <- tryCatch(
    AIC(fit),
    error = function(e) NA_real_
  )
  
  value <- as.numeric(
    value
  )
  
  if (length(value) == 0) {
    return(NA_real_)
  }
  
  value[1]
}


get_bic <- function(fit) {
  
  value <- tryCatch(
    BIC(fit),
    error = function(e) NA_real_
  )
  
  value <- as.numeric(
    value
  )
  
  if (length(value) == 0) {
    return(NA_real_)
  }
  
  value[1]
}


# ==============================================================================
# 10. Run K=3 DMM 20 times
# ==============================================================================

fit_list <- vector(
  "list",
  N_REPEATS
)

fit_summary <- data.frame(
  Repeat = seq_len(N_REPEATS),
  Seed = FIT_SEEDS,
  Converged = FALSE,
  Laplace = NA_real_,
  AIC = NA_real_,
  BIC = NA_real_,
  Error_message = NA_character_,
  stringsAsFactors = FALSE
)

cat("\n")
cat("============================================================\n")
cat("RUNNING RAREFIED K3 DMM REPEATS\n")
cat("============================================================\n")

for (i in seq_len(N_REPEATS)) {
  
  seed_i <- FIT_SEEDS[i]
  
  cat(
    sprintf(
      "Repeat %02d/%02d | seed = %d ... ",
      i,
      N_REPEATS,
      seed_i
    )
  )
  
  set.seed(
    seed_i
  )
  
  fit_i <- tryCatch(
    DirichletMultinomial::dmn(
      rare_mat,
      k = K_FIXED,
      verbose = FALSE
    ),
    error = function(e) e
  )
  
  if (inherits(fit_i, "error")) {
    
    fit_summary$Error_message[i] <- conditionMessage(
      fit_i
    )
    
    cat(
      "FAILED\n"
    )
    
  } else {
    
    lap_i <- get_laplace(
      fit_i
    )
    
    aic_i <- get_aic(
      fit_i
    )
    
    bic_i <- get_bic(
      fit_i
    )
    
    fit_list[i] <- list(
      fit_i
    )
    
    fit_summary$Converged[i] <- is.finite(
      lap_i
    )
    
    fit_summary$Laplace[i] <- lap_i
    fit_summary$AIC[i] <- aic_i
    fit_summary$BIC[i] <- bic_i
    
    cat(
      sprintf(
        "Laplace = %.6f\n",
        lap_i
      )
    )
  }
}


# ==============================================================================
# 11. Save repeat-level summary
# ==============================================================================

fwrite(
  fit_summary,
  file.path(
    TABLE_DIR,
    "DMM_rarefied_K3_all_repeats_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 12. Select best model by minimum finite Laplace
# ==============================================================================

valid_idx <- which(
  fit_summary$Converged &
    is.finite(fit_summary$Laplace)
)

if (length(valid_idx) == 0) {
  stop(
    "No valid K=3 DMM fit was obtained."
  )
}

best_idx <- valid_idx[
  which.min(
    fit_summary$Laplace[
      valid_idx
    ]
  )
]

best_fit <- fit_list[
  best_idx
]

if (length(best_fit) != 1) {
  stop("Could not recover best DMM fit.")
}

best_fit <- best_fit[[1]]

best_repeat <- fit_summary$Repeat[
  best_idx
]

best_seed <- fit_summary$Seed[
  best_idx
]

best_laplace <- fit_summary$Laplace[
  best_idx
]

best_aic <- fit_summary$AIC[
  best_idx
]

best_bic <- fit_summary$BIC[
  best_idx
]


# ==============================================================================
# 13. Save best model
# ==============================================================================

saveRDS(
  best_fit,
  file.path(
    MODEL_DIR,
    "DMM_rarefied_best_model_K3_7KB.rds"
  )
)


# ==============================================================================
# 14. Extract posterior probabilities
# ==============================================================================

post <- DirichletMultinomial::mixture(
  best_fit
)

post <- as.matrix(
  post
)

if (nrow(post) != nrow(rare_mat)) {
  stop(
    paste0(
      "Posterior row count mismatch: ",
      nrow(post),
      " versus ",
      nrow(rare_mat)
    )
  )
}

if (ncol(post) != K_FIXED) {
  stop(
    paste0(
      "Expected ",
      K_FIXED,
      " posterior components, found ",
      ncol(post)
    )
  )
}

rare_component_num <- max.col(
  post,
  ties.method = "first"
)

rare_component_raw <- paste0(
  "R",
  rare_component_num
)


# ==============================================================================
# 15. Match native assignments
# ==============================================================================

native_match_idx <- match(
  rownames(rare_mat),
  assign_df$SampleID
)

if (any(is.na(native_match_idx))) {
  
  missing_samples <- rownames(rare_mat)[
    is.na(native_match_idx)
  ]
  
  stop(
    paste0(
      "Samples missing from native assignment table: ",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  )
}

native_cluster <- as.character(
  assign_df$DMM_cluster[
    native_match_idx
  ]
)

if (!all(native_cluster %in% c("C1", "C2", "C3"))) {
  stop("Unexpected native cluster labels.")
}


# ==============================================================================
# 16. Adjusted Rand Index
#
# ARI is invariant to arbitrary component labels.
# ==============================================================================

adjusted_rand_index <- function(x, y) {
  
  tab <- table(
    x,
    y
  )
  
  choose2 <- function(z) {
    z * (z - 1) / 2
  }
  
  n <- sum(
    tab
  )
  
  sum_cells <- sum(
    choose2(tab)
  )
  
  row_sum <- rowSums(
    tab
  )
  
  col_sum <- colSums(
    tab
  )
  
  sum_rows <- sum(
    choose2(row_sum)
  )
  
  sum_cols <- sum(
    choose2(col_sum)
  )
  
  total_pairs <- choose2(
    n
  )
  
  expected <- (
    sum_rows * sum_cols
  ) / total_pairs
  
  max_index <- 0.5 * (
    sum_rows + sum_cols
  )
  
  denominator <- max_index - expected
  
  if (!is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  
  (
    sum_cells - expected
  ) / denominator
}


ari <- adjusted_rand_index(
  native_cluster,
  rare_component_raw
)


# ==============================================================================
# 17. Map rarefied components to published C1/C2/C3
#
# Test all 6 permutations and maximize exact agreement with native assignments.
# This mapping does NOT affect ARI.
# ==============================================================================

permutations <- rbind(
  c(1, 2, 3),
  c(1, 3, 2),
  c(2, 1, 3),
  c(2, 3, 1),
  c(3, 1, 2),
  c(3, 2, 1)
)

mapping_eval <- data.frame(
  Permutation = character(0),
  N_exact_match = integer(0),
  Exact_agreement_rate = numeric(0),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(permutations))) {
  
  perm_i <- permutations[
    i,
    ,
    drop = TRUE
  ]
  
  mapped_i <- paste0(
    "C",
    perm_i[
      rare_component_num
    ]
  )
  
  n_match_i <- sum(
    mapped_i == native_cluster
  )
  
  mapping_eval <- rbind(
    mapping_eval,
    data.frame(
      Permutation = paste(
        perm_i,
        collapse = "-"
      ),
      N_exact_match = n_match_i,
      Exact_agreement_rate = n_match_i /
        length(native_cluster),
      stringsAsFactors = FALSE
    )
  )
}

best_mapping_idx <- which.max(
  mapping_eval$N_exact_match
)

best_perm <- permutations[
  best_mapping_idx,
  ,
  drop = TRUE
]

rare_cluster_mapped <- paste0(
  "C",
  best_perm[
    rare_component_num
  ]
)

exact_match_vector <- rare_cluster_mapped ==
  native_cluster

exact_matches <- sum(
  exact_match_vector
)

exact_agreement <- exact_matches /
  length(native_cluster)


# ==============================================================================
# 18. Contingency table
# ==============================================================================

cluster_levels <- c(
  "C1",
  "C2",
  "C3"
)

native_factor <- factor(
  native_cluster,
  levels = cluster_levels
)

rare_factor <- factor(
  rare_cluster_mapped,
  levels = cluster_levels
)

cont <- table(
  Native = native_factor,
  Rarefied = rare_factor
)


# ==============================================================================
# 19. Cluster sizes
# ==============================================================================

native_sizes <- table(
  factor(
    native_cluster,
    levels = cluster_levels
  )
)

rare_sizes <- table(
  factor(
    rare_cluster_mapped,
    levels = cluster_levels
  )
)

cluster_size_out <- data.frame(
  Cluster = cluster_levels,
  Native_N = as.integer(
    native_sizes
  ),
  Rarefied_N = as.integer(
    rare_sizes
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 20. Save sample-level assignments
# ==============================================================================

assignment_out <- data.frame(
  SampleID = rownames(rare_mat),
  Native_cluster = native_cluster,
  Rarefied_component_raw = rare_component_raw,
  Rarefied_cluster_mapped = rare_cluster_mapped,
  Exact_agreement = exact_match_vector,
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

for (j in seq_len(K_FIXED)) {
  
  new_col_name <- paste0(
    "Posterior_R",
    j
  )
  
  assignment_out[
    ,
    new_col_name
  ] <- post[
    ,
    j
  ]
}

fwrite(
  assignment_out,
  file.path(
    TABLE_DIR,
    "DMM_rarefied_bestK3_vs_native_assignments_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 21. Save disagreement samples
# ==============================================================================

disagreement_out <- assignment_out[
  !assignment_out$Exact_agreement,
  ,
  drop = FALSE
]

fwrite(
  disagreement_out,
  file.path(
    TABLE_DIR,
    "DMM_rarefied_bestK3_discordant_samples_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 22. Save contingency table
# ==============================================================================

cont_out <- data.frame(
  Native_cluster = cluster_levels,
  C1 = as.integer(
    cont[
      ,
      "C1"
    ]
  ),
  C2 = as.integer(
    cont[
      ,
      "C2"
    ]
  ),
  C3 = as.integer(
    cont[
      ,
      "C3"
    ]
  ),
  stringsAsFactors = FALSE
)

fwrite(
  cont_out,
  file.path(
    TABLE_DIR,
    "DMM_native_vs_rarefied_bestK3_contingency_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 23. Save mapping audit
# ==============================================================================

fwrite(
  mapping_eval,
  file.path(
    TABLE_DIR,
    "DMM_rarefied_bestK3_component_mapping_audit_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 24. Save cluster sizes
# ==============================================================================

fwrite(
  cluster_size_out,
  file.path(
    TABLE_DIR,
    "DMM_native_vs_rarefied_bestK3_cluster_sizes_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 25. Save best-repeat summary
# ==============================================================================

best_repeat_out <- data.frame(
  K = K_FIXED,
  Repeat = best_repeat,
  Seed = best_seed,
  Laplace = best_laplace,
  AIC = best_aic,
  BIC = best_bic,
  stringsAsFactors = FALSE
)

fwrite(
  best_repeat_out,
  file.path(
    TABLE_DIR,
    "DMM_rarefied_K3_best_repeat_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 26. Save manuscript robustness summary
# ==============================================================================

summary_out <- data.frame(
  Parameter = c(
    "Number_of_samples",
    "Number_of_genera",
    "Fixed_K",
    "Number_of_DMM_repeats",
    "Rarefaction_seed",
    "Uniform_rarefaction_depth",
    "Best_repeat",
    "Best_fit_seed",
    "Best_Laplace",
    "Best_AIC",
    "Best_BIC",
    "Adjusted_Rand_index",
    "Adjusted_Rand_index_3dp",
    "Exact_matches",
    "Exact_agreement_rate",
    "Exact_agreement_percent"
  ),
  Value = c(
    nrow(rare_mat),
    ncol(rare_mat),
    K_FIXED,
    N_REPEATS,
    RAREFACTION_SEED,
    rare_depth,
    best_repeat,
    best_seed,
    best_laplace,
    best_aic,
    best_bic,
    ari,
    round(
      ari,
      3
    ),
    exact_matches,
    exact_agreement,
    100 * exact_agreement
  ),
  stringsAsFactors = FALSE
)

fwrite(
  summary_out,
  file.path(
    TABLE_DIR,
    "DMM_rarefaction_K3_multirepeat_robustness_summary_7KB.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 27. Console report
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("BEST RAREFIED K3 MODEL\n")
cat("============================================================\n")

cat(
  sprintf(
    "Best repeat = %d\n",
    best_repeat
  )
)

cat(
  sprintf(
    "Best seed = %d\n",
    best_seed
  )
)

cat(
  sprintf(
    "Best Laplace = %.12f\n",
    best_laplace
  )
)

cat(
  sprintf(
    "Best AIC = %.12f\n",
    best_aic
  )
)

cat(
  sprintf(
    "Best BIC = %.12f\n",
    best_bic
  )
)


cat("\n")
cat("============================================================\n")
cat("FINAL RAREFACTION ROBUSTNESS RESULT\n")
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
  "Best component mapping: ",
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

for (cluster_i in cluster_levels) {
  
  vals_i <- cont[
    cluster_i,
    cluster_levels
  ]
  
  cat(
    sprintf(
      "%s: C1=%d  C2=%d  C3=%d\n",
      cluster_i,
      vals_i["C1"],
      vals_i["C2"],
      vals_i["C3"]
    )
  )
}


cat("\n")
cat("Native cluster sizes:\n")

for (cluster_i in cluster_levels) {
  
  cat(
    sprintf(
      "%s = %d\n",
      cluster_i,
      native_sizes[
        cluster_i
      ]
    )
  )
}


cat("\n")
cat("Rarefied cluster sizes:\n")

for (cluster_i in cluster_levels) {
  
  cat(
    sprintf(
      "%s = %d\n",
      cluster_i,
      rare_sizes[
        cluster_i
      ]
    )
  )
}


cat("\n")
cat("Discordant samples:\n")

if (nrow(disagreement_out) == 0) {
  
  cat("NONE\n")
  
} else {
  
  for (i in seq_len(nrow(disagreement_out))) {
    
    cat(
      sprintf(
        "%s: %s -> %s | posterior = %.6f\n",
        disagreement_out$SampleID[i],
        disagreement_out$Native_cluster[i],
        disagreement_out$Rarefied_cluster_mapped[i],
        disagreement_out$Rarefied_maximum_posterior[i]
      )
    )
  }
}


cat("\n")
cat("============================================================\n")
cat("OUTPUT DIRECTORY\n")
cat("============================================================\n")

cat(
  OUT_DIR,
  "\n"
)


# ==============================================================================
# 28. Session information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_DIR,
    "sessionInfo_DMM_rarefaction_K3_multirepeat_7KB.txt"
  )
)

cat("\nDone.\n")