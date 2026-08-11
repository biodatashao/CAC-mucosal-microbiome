#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_04_sequencing_depth_rarefaction_sensitivity.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Sequencing depth summary and alpha-diversity sensitivity after rarefaction
## to the minimum library size. Source for the rarefaction supplementary table.
############################################################



options(stringsAsFactors = FALSE)
options(width = 240)

suppressPackageStartupMessages({
  library(data.table)
})

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop("Package 'vegan' is required.")
}

# ==============================================================================
# 1. Paths
# ==============================================================================

CLEAN_DIR <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "00_clean_data/progression127"
)

COUNT_FILE <- file.path(
  CLEAN_DIR,
  "asv_count_7KB_progression127.tsv"
)

META_FILE <- file.path(
  CLEAN_DIR,
  "metadata_7KB_progression127.tsv"
)

ALPHA_SCRIPT <- paste0(
  file.path(PROJECT_ROOT, "script/analysis/"),
  "02_alpha_beta_statistics_7KB.R"
)

OUT_DIR <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "05_manuscript_audit/sequencing_depth_alpha_7KB"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

GROUP_LEVELS <- c(
  "Polyp",
  "UC remission",
  "UC active",
  "Dysplasia",
  "CAC"
)

RAREFACTION_SEED <- 20260726L


# ==============================================================================
# 2. Input checks
# ==============================================================================

if (!file.exists(COUNT_FILE)) {
  stop(
    paste0(
      "Count file not found:\n",
      COUNT_FILE
    )
  )
}

if (!file.exists(META_FILE)) {
  stop(
    paste0(
      "Metadata file not found:\n",
      META_FILE
    )
  )
}


# ==============================================================================
# 3. Read progression127 data
# ==============================================================================

count_df <- fread(
  COUNT_FILE,
  data.table = FALSE,
  check.names = FALSE
)

meta <- fread(
  META_FILE,
  data.table = FALSE,
  check.names = FALSE
)


# ==============================================================================
# 4. Detect count-table orientation
# ==============================================================================

sample_id_candidates <- c(
  "SampleID",
  "Sample_ID",
  "sample_id",
  "sample"
)

sample_col_count <- intersect(
  sample_id_candidates,
  colnames(count_df)
)

sample_col_meta <- intersect(
  sample_id_candidates,
  colnames(meta)
)

if (length(sample_col_meta) == 0) {
  stop(
    paste0(
      "Could not identify SampleID column in metadata.\nColumns: ",
      paste(
        colnames(meta),
        collapse = " | "
      )
    )
  )
}

sample_col_meta <- sample_col_meta[1]

meta_sample_ids <- as.character(
  meta[, sample_col_meta]
)

# Case A:
# samples are rows, one SampleID column + ASV columns
if (length(sample_col_count) > 0) {
  
  sample_col_count <- sample_col_count[1]
  
  count_sample_ids <- as.character(
    count_df[, sample_col_count]
  )
  
  feature_cols <- setdiff(
    colnames(count_df),
    sample_col_count
  )
  
  count_mat <- as.matrix(
    count_df[
      ,
      feature_cols,
      drop = FALSE
    ]
  )
  
  storage.mode(count_mat) <- "numeric"
  
  rownames(count_mat) <- count_sample_ids
  
  # Case B:
  # ASVs are rows and samples are columns
} else {
  
  first_col <- colnames(count_df)[1]
  
  overlap_columns <- sum(
    colnames(count_df) %in% meta_sample_ids
  )
  
  if (overlap_columns < 100) {
    stop(
      paste0(
        "Could not determine count-table orientation.\n",
        "Sample-column overlap with metadata = ",
        overlap_columns
      )
    )
  }
  
  asv_ids <- as.character(
    count_df[, first_col]
  )
  
  sample_cols <- intersect(
    meta_sample_ids,
    colnames(count_df)
  )
  
  tmp <- as.matrix(
    count_df[
      ,
      sample_cols,
      drop = FALSE
    ]
  )
  
  storage.mode(tmp) <- "numeric"
  
  rownames(tmp) <- asv_ids
  
  count_mat <- t(
    tmp
  )
}


# ==============================================================================
# 5. Match metadata and counts
# ==============================================================================

common_samples <- intersect(
  meta_sample_ids,
  rownames(count_mat)
)

if (length(common_samples) != 127) {
  stop(
    paste0(
      "Expected 127 progression samples; matched ",
      length(common_samples),
      "."
    )
  )
}

count_mat <- count_mat[
  common_samples,
  ,
  drop = FALSE
]

meta_idx <- match(
  common_samples,
  meta_sample_ids
)

meta <- meta[
  meta_idx,
  ,
  drop = FALSE
]

if (!all(
  as.character(meta[, sample_col_meta]) ==
  rownames(count_mat)
)) {
  stop("Metadata/count sample order check failed.")
}

if (any(!is.finite(count_mat))) {
  stop("Non-finite count values detected.")
}

if (any(count_mat < 0)) {
  stop("Negative count values detected.")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop("Non-integer values detected in ASV count table.")
}

count_mat <- round(
  count_mat
)


# ==============================================================================
# 6. Detect progression group
# ==============================================================================

candidate_group_cols <- c(
  "Progression5",
  "progression5",
  "Group",
  "group",
  "Group_main",
  "Group_raw"
)

group_col <- intersect(
  candidate_group_cols,
  colnames(meta)
)

if (length(group_col) == 0) {
  stop(
    paste0(
      "Could not identify progression group column.\nColumns: ",
      paste(
        colnames(meta),
        collapse = " | "
      )
    )
  )
}

group_col <- group_col[1]

group_raw <- trimws(
  as.character(
    meta[, group_col]
  )
)

standardize_group <- function(x) {
  
  key <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  out <- rep(
    NA_character_,
    length(key)
  )
  
  out[
    key %in% c(
      "polyp"
    )
  ] <- "Polyp"
  
  out[
    key %in% c(
      "uc_remission",
      "uc remission",
      "remission",
      "mayo0"
    )
  ] <- "UC remission"
  
  out[
    key %in% c(
      "uc_active",
      "uc active",
      "active",
      "mayo1"
    )
  ] <- "UC active"
  
  out[
    key %in% c(
      "dysplasia"
    )
  ] <- "Dysplasia"
  
  out[
    key %in% c(
      "cac",
      "ca",
      "cancer"
    )
  ] <- "CAC"
  
  out
}

group <- standardize_group(
  group_raw
)

if (any(is.na(group))) {
  
  bad <- unique(
    group_raw[
      is.na(group)
    ]
  )
  
  stop(
    paste0(
      "Unrecognized progression group value(s): ",
      paste(
        bad,
        collapse = ", "
      )
    )
  )
}

group <- factor(
  group,
  levels = GROUP_LEVELS
)


# ==============================================================================
# 7. Confirm sample counts
# ==============================================================================

group_counts <- table(
  group
)

expected_counts <- c(
  "Polyp" = 26,
  "UC remission" = 36,
  "UC active" = 25,
  "Dysplasia" = 17,
  "CAC" = 23
)

if (!all(
  as.integer(group_counts[names(expected_counts)]) ==
  as.integer(expected_counts)
)) {
  
  print(group_counts)
  
  stop(
    "Progression group counts do not match expected 26/36/25/17/23."
  )
}


# ==============================================================================
# 8. Sequencing depth per sample
# ==============================================================================

depth <- rowSums(
  count_mat
)

if (any(depth <= 0)) {
  stop("At least one progression sample has zero sequencing depth.")
}

depth_sample <- data.frame(
  SampleID = rownames(count_mat),
  Group = as.character(group),
  Sequencing_depth = depth,
  stringsAsFactors = FALSE
)


# ==============================================================================
# 9. Group sequencing-depth descriptive statistics
# ==============================================================================

depth_summary_list <- list()

for (g in GROUP_LEVELS) {
  
  x <- depth[
    group == g
  ]
  
  tmp <- data.frame(
    Group = g,
    N = length(x),
    Min = min(x),
    Q1 = unname(
      quantile(
        x,
        0.25
      )
    ),
    Median = median(x),
    Q3 = unname(
      quantile(
        x,
        0.75
      )
    ),
    Max = max(x),
    IQR = IQR(x),
    Mean = mean(x),
    SD = sd(x),
    stringsAsFactors = FALSE
  )
  
  depth_summary_list[length(depth_summary_list) + 1] <- list(
    tmp
  )
}

depth_summary <- do.call(
  rbind,
  depth_summary_list
)


# ==============================================================================
# 10. Test sequencing-depth differences across five groups
# ==============================================================================

depth_kw <- kruskal.test(
  depth ~ group
)

depth_test <- data.frame(
  Test = "Kruskal-Wallis",
  N = length(depth),
  Groups = 5,
  Statistic = unname(
    depth_kw$statistic
  ),
  df = unname(
    depth_kw$parameter
  ),
  P_value = depth_kw$p.value,
  stringsAsFactors = FALSE
)


# ==============================================================================
# 11. CAC versus pooled other progression groups
#
# This is only a technical depth check, not the biological primary comparison.
# ==============================================================================

cac_status <- ifelse(
  group == "CAC",
  "CAC",
  "Other four groups"
)

depth_cac_wilcox <- wilcox.test(
  depth ~ cac_status,
  exact = FALSE
)

depth_cac_test <- data.frame(
  Comparison = "CAC vs other four progression groups",
  CAC_N = sum(cac_status == "CAC"),
  Other_N = sum(cac_status == "Other four groups"),
  CAC_median_depth = median(
    depth[
      cac_status == "CAC"
    ]
  ),
  Other_median_depth = median(
    depth[
      cac_status == "Other four groups"
    ]
  ),
  Wilcoxon_P = depth_cac_wilcox$p.value,
  stringsAsFactors = FALSE
)


# ==============================================================================
# 12. Calculate primary alpha directly on non-rarefied counts
#
# These are calculated here independently as an audit.
# ==============================================================================

alpha_native <- data.frame(
  SampleID = rownames(count_mat),
  Group = as.character(group),
  Sequencing_depth = depth,
  Observed_ASVs = rowSums(
    count_mat > 0
  ),
  Shannon = vegan::diversity(
    count_mat,
    index = "shannon"
  ),
  Simpson = vegan::diversity(
    count_mat,
    index = "simpson"
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 13. Native alpha five-group Kruskal-Wallis
# ==============================================================================

native_alpha_tests <- data.frame()

for (metric in c(
  "Observed_ASVs",
  "Shannon",
  "Simpson"
)) {
  
  test_i <- kruskal.test(
    alpha_native[, metric] ~ group
  )
  
  tmp <- data.frame(
    Dataset = "Native_non_rarefied",
    Metric = metric,
    Statistic = unname(
      test_i$statistic
    ),
    df = unname(
      test_i$parameter
    ),
    P_value = test_i$p.value,
    stringsAsFactors = FALSE
  )
  
  native_alpha_tests <- rbind(
    native_alpha_tests,
    tmp
  )
}


# ==============================================================================
# 14. Inspect actual alpha/beta analysis script
#
# This determines whether the existing submitted analysis explicitly rarefied
# counts before alpha-diversity calculation.
# ==============================================================================

alpha_script_audit <- data.frame(
  Script_exists = file.exists(ALPHA_SCRIPT),
  Script = ALPHA_SCRIPT,
  Contains_rrarefy = FALSE,
  Contains_rarefy_even_depth = FALSE,
  Contains_rarefaction_word = FALSE,
  Contains_vegan_diversity = FALSE,
  Contains_estimate_richness = FALSE,
  stringsAsFactors = FALSE
)

rarefaction_matching_lines <- character(0)
alpha_matching_lines <- character(0)

if (file.exists(ALPHA_SCRIPT)) {
  
  script_lines <- readLines(
    ALPHA_SCRIPT,
    warn = FALSE
  )
  
  alpha_script_audit$Contains_rrarefy <- any(
    grepl(
      "rrarefy",
      script_lines,
      ignore.case = TRUE
    )
  )
  
  alpha_script_audit$Contains_rarefy_even_depth <- any(
    grepl(
      "rarefy_even_depth",
      script_lines,
      ignore.case = TRUE
    )
  )
  
  alpha_script_audit$Contains_rarefaction_word <- any(
    grepl(
      "raref|subsampl",
      script_lines,
      ignore.case = TRUE
    )
  )
  
  alpha_script_audit$Contains_vegan_diversity <- any(
    grepl(
      "diversity\\s*\\(",
      script_lines,
      ignore.case = TRUE
    )
  )
  
  alpha_script_audit$Contains_estimate_richness <- any(
    grepl(
      "estimate_richness",
      script_lines,
      ignore.case = TRUE
    )
  )
  
  rare_idx <- grep(
    "raref|rrarefy|subsampl",
    script_lines,
    ignore.case = TRUE
  )
  
  if (length(rare_idx) > 0) {
    
    rarefaction_matching_lines <- paste0(
      rare_idx,
      ": ",
      script_lines[
        rare_idx
      ]
    )
  }
  
  alpha_idx <- grep(
    "Shannon|Simpson|Observed|diversity\\s*\\(|estimate_richness",
    script_lines,
    ignore.case = TRUE
  )
  
  if (length(alpha_idx) > 0) {
    
    alpha_matching_lines <- paste0(
      alpha_idx,
      ": ",
      script_lines[
        alpha_idx
      ]
    )
  }
}


# ==============================================================================
# 15. Rarefaction sensitivity analysis
#
# Rarefy the SAME progression127 ASV count matrix to the minimum observed
# library size across the 127 samples.
# ==============================================================================

rarefaction_depth <- min(
  depth
)

set.seed(
  RAREFACTION_SEED
)

rare_mat <- suppressWarnings(
  vegan::rrarefy(
    count_mat,
    sample = rarefaction_depth
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

if (!all(
  rowSums(rare_mat) == rarefaction_depth
)) {
  stop("Rarefaction depth validation failed.")
}


# ==============================================================================
# 16. Rarefied alpha
# ==============================================================================

alpha_rare <- data.frame(
  SampleID = rownames(rare_mat),
  Group = as.character(group),
  Rarefaction_depth = rowSums(
    rare_mat
  ),
  Observed_ASVs = rowSums(
    rare_mat > 0
  ),
  Shannon = vegan::diversity(
    rare_mat,
    index = "shannon"
  ),
  Simpson = vegan::diversity(
    rare_mat,
    index = "simpson"
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 17. Rarefied alpha five-group tests
# ==============================================================================

rare_alpha_tests <- data.frame()

for (metric in c(
  "Observed_ASVs",
  "Shannon",
  "Simpson"
)) {
  
  test_i <- kruskal.test(
    alpha_rare[, metric] ~ group
  )
  
  tmp <- data.frame(
    Dataset = "Rarefied",
    Metric = metric,
    Statistic = unname(
      test_i$statistic
    ),
    df = unname(
      test_i$parameter
    ),
    P_value = test_i$p.value,
    stringsAsFactors = FALSE
  )
  
  rare_alpha_tests <- rbind(
    rare_alpha_tests,
    tmp
  )
}


# ==============================================================================
# 18. Rarefied group medians
# ==============================================================================

rare_desc <- data.frame()

for (g in GROUP_LEVELS) {
  
  idx <- group == g
  
  for (metric in c(
    "Observed_ASVs",
    "Shannon",
    "Simpson"
  )) {
    
    x <- alpha_rare[
      idx,
      metric
    ]
    
    tmp <- data.frame(
      Group = g,
      Metric = metric,
      N = length(x),
      Q1 = unname(
        quantile(
          x,
          0.25
        )
      ),
      Median = median(x),
      Q3 = unname(
        quantile(
          x,
          0.75
        )
      ),
      stringsAsFactors = FALSE
    )
    
    rare_desc <- rbind(
      rare_desc,
      tmp
    )
  }
}


# ==============================================================================
# 19. Compare CAC Shannon under native and rarefied analyses
# ==============================================================================

cac_native_shannon <- alpha_native$Shannon[
  group == "CAC"
]

cac_rare_shannon <- alpha_rare$Shannon[
  group == "CAC"
]

cac_shannon_summary <- data.frame(
  Dataset = c(
    "Native_non_rarefied",
    "Rarefied"
  ),
  N = c(
    length(cac_native_shannon),
    length(cac_rare_shannon)
  ),
  Q1 = c(
    unname(
      quantile(
        cac_native_shannon,
        0.25
      )
    ),
    unname(
      quantile(
        cac_rare_shannon,
        0.25
      )
    )
  ),
  Median = c(
    median(
      cac_native_shannon
    ),
    median(
      cac_rare_shannon
    )
  ),
  Q3 = c(
    unname(
      quantile(
        cac_native_shannon,
        0.75
      )
    ),
    unname(
      quantile(
        cac_rare_shannon,
        0.75
      )
    )
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 20. Depth-Shannon association
#
# Technical diagnostic only.
# ==============================================================================

depth_shannon_spearman <- cor.test(
  alpha_native$Sequencing_depth,
  alpha_native$Shannon,
  method = "spearman",
  exact = FALSE
)

depth_shannon_result <- data.frame(
  Analysis = "Sequencing depth vs native Shannon",
  Spearman_rho = unname(
    depth_shannon_spearman$estimate
  ),
  P_value = depth_shannon_spearman$p.value,
  stringsAsFactors = FALSE
)


# ==============================================================================
# 21. Console
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("PROGRESSION127 SEQUENCING DEPTH AUDIT\n")
cat("============================================================\n")

cat(
  "Samples = ",
  nrow(count_mat),
  "\n",
  sep = ""
)

cat(
  "ASVs = ",
  ncol(count_mat),
  "\n",
  sep = ""
)

cat(
  "Minimum sequencing depth = ",
  min(depth),
  "\n",
  sep = ""
)

cat(
  "Median sequencing depth = ",
  median(depth),
  "\n",
  sep = ""
)

cat(
  "Maximum sequencing depth = ",
  max(depth),
  "\n",
  sep = ""
)


cat("\n")
cat("============================================================\n")
cat("SEQUENCING DEPTH BY GROUP\n")
cat("============================================================\n")

print(
  depth_summary,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("SEQUENCING DEPTH FIVE-GROUP TEST\n")
cat("============================================================\n")

print(
  depth_test,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("CAC VS OTHER FOUR GROUPS DEPTH CHECK\n")
cat("============================================================\n")

print(
  depth_cac_test,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("PRIMARY/NATIVE ALPHA RECALCULATION\n")
cat("============================================================\n")

print(
  native_alpha_tests,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("EXISTING 02_ALPHA_BETA SCRIPT AUDIT\n")
cat("============================================================\n")

print(
  alpha_script_audit,
  row.names = FALSE
)

cat("\nRarefaction-related script lines:\n")

if (length(rarefaction_matching_lines) == 0) {
  
  cat("NONE\n")
  
} else {
  
  cat(
    paste(
      rarefaction_matching_lines,
      collapse = "\n"
    ),
    "\n"
  )
}

cat("\nAlpha-related script lines:\n")

if (length(alpha_matching_lines) == 0) {
  
  cat("NONE\n")
  
} else {
  
  cat(
    paste(
      alpha_matching_lines,
      collapse = "\n"
    ),
    "\n"
  )
}


cat("\n")
cat("============================================================\n")
cat("RAREFIED ALPHA SENSITIVITY\n")
cat("============================================================\n")

cat(
  "Rarefaction depth = ",
  rarefaction_depth,
  "\n",
  sep = ""
)

cat(
  "Rarefaction seed = ",
  RAREFACTION_SEED,
  "\n",
  sep = ""
)

print(
  rare_alpha_tests,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("RAREFIED ALPHA GROUP MEDIANS\n")
cat("============================================================\n")

print(
  rare_desc,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("CAC SHANNON NATIVE VS RAREFIED\n")
cat("============================================================\n")

print(
  cac_shannon_summary,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("SEQUENCING DEPTH VS SHANNON\n")
cat("============================================================\n")

print(
  depth_shannon_result,
  row.names = FALSE
)


# ==============================================================================
# 22. Save outputs
# ==============================================================================

fwrite(
  depth_sample,
  file.path(
    OUT_DIR,
    "sequencing_depth_per_sample_progression127_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  depth_summary,
  file.path(
    OUT_DIR,
    "sequencing_depth_by_group_progression127_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  depth_test,
  file.path(
    OUT_DIR,
    "sequencing_depth_five_group_Kruskal_Wallis_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  depth_cac_test,
  file.path(
    OUT_DIR,
    "sequencing_depth_CAC_vs_other_groups_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  alpha_native,
  file.path(
    OUT_DIR,
    "alpha_native_recalculated_progression127_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  native_alpha_tests,
  file.path(
    OUT_DIR,
    "alpha_native_Kruskal_Wallis_audit_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  alpha_script_audit,
  file.path(
    OUT_DIR,
    "alpha_analysis_script_rarefaction_audit_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  alpha_rare,
  file.path(
    OUT_DIR,
    "alpha_rarefied_progression127_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  rare_alpha_tests,
  file.path(
    OUT_DIR,
    "alpha_rarefied_Kruskal_Wallis_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  rare_desc,
  file.path(
    OUT_DIR,
    "alpha_rarefied_group_descriptive_statistics_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  cac_shannon_summary,
  file.path(
    OUT_DIR,
    "CAC_Shannon_native_vs_rarefied_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  depth_shannon_result,
  file.path(
    OUT_DIR,
    "sequencing_depth_vs_Shannon_Spearman_7KB.tsv"
  ),
  sep = "\t"
)

if (length(rarefaction_matching_lines) > 0) {
  
  writeLines(
    rarefaction_matching_lines,
    file.path(
      OUT_DIR,
      "02_alpha_beta_rarefaction_matching_lines_7KB.txt"
    )
  )
}

if (length(alpha_matching_lines) > 0) {
  
  writeLines(
    alpha_matching_lines,
    file.path(
      OUT_DIR,
      "02_alpha_beta_alpha_matching_lines_7KB.txt"
    )
  )
}


# ==============================================================================
# 23. Session info
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_DIR,
    "sessionInfo_sequencing_depth_alpha_audit_7KB.txt"
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