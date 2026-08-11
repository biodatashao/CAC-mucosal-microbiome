#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 05_02_core_analysis_Figure7_SuppFig3.R
##
## Module 05 - Clinical association and recurrence (Figure 7, Supplementary Figure 3)
##
## Purpose
############################################################
# Recalculate all core statistics/source tables required for:
#
# Figure 7:
#   A Tumor Shannon diversity by tumor stage
#   B Dysplasia Shannon diversity by grade
#   C Tumor Shannon diversity by recurrence
#   D Bray-Curtis PCoA by recurrence + PERMANOVA
#   E DFS according to tumor Shannon diversity
#   F Lactococcus abundance according to recurrence
#   G UCG-005 abundance according to recurrence
#   H Two-genus microbial index
#   I Strict LOOCV ROC
#
# Supplementary Figure 3:
#   recurrence abundance comparisons for the remaining 8 representative genera
#
# IMPORTANT
# ---------
# - Uses fixed progression127 7KB data.
# - CAC = the fixed 23 CAC samples in progression127.
# - Dysplasia = the fixed 17 dysplasia samples in progression127.
# - Clinical groupings come from the already audited clinical metadata tables.
# - No plotting in this script.
# ==============================================================================


options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "stringr",
  "vegan",
  "survival",
  "pROC"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(vegan)
  library(survival)
  library(pROC)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

rerun_root <- file.path(
  project_root,
  "output",
  "analysis"
)

progression_dir <- file.path(
  rerun_root,
  "00_clean_data",
  "progression127"
)

asv_count_file <- file.path(
  progression_dir,
  "asv_count_7KB_progression127.tsv"
)

taxonomy_file <- file.path(
  progression_dir,
  "taxonomy_7KB_progression127.tsv"
)

metadata_file <- file.path(
  progression_dir,
  "metadata_7KB_progression127.tsv"
)

clinical_audit_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "00_clinical_metadata_audit"
)

ca_clinical_file <- file.path(
  clinical_audit_dir,
  "CAC23_clinical_metadata_audit_7KB.tsv"
)

dys_clinical_file <- file.path(
  clinical_audit_dir,
  "Dysplasia17_clinical_metadata_audit_7KB.tsv"
)

output_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "01_core_analysis"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Output subdirectories
# ==============================================================================

alpha_dir <- file.path(
  output_dir,
  "alpha_diversity"
)

beta_dir <- file.path(
  output_dir,
  "beta_diversity"
)

survival_dir <- file.path(
  output_dir,
  "survival"
)

marker_dir <- file.path(
  output_dir,
  "marker_abundance"
)

rmi_dir <- file.path(
  output_dir,
  "two_genus_index_LOOCV"
)

supp_dir <- file.path(
  output_dir,
  "Supplementary_Figure3"
)

dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(beta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(survival_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rmi_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 3. Input checks
# ==============================================================================

required_files <- c(
  asv_count_file,
  taxonomy_file,
  metadata_file,
  ca_clinical_file,
  dys_clinical_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required file(s):\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 4. Helpers
# ==============================================================================

section <- function(title) {
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
}


extract_genus_strict <- function(x) {
  
  x <- as.character(x)
  x <- stringr::str_trim(x)
  
  genus <- rep(
    NA_character_,
    length(x)
  )
  
  has_g <- stringr::str_detect(
    x,
    "g__"
  )
  
  if (any(has_g, na.rm = TRUE)) {
    
    g_part <- stringr::str_extract(
      x[has_g],
      "g__[^;]+"
    )
    
    g_part <- stringr::str_replace(
      g_part,
      "^g__",
      ""
    )
    
    g_part <- stringr::str_trim(
      g_part
    )
    
    g_part[
      is.na(g_part) |
        g_part == ""
    ] <- NA_character_
    
    genus[has_g] <- g_part
  }
  
  no_g_idx <- which(
    !has_g |
      is.na(has_g)
  )
  
  if (length(no_g_idx) > 0) {
    
    candidate <- x[no_g_idx]
    
    candidate <- stringr::str_replace(
      candidate,
      "^g__",
      ""
    )
    
    candidate <- stringr::str_trim(
      candidate
    )
    
    looks_higher_tax <- stringr::str_detect(
      candidate,
      stringr::regex(
        "(^|;)\\s*[dkpcofs]__|k__|p__|c__|o__|f__|s__",
        ignore_case = TRUE
      )
    )
    
    candidate[looks_higher_tax] <- NA_character_
    candidate[candidate == ""] <- NA_character_
    
    genus[no_g_idx] <- candidate
  }
  
  genus
}


is_bad_genus <- function(x) {
  
  x0 <- tolower(
    as.character(x)
  )
  
  is.na(x0) |
    x0 == "" |
    x0 == "na" |
    x0 == "unassigned" |
    x0 == "unclassified" |
    x0 == "uncultured" |
    x0 == "unknown" |
    x0 == "metagenome" |
    x0 == "norank" |
    x0 == "no_rank" |
    x0 == "ambiguous" |
    stringr::str_detect(x0, "uncultured") |
    stringr::str_detect(x0, "unclassified") |
    stringr::str_detect(x0, "unknown") |
    stringr::str_detect(x0, "metagenome") |
    stringr::str_detect(x0, "norank") |
    stringr::str_detect(x0, "no_rank") |
    stringr::str_detect(x0, "ambiguous") |
    stringr::str_detect(x0, "^[dkpcofs]__") |
    stringr::str_detect(x0, "k__|p__|c__|o__|f__|s__")
}


canonicalize_genus <- function(x) {
  
  x <- as.character(x)
  
  x[x == "UCG-005"] <- "UCG_005"
  x[x == "UCG-009"] <- "UCG_009"
  x[x == "Escherichia-Shigella"] <- "Escherichia_Shigella"
  
  x[
    x == "Christensenellaceae_R-7_group"
  ] <- "Christensenellaceae_R_7_group"
  
  x[
    x == "f__Eubacterium__eligens_group"
  ] <- "Eubacterium_eligens_group"
  
  x[
    x == "[Eubacterium]_eligens_group"
  ] <- "Eubacterium_eligens_group"
  
  x[
    x == "[Eubacterium] eligens group"
  ] <- "Eubacterium_eligens_group"
  
  x
}


safe_wilcox <- function(value, group) {
  
  keep <- is.finite(value) &
    !is.na(group)
  
  value <- value[keep]
  group <- droplevels(
    factor(group[keep])
  )
  
  if (length(levels(group)) != 2) {
    return(NA_real_)
  }
  
  suppressWarnings(
    stats::wilcox.test(
      value ~ group,
      exact = FALSE
    )$p.value
  )
}


safe_kruskal <- function(value, group) {
  
  keep <- is.finite(value) &
    !is.na(group)
  
  value <- value[keep]
  group <- droplevels(
    factor(group[keep])
  )
  
  if (length(levels(group)) < 2) {
    return(NA_real_)
  }
  
  stats::kruskal.test(
    value ~ group
  )$p.value
}


zscore_using_training <- function(
    train_values,
    target_values
) {
  
  train_mean <- mean(
    train_values,
    na.rm = TRUE
  )
  
  train_sd <- stats::sd(
    train_values,
    na.rm = TRUE
  )
  
  if (
    !is.finite(train_sd) ||
    train_sd == 0
  ) {
    stop(
      "Training-set SD is zero or invalid during LOOCV.",
      call. = FALSE
    )
  }
  
  (
    target_values -
      train_mean
  ) /
    train_sd
}


# ==============================================================================
# 5. Read progression127
# ==============================================================================

section(
  "READ PROGRESSION127"
)

asv_count <- data.table::fread(
  asv_count_file,
  data.table = FALSE,
  check.names = FALSE
)

taxonomy <- data.table::fread(
  taxonomy_file,
  data.table = FALSE,
  check.names = FALSE
)

metadata <- data.table::fread(
  metadata_file,
  data.table = FALSE,
  check.names = FALSE
)


if (!"ASV" %in% colnames(asv_count)) {
  stop(
    "ASV count table does not contain ASV.",
    call. = FALSE
  )
}

if (
  !all(
    c(
      "ASV",
      "Taxonomy"
    ) %in%
    colnames(taxonomy)
  )
) {
  stop(
    "Taxonomy table must contain ASV and Taxonomy.",
    call. = FALSE
  )
}

if (!"SampleID" %in% colnames(metadata)) {
  stop(
    "Metadata does not contain SampleID.",
    call. = FALSE
  )
}

if (!"Progression5" %in% colnames(metadata)) {
  stop(
    "Metadata does not contain Progression5.",
    call. = FALSE
  )
}


metadata$SampleID <- as.character(
  metadata$SampleID
)

metadata$Progression5 <- as.character(
  metadata$Progression5
)


cac_samples <- metadata$SampleID[
  metadata$Progression5 %in%
    c(
      "CA",
      "CAC"
    )
]

dysplasia_samples <- metadata$SampleID[
  metadata$Progression5 == "Dysplasia"
]


if (length(cac_samples) != 23) {
  stop(
    paste0(
      "Expected CAC23, found ",
      length(cac_samples),
      "."
    ),
    call. = FALSE
  )
}

if (length(dysplasia_samples) != 17) {
  stop(
    paste0(
      "Expected Dysplasia17, found ",
      length(dysplasia_samples),
      "."
    ),
    call. = FALSE
  )
}


all_target_samples <- c(
  cac_samples,
  dysplasia_samples
)


missing_target_samples <- setdiff(
  all_target_samples,
  colnames(asv_count)
)

if (length(missing_target_samples) > 0) {
  stop(
    paste0(
      "Target samples missing from count table:\n",
      paste(missing_target_samples, collapse = "\n")
    ),
    call. = FALSE
  )
}


cat(
  "CAC samples: ",
  length(cac_samples),
  "\n",
  sep = ""
)

cat(
  "Dysplasia samples: ",
  length(dysplasia_samples),
  "\n",
  sep = ""
)


# ==============================================================================
# 6. ASV raw count matrix
# ==============================================================================

section(
  "BUILD RAW ASV COUNT MATRIX"
)

asv_matrix <- as.matrix(
  asv_count[
    ,
    all_target_samples,
    drop = FALSE
  ]
)

storage.mode(
  asv_matrix
) <- "numeric"

rownames(
  asv_matrix
) <- as.character(
  asv_count$ASV
)


if (
  any(
    !is.finite(
      asv_matrix
    )
  )
) {
  stop(
    "Non-finite values detected in ASV count matrix.",
    call. = FALSE
  )
}

if (
  any(
    asv_matrix < 0
  )
) {
  stop(
    "Negative ASV counts detected.",
    call. = FALSE
  )
}


sample_depth <- colSums(
  asv_matrix
)

if (
  any(
    sample_depth <= 0
  )
) {
  stop(
    "At least one target sample has zero total counts.",
    call. = FALSE
  )
}


# ==============================================================================
# 7. Shannon diversity
#
# Shannon is calculated directly from ASV raw counts.
# ==============================================================================

section(
  "CALCULATE ASV-LEVEL SHANNON DIVERSITY"
)

shannon_values <- vegan::diversity(
  t(
    asv_matrix
  ),
  index = "shannon"
)

shannon_table <- data.frame(
  SampleID = names(
    shannon_values
  ),
  Shannon = as.numeric(
    shannon_values
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 8. Read audited clinical metadata
# ==============================================================================

section(
  "READ AUDITED CLINICAL METADATA"
)

ca_clinical <- data.table::fread(
  ca_clinical_file,
  data.table = FALSE,
  check.names = FALSE
)

dys_clinical <- data.table::fread(
  dys_clinical_file,
  data.table = FALSE,
  check.names = FALSE
)


required_ca_columns <- c(
  "SampleID",
  "Tumor_stage_plot",
  "Recurrence_plot",
  "DFS_time_months_raw",
  "DFS_status_raw",
  "DFS_valid"
)

missing_ca_columns <- setdiff(
  required_ca_columns,
  colnames(ca_clinical)
)

if (length(missing_ca_columns) > 0) {
  stop(
    paste0(
      "CAC audit table missing: ",
      paste(missing_ca_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}


required_dys_columns <- c(
  "SampleID",
  "Dysplasia_grade_plot"
)

missing_dys_columns <- setdiff(
  required_dys_columns,
  colnames(dys_clinical)
)

if (length(missing_dys_columns) > 0) {
  stop(
    paste0(
      "Dysplasia audit table missing: ",
      paste(missing_dys_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 9. CAC23 Shannon clinical table
# ==============================================================================

ca_shannon <- merge(
  ca_clinical,
  shannon_table,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

if (
  nrow(ca_shannon) != 23 ||
  anyNA(ca_shannon$Shannon)
) {
  stop(
    "CAC23 Shannon merge failed.",
    call. = FALSE
  )
}


ca_shannon$Tumor_stage_plot <- factor(
  ca_shannon$Tumor_stage_plot,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III"
  )
)

ca_shannon$Recurrence_plot <- factor(
  ca_shannon$Recurrence_plot,
  levels = c(
    "No recurrence",
    "Recurrence"
  )
)


data.table::fwrite(
  ca_shannon,
  file.path(
    alpha_dir,
    "Figure7_CAC23_Shannon_clinical_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 10. Dysplasia17 Shannon table
# ==============================================================================

dys_shannon <- merge(
  dys_clinical,
  shannon_table,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

if (
  nrow(dys_shannon) != 17 ||
  anyNA(dys_shannon$Shannon)
) {
  stop(
    "Dysplasia17 Shannon merge failed.",
    call. = FALSE
  )
}


dys_shannon$Dysplasia_grade_plot <- factor(
  dys_shannon$Dysplasia_grade_plot,
  levels = c(
    "High-grade dysplasia",
    "Low-grade dysplasia"
  )
)


data.table::fwrite(
  dys_shannon,
  file.path(
    alpha_dir,
    "Figure7_Dysplasia17_Shannon_grade_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 11. Figure 7A/B/C statistics
# ==============================================================================

section(
  "FIGURE 7 A/B/C SHANNON TESTS"
)


p_stage <- safe_kruskal(
  ca_shannon$Shannon,
  ca_shannon$Tumor_stage_plot
)

p_grade <- safe_wilcox(
  dys_shannon$Shannon,
  dys_shannon$Dysplasia_grade_plot
)

p_recurrence_shannon <- safe_wilcox(
  ca_shannon$Shannon,
  ca_shannon$Recurrence_plot
)


alpha_stats <- data.frame(
  Panel = c(
    "Figure7A",
    "Figure7B",
    "Figure7C"
  ),
  
  Analysis = c(
    "CAC Shannon by tumor stage",
    "Dysplasia Shannon by grade",
    "CAC Shannon by recurrence"
  ),
  
  Test = c(
    "Kruskal-Wallis",
    "Wilcoxon rank-sum",
    "Wilcoxon rank-sum"
  ),
  
  P_value = c(
    p_stage,
    p_grade,
    p_recurrence_shannon
  ),
  
  stringsAsFactors = FALSE
)


print(
  alpha_stats,
  row.names = FALSE
)


data.table::fwrite(
  alpha_stats,
  file.path(
    alpha_dir,
    "Figure7_ABC_Shannon_statistics_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 12. Recurrence Bray-Curtis PCoA
#
# CAC23 only.
# Bray-Curtis is calculated from sample-level ASV relative abundance.
# ==============================================================================

section(
  "FIGURE 7D RECURRENCE BRAY-CURTIS"
)


cac_asv_matrix <- asv_matrix[
  ,
  cac_samples,
  drop = FALSE
]


cac_relative <- sweep(
  cac_asv_matrix,
  2,
  colSums(cac_asv_matrix),
  "/"
)


cac_relative_sample <- t(
  cac_relative
)


recurrence_by_sample <- ca_shannon$Recurrence_plot[
  match(
    rownames(cac_relative_sample),
    ca_shannon$SampleID
  )
]


if (
  anyNA(
    recurrence_by_sample
  )
) {
  stop(
    "Recurrence labels could not be aligned to CAC Bray-Curtis matrix.",
    call. = FALSE
  )
}


bray_distance <- vegan::vegdist(
  cac_relative_sample,
  method = "bray"
)


pcoa_result <- stats::cmdscale(
  bray_distance,
  k = 2,
  eig = TRUE,
  add = TRUE
)


pcoa_coordinates <- as.data.frame(
  pcoa_result$points
)

colnames(
  pcoa_coordinates
) <- c(
  "PCoA1",
  "PCoA2"
)

pcoa_coordinates$SampleID <- rownames(
  pcoa_coordinates
)

pcoa_coordinates$Recurrence_plot <- recurrence_by_sample[
  match(
    pcoa_coordinates$SampleID,
    rownames(cac_relative_sample)
  )
]


positive_eigenvalues <- pcoa_result$eig[
  pcoa_result$eig > 0
]

pc1_percent <- 100 *
  pcoa_result$eig[1] /
  sum(
    positive_eigenvalues
  )

pc2_percent <- 100 *
  pcoa_result$eig[2] /
  sum(
    positive_eigenvalues
  )


set.seed(
  20260726
)

permanova_result <- vegan::adonis2(
  bray_distance ~ recurrence_by_sample,
  permutations = 9999
)


permanova_r2 <- permanova_result$R2[1]
permanova_p <- permanova_result$`Pr(>F)`[1]


beta_summary <- data.frame(
  Metric = c(
    "PCoA1_variance_percent",
    "PCoA2_variance_percent",
    "PERMANOVA_R2",
    "PERMANOVA_P",
    "PERMANOVA_permutations"
  ),
  
  Value = c(
    pc1_percent,
    pc2_percent,
    permanova_r2,
    permanova_p,
    9999
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  pcoa_coordinates,
  file.path(
    beta_dir,
    "Figure7D_recurrence_Bray_PCoA_coordinates_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  beta_summary,
  file.path(
    beta_dir,
    "Figure7D_recurrence_Bray_PERMANOVA_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 13. DFS analysis according to CAC Shannon
#
# Median split across all 23 CAC samples.
# ==============================================================================

section(
  "FIGURE 7E DFS BY TUMOR SHANNON"
)


dfs_data <- ca_shannon[
  ca_shannon$DFS_valid %in%
    c(
      TRUE,
      "TRUE"
    ),
  ,
  drop = FALSE
]


if (
  nrow(
    dfs_data
  ) != 23
) {
  stop(
    paste0(
      "Expected 23 valid DFS records, found ",
      nrow(dfs_data),
      "."
    ),
    call. = FALSE
  )
}


shannon_cutoff <- stats::median(
  dfs_data$Shannon,
  na.rm = TRUE
)


dfs_data$Shannon_group <- ifelse(
  dfs_data$Shannon <=
    shannon_cutoff,
  "Low Shannon",
  "High Shannon"
)


dfs_data$Shannon_group <- factor(
  dfs_data$Shannon_group,
  levels = c(
    "Low Shannon",
    "High Shannon"
  )
)


dfs_status_text <- tolower(
  trimws(
    as.character(
      dfs_data$DFS_status_raw
    )
  )
)


dfs_event <- rep(
  NA_integer_,
  length(
    dfs_status_text
  )
)


dfs_event[
  dfs_status_text %in%
    c(
      "1",
      "recurrence",
      "recurred",
      "event",
      "yes"
    )
] <- 1L


dfs_event[
  dfs_status_text %in%
    c(
      "0",
      "no recurrence",
      "non-recurrence",
      "nonrecurrence",
      "censored",
      "no"
    )
] <- 0L


fallback_event <- as.character(
  dfs_data$Recurrence_plot
) == "Recurrence"


dfs_event[
  is.na(
    dfs_event
  )
] <- as.integer(
  fallback_event[
    is.na(
      dfs_event
    )
  ]
)


if (
  anyNA(
    dfs_event
  )
) {
  stop(
    "Could not define DFS event for all CAC samples.",
    call. = FALSE
  )
}


dfs_data$DFS_event <- dfs_event

dfs_data$DFS_time_months <- as.numeric(
  dfs_data$DFS_time_months_raw
)


if (
  any(
    !is.finite(
      dfs_data$DFS_time_months
    )
  )
) {
  stop(
    "Invalid DFS time detected.",
    call. = FALSE
  )
}


surv_object <- survival::Surv(
  time = dfs_data$DFS_time_months,
  event = dfs_data$DFS_event
)


km_fit <- survival::survfit(
  surv_object ~ Shannon_group,
  data = dfs_data
)


logrank_test <- survival::survdiff(
  surv_object ~ Shannon_group,
  data = dfs_data
)


logrank_p <- stats::pchisq(
  logrank_test$chisq,
  df = length(logrank_test$n) - 1,
  lower.tail = FALSE
)


dfs_summary <- data.frame(
  Metric = c(
    "Shannon_median_cutoff",
    "Low_Shannon_N",
    "High_Shannon_N",
    "Low_Shannon_events",
    "High_Shannon_events",
    "Logrank_chisq",
    "Logrank_P"
  ),
  
  Value = c(
    shannon_cutoff,
    sum(
      dfs_data$Shannon_group ==
        "Low Shannon"
    ),
    sum(
      dfs_data$Shannon_group ==
        "High Shannon"
    ),
    sum(
      dfs_data$DFS_event[
        dfs_data$Shannon_group ==
          "Low Shannon"
      ]
    ),
    sum(
      dfs_data$DFS_event[
        dfs_data$Shannon_group ==
          "High Shannon"
      ]
    ),
    logrank_test$chisq,
    logrank_p
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  dfs_data,
  file.path(
    survival_dir,
    "Figure7E_DFS_Shannon_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  dfs_summary,
  file.path(
    survival_dir,
    "Figure7E_DFS_Shannon_statistics_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


saveRDS(
  km_fit,
  file.path(
    survival_dir,
    "Figure7E_DFS_Shannon_survfit_7KB.rds"
  )
)


# ==============================================================================
# 14. Collapse ASV counts to strict genus level
# ==============================================================================

section(
  "COLLAPSE ASV COUNTS TO GENUS"
)


taxonomy_match <- match(
  as.character(
    asv_count$ASV
  ),
  as.character(
    taxonomy$ASV
  )
)


if (
  anyNA(
    taxonomy_match
  )
) {
  stop(
    "Some ASVs lack taxonomy.",
    call. = FALSE
  )
}


taxonomy_aligned <- taxonomy[
  taxonomy_match,
  ,
  drop = FALSE
]


genus_vector <- extract_genus_strict(
  taxonomy_aligned$Taxonomy
)

genus_vector <- canonicalize_genus(
  genus_vector
)

genus_vector[
  is_bad_genus(
    genus_vector
  )
] <- NA_character_


valid_genus_rows <- !is.na(
  genus_vector
)


genus_count_matrix <- rowsum(
  asv_matrix[
    valid_genus_rows,
    ,
    drop = FALSE
  ],
  group = genus_vector[
    valid_genus_rows
  ],
  reorder = FALSE
)


genus_relative_matrix <- sweep(
  genus_count_matrix,
  2,
  colSums(genus_count_matrix),
  "/"
)


cat(
  "Strict genera available: ",
  nrow(genus_count_matrix),
  "\n",
  sep = ""
)


# ==============================================================================
# 15. Target genera
# ==============================================================================

figure7_marker_genera <- c(
  "Lactococcus",
  "UCG_005"
)


supplementary_genera <- c(
  "Mediterraneibacter",
  "Candidatus_Soleaferrea",
  "Peptoclostridium",
  "Monoglobus",
  "Dorea",
  "Atopostipes",
  "Eubacterium_eligens_group",
  "Desulfovibrio"
)


all_marker_genera <- c(
  figure7_marker_genera,
  supplementary_genera
)


missing_marker_genera <- setdiff(
  all_marker_genera,
  rownames(
    genus_relative_matrix
  )
)


if (
  length(
    missing_marker_genera
  ) > 0
) {
  
  stop(
    paste0(
      "Required genus/genus labels absent from the 7KB strict-genus matrix:\n",
      paste(
        missing_marker_genera,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 16. Build CAC23 genus abundance source
# ==============================================================================

cac_genus_relative <- genus_relative_matrix[
  all_marker_genera,
  cac_samples,
  drop = FALSE
]


cac_marker_long <- as.data.frame(
  as.table(
    cac_genus_relative
  ),
  stringsAsFactors = FALSE
)


colnames(
  cac_marker_long
) <- c(
  "Genus",
  "SampleID",
  "Relative_abundance"
)


cac_marker_long$Relative_abundance <- as.numeric(
  cac_marker_long$Relative_abundance
)


cac_marker_long <- merge(
  cac_marker_long,
  ca_shannon[
    ,
    c(
      "SampleID",
      "Recurrence_plot"
    ),
    drop = FALSE
  ],
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)


cac_marker_long$Relative_abundance_percent <-
  100 *
  cac_marker_long$Relative_abundance


if (
  anyNA(
    cac_marker_long$Recurrence_plot
  )
) {
  stop(
    "Recurrence labels missing after genus-abundance merge.",
    call. = FALSE
  )
}


data.table::fwrite(
  cac_marker_long,
  file.path(
    marker_dir,
    "Figure7_SuppFig3_CAC23_marker_abundance_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 17. Figure 7F/G abundance tests
# ==============================================================================

section(
  "FIGURE 7F/G MARKER ABUNDANCE"
)


fg_stats <- do.call(
  rbind,
  lapply(
    figure7_marker_genera,
    function(current_genus) {
      
      current <- cac_marker_long[
        cac_marker_long$Genus ==
          current_genus,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Genus = current_genus,
        
        No_recurrence_N = sum(
          current$Recurrence_plot ==
            "No recurrence"
        ),
        
        Recurrence_N = sum(
          current$Recurrence_plot ==
            "Recurrence"
        ),
        
        No_recurrence_median_percent = stats::median(
          current$Relative_abundance_percent[
            current$Recurrence_plot ==
              "No recurrence"
          ]
        ),
        
        Recurrence_median_percent = stats::median(
          current$Relative_abundance_percent[
            current$Recurrence_plot ==
              "Recurrence"
          ]
        ),
        
        Wilcoxon_P = safe_wilcox(
          current$Relative_abundance_percent,
          current$Recurrence_plot
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


print(
  fg_stats,
  row.names = FALSE
)


data.table::fwrite(
  fg_stats,
  file.path(
    marker_dir,
    "Figure7_FG_marker_abundance_statistics_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 18. Supplementary Figure 3 abundance tests
# ==============================================================================

section(
  "SUPPLEMENTARY FIGURE 3 MARKER ABUNDANCE"
)


supp_stats <- do.call(
  rbind,
  lapply(
    supplementary_genera,
    function(current_genus) {
      
      current <- cac_marker_long[
        cac_marker_long$Genus ==
          current_genus,
        ,
        drop = FALSE
      ]
      
      data.frame(
        Genus = current_genus,
        
        No_recurrence_N = sum(
          current$Recurrence_plot ==
            "No recurrence"
        ),
        
        Recurrence_N = sum(
          current$Recurrence_plot ==
            "Recurrence"
        ),
        
        No_recurrence_median_percent = stats::median(
          current$Relative_abundance_percent[
            current$Recurrence_plot ==
              "No recurrence"
          ]
        ),
        
        Recurrence_median_percent = stats::median(
          current$Relative_abundance_percent[
            current$Recurrence_plot ==
              "Recurrence"
          ]
        ),
        
        Wilcoxon_P = safe_wilcox(
          current$Relative_abundance_percent,
          current$Recurrence_plot
        ),
        
        stringsAsFactors = FALSE
      )
    }
  )
)


print(
  supp_stats,
  row.names = FALSE
)


data.table::fwrite(
  supp_stats,
  file.path(
    supp_dir,
    "Supplementary_Figure3_marker_statistics_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  cac_marker_long[
    cac_marker_long$Genus %in%
      supplementary_genera,
    ,
    drop = FALSE
  ],
  file.path(
    supp_dir,
    "Supplementary_Figure3_marker_abundance_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 19. Two-genus microbial index
#
# Global plotting index:
#   z(UCG_005) - z(Lactococcus)
#
# This global index is used only for Figure 7H visualization.
# Strict LOOCV below recalculates z scores from TRAINING samples only.
# ==============================================================================

section(
  "FIGURE 7H TWO-GENUS MICROBIAL INDEX"
)


lactococcus_vector <- genus_relative_matrix[
  "Lactococcus",
  cac_samples
]

ucg005_vector <- genus_relative_matrix[
  "UCG_005",
  cac_samples
]


rmi_table <- data.frame(
  SampleID = cac_samples,
  
  Recurrence_plot = as.character(
    ca_shannon$Recurrence_plot[
      match(
        cac_samples,
        ca_shannon$SampleID
      )
    ]
  ),
  
  Lactococcus = as.numeric(
    lactococcus_vector
  ),
  
  UCG_005 = as.numeric(
    ucg005_vector
  ),
  
  stringsAsFactors = FALSE
)


rmi_table$Lactococcus_z_global <- as.numeric(
  scale(
    rmi_table$Lactococcus
  )
)

rmi_table$UCG_005_z_global <- as.numeric(
  scale(
    rmi_table$UCG_005
  )
)

rmi_table$RMI_two_genus_global <-
  rmi_table$UCG_005_z_global -
  rmi_table$Lactococcus_z_global


rmi_table$Recurrence_plot <- factor(
  rmi_table$Recurrence_plot,
  levels = c(
    "No recurrence",
    "Recurrence"
  )
)


rmi_wilcoxon_p <- safe_wilcox(
  rmi_table$RMI_two_genus_global,
  rmi_table$Recurrence_plot
)


rmi_stats <- data.frame(
  Analysis = "Two-genus microbial index by recurrence",
  
  Definition = "z(UCG_005) - z(Lactococcus)",
  
  Wilcoxon_P = rmi_wilcoxon_p,
  
  No_recurrence_median = stats::median(
    rmi_table$RMI_two_genus_global[
      rmi_table$Recurrence_plot ==
        "No recurrence"
    ]
  ),
  
  Recurrence_median = stats::median(
    rmi_table$RMI_two_genus_global[
      rmi_table$Recurrence_plot ==
        "Recurrence"
    ]
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  rmi_table,
  file.path(
    rmi_dir,
    "Figure7H_two_genus_index_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  rmi_stats,
  file.path(
    rmi_dir,
    "Figure7H_two_genus_index_statistics_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 20. Strict LOOCV logistic prediction
#
# For each held-out CAC sample:
#
#   1. Use remaining 22 samples as training data.
#   2. Calculate training mean/SD for UCG_005 and Lactococcus.
#   3. Standardize both training and held-out sample using TRAINING parameters.
#   4. Training RMI = z(UCG_005) - z(Lactococcus).
#   5. Fit logistic regression on the training set.
#   6. Predict recurrence probability for the held-out sample.
#
# This avoids leakage from global z scoring.
# ==============================================================================

section(
  "FIGURE 7I STRICT LOOCV ROC"
)


loocv_input <- rmi_table

loocv_input$Outcome <- as.integer(
  loocv_input$Recurrence_plot ==
    "Recurrence"
)


n_samples <- nrow(
  loocv_input
)


loocv_result <- data.frame(
  SampleID = loocv_input$SampleID,
  Outcome = loocv_input$Outcome,
  Recurrence_plot = as.character(
    loocv_input$Recurrence_plot
  ),
  LOOCV_RMI = NA_real_,
  LOOCV_predicted_probability = NA_real_,
  stringsAsFactors = FALSE
)


for (i in seq_len(n_samples)) {
  
  train_index <- setdiff(
    seq_len(n_samples),
    i
  )
  
  test_index <- i
  
  
  train_lacto <- loocv_input$Lactococcus[
    train_index
  ]
  
  train_ucg <- loocv_input$UCG_005[
    train_index
  ]
  
  
  train_lacto_z <- zscore_using_training(
    train_values = train_lacto,
    target_values = train_lacto
  )
  
  train_ucg_z <- zscore_using_training(
    train_values = train_ucg,
    target_values = train_ucg
  )
  
  
  test_lacto_z <- zscore_using_training(
    train_values = train_lacto,
    target_values = loocv_input$Lactococcus[
      test_index
    ]
  )
  
  test_ucg_z <- zscore_using_training(
    train_values = train_ucg,
    target_values = loocv_input$UCG_005[
      test_index
    ]
  )
  
  
  train_rmi <-
    train_ucg_z -
    train_lacto_z
  
  
  test_rmi <-
    test_ucg_z -
    test_lacto_z
  
  
  train_df <- data.frame(
    Outcome = loocv_input$Outcome[
      train_index
    ],
    RMI = train_rmi
  )
  
  
  fit <- suppressWarnings(
    stats::glm(
      Outcome ~ RMI,
      data = train_df,
      family = stats::binomial()
    )
  )
  
  
  prediction <- suppressWarnings(
    stats::predict(
      fit,
      newdata = data.frame(
        RMI = test_rmi
      ),
      type = "response"
    )
  )
  
  
  loocv_result$LOOCV_RMI[i] <- test_rmi
  loocv_result$LOOCV_predicted_probability[i] <- prediction
}


if (
  any(
    !is.finite(
      loocv_result$LOOCV_predicted_probability
    )
  )
) {
  stop(
    "Invalid LOOCV predicted probability detected.",
    call. = FALSE
  )
}


roc_object <- pROC::roc(
  response = loocv_result$Outcome,
  predictor = loocv_result$LOOCV_predicted_probability,
  levels = c(
    0,
    1
  ),
  direction = "<",
  quiet = TRUE
)


loocv_auc <- as.numeric(
  pROC::auc(
    roc_object
  )
)


loocv_auc_ci <- as.numeric(
  pROC::ci.auc(
    roc_object,
    method = "delong"
  )
)


roc_coordinates <- pROC::coords(
  roc_object,
  x = "all",
  ret = c(
    "specificity",
    "sensitivity",
    "threshold"
  ),
  transpose = FALSE
)


roc_coordinates <- as.data.frame(
  roc_coordinates
)


roc_coordinates$One_minus_specificity <-
  1 -
  roc_coordinates$specificity


loocv_summary <- data.frame(
  N = n_samples,
  
  No_recurrence_N = sum(
    loocv_result$Outcome == 0
  ),
  
  Recurrence_N = sum(
    loocv_result$Outcome == 1
  ),
  
  AUC = loocv_auc,
  
  AUC_CI_lower = loocv_auc_ci[1],
  
  AUC_CI_median = loocv_auc_ci[2],
  
  AUC_CI_upper = loocv_auc_ci[3],
  
  CI_method = "DeLong",
  
  Predictor = "Strict LOOCV two-genus microbial index",
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  loocv_result,
  file.path(
    rmi_dir,
    "Figure7I_strict_LOOCV_predictions_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  roc_coordinates,
  file.path(
    rmi_dir,
    "Figure7I_strict_LOOCV_ROC_coordinates_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  loocv_summary,
  file.path(
    rmi_dir,
    "Figure7I_strict_LOOCV_AUC_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


saveRDS(
  roc_object,
  file.path(
    rmi_dir,
    "Figure7I_strict_LOOCV_ROC_7KB.rds"
  )
)


# ==============================================================================
# 21. Master result summary
# ==============================================================================

section(
  "MASTER RESULT SUMMARY"
)


master_summary <- data.frame(
  Result = c(
    "Figure7A_stage_Shannon_P",
    "Figure7B_grade_Shannon_P",
    "Figure7C_recurrence_Shannon_P",
    "Figure7D_PERMANOVA_R2",
    "Figure7D_PERMANOVA_P",
    "Figure7D_PCoA1_percent",
    "Figure7D_PCoA2_percent",
    "Figure7E_Shannon_median_cutoff",
    "Figure7E_Logrank_P",
    "Figure7F_Lactococcus_P",
    "Figure7G_UCG005_P",
    "Figure7H_RMI_P",
    "Figure7I_LOOCV_AUC",
    "Figure7I_AUC_CI_lower",
    "Figure7I_AUC_CI_upper"
  ),
  
  Value = c(
    p_stage,
    p_grade,
    p_recurrence_shannon,
    permanova_r2,
    permanova_p,
    pc1_percent,
    pc2_percent,
    shannon_cutoff,
    logrank_p,
    fg_stats$Wilcoxon_P[
      fg_stats$Genus ==
        "Lactococcus"
    ],
    fg_stats$Wilcoxon_P[
      fg_stats$Genus ==
        "UCG_005"
    ],
    rmi_wilcoxon_p,
    loocv_auc,
    loocv_auc_ci[1],
    loocv_auc_ci[3]
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  master_summary,
  file.path(
    output_dir,
    "Figure7_master_statistics_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


print(
  master_summary,
  row.names = FALSE
)


cat("\nSupplementary Figure 3 statistics:\n")

print(
  supp_stats,
  row.names = FALSE
)


cat("\n")
cat("Output directory:\n")
cat(
  output_dir,
  "\n"
)