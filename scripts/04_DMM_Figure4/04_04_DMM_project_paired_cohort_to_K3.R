#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_04_DMM_project_paired_cohort_to_K3.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Purpose
############################################################
# Project the fixed paired nonCAC23 + CAC23 samples onto the already fitted
# progression127 K=3 DMM model.
#
# This script DOES NOT refit the DMM.
#
# Workflow
# --------
# 1. Read the fixed 7KB K=3 DMM model.
# 2. Extract fixed Dirichlet-multinomial alpha parameters and mixture weights.
# 3. Read the fixed paired 23 nonCAC + 23 CAC dataset.
# 4. Collapse ASVs to strict genus level using the same genus parser used for
#    the primary 7KB DMM analysis.
# 5. Align paired genus counts to the exact fixed-model genus feature space.
# 6. Calculate fixed-model posterior probabilities.
# 7. Recover internal model component -> C1/C2/C3 mapping from the original
#    127 training samples.
# 8. Assign all paired samples to C1/C2/C3.
# 9. Validate the 23 CAC projections against their original primary K=3 states.
# 10. Generate the paired nonCAC -> CAC transition table for Figure 4E.
#
# IMPORTANT
# ---------
# - Raw integer counts are used.
# - No DMM refitting is performed.
# - Missing model genera in paired samples are filled with zero.
# - Paired genera outside the fixed model feature space are excluded only from
#   projection.
# ==============================================================================


options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "stringr",
  "DirichletMultinomial"
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
  library(stringr)
  library(DirichletMultinomial)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

primary_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "DMM_progression127_primary_prevalence10_7KB"
)

paired_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "00_clean_data",
  "CA23_nonCA23_paired"
)

model_file <- file.path(
  primary_dir,
  "models",
  "DMM_best_model_K3_7KB.rds"
)

primary_assignment_file <- file.path(
  primary_dir,
  "tables",
  "DMM_sample_assignments_optimalK3_native_7KB.tsv"
)

asv_count_file <- file.path(
  paired_dir,
  "asv_count_7KB_CA23_nonCA23_paired.tsv"
)

taxonomy_file <- file.path(
  paired_dir,
  "taxonomy_7KB_CA23_nonCA23_paired.tsv"
)

metadata_file <- file.path(
  paired_dir,
  "metadata_7KB_CA23_nonCA23_paired.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "projection_CA23_nonCA23_fixed_K3_7KB"
)

table_dir <- file.path(
  output_dir,
  "tables"
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Output files
# ==============================================================================

projection_all_file <- file.path(
  table_dir,
  "DMM_fixed_K3_projection_all_CA_nonCA_7KB.tsv"
)

projection_compatibility_file <- file.path(
  table_dir,
  "DMM_fixed_K3_projection_all_CA_nonCA.tsv"
)

projection_nonca_file <- file.path(
  table_dir,
  "DMM_fixed_K3_projection_nonCA23_7KB.tsv"
)

projection_ca_file <- file.path(
  table_dir,
  "DMM_fixed_K3_projection_CA23_7KB.tsv"
)

feature_alignment_file <- file.path(
  table_dir,
  "DMM_fixed_K3_feature_alignment_summary_7KB.tsv"
)

sample_retention_file <- file.path(
  table_dir,
  "DMM_fixed_K3_sample_count_retention_7KB.tsv"
)

component_mapping_file <- file.path(
  table_dir,
  "DMM_internal_component_to_C1_C2_C3_7KB.tsv"
)

component_mapping_contingency_file <- file.path(
  table_dir,
  "DMM_component_mapping_contingency_7KB.tsv"
)

cluster_distribution_file <- file.path(
  table_dir,
  "DMM_fixed_K3_cluster_distribution_by_group_7KB.tsv"
)

confidence_summary_file <- file.path(
  table_dir,
  "DMM_fixed_K3_assignment_confidence_summary_7KB.tsv"
)

ca_validation_file <- file.path(
  table_dir,
  "DMM_fixed_K3_CA_validation_against_primary_7KB.tsv"
)

ca_validation_summary_file <- file.path(
  table_dir,
  "DMM_fixed_K3_CA_validation_summary_7KB.tsv"
)

pair_transition_file <- file.path(
  table_dir,
  "DMM_fixed_K3_paired_state_transition_source_7KB.tsv"
)

pair_transition_matrix_file <- file.path(
  table_dir,
  "DMM_fixed_K3_paired_state_transition_matrix_7KB.tsv"
)


# ==============================================================================
# 3. Helper functions
# ==============================================================================

print_section <- function(title) {
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
}


normalize_text <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x[
    is.na(x) |
      x %in% c("", "NA", "NaN", "NULL")
  ] <- NA_character_
  
  x
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


standardize_cluster_label <- function(x) {
  
  x <- normalize_text(x)
  
  output <- rep(
    NA_character_,
    length(x)
  )
  
  valid <- !is.na(x) &
    grepl(
      "[123]",
      x
    )
  
  output[valid] <- paste0(
    "C",
    sub(
      ".*?([123]).*",
      "\\1",
      x[valid]
    )
  )
  
  output[
    !output %in% c("C1", "C2", "C3")
  ] <- NA_character_
  
  output
}


all_permutations_3 <- function() {
  
  rbind(
    c("C1", "C2", "C3"),
    c("C1", "C3", "C2"),
    c("C2", "C1", "C3"),
    c("C2", "C3", "C1"),
    c("C3", "C1", "C2"),
    c("C3", "C2", "C1")
  )
}


derive_component_mapping <- function(
    internal_component,
    published_cluster
) {
  
  mapping_input <- data.frame(
    Internal_component = as.integer(
      internal_component
    ),
    Published_cluster = standardize_cluster_label(
      published_cluster
    ),
    stringsAsFactors = FALSE
  )
  
  mapping_input <- mapping_input[
    !is.na(mapping_input$Internal_component) &
      !is.na(mapping_input$Published_cluster),
    ,
    drop = FALSE
  ]
  
  if (nrow(mapping_input) == 0) {
    stop(
      "No valid training assignments available for component mapping.",
      call. = FALSE
    )
  }
  
  contingency <- table(
    factor(
      mapping_input$Internal_component,
      levels = 1:3
    ),
    factor(
      mapping_input$Published_cluster,
      levels = c("C1", "C2", "C3")
    )
  )
  
  permutations <- all_permutations_3()
  
  mapping_scores <- vapply(
    seq_len(nrow(permutations)),
    function(i) {
      
      current_mapping <- permutations[i, ]
      
      sum(
        contingency[
          cbind(
            1:3,
            match(
              current_mapping,
              colnames(contingency)
            )
          )
        ]
      )
    },
    FUN.VALUE = numeric(1)
  )
  
  best_index <- which.max(
    mapping_scores
  )
  
  best_mapping <- as.character(
    permutations[best_index, ]
  )
  
  matched_by_component <- vapply(
    1:3,
    function(component_number) {
      
      as.integer(
        contingency[
          component_number,
          best_mapping[component_number]
        ]
      )
    },
    FUN.VALUE = integer(1)
  )
  
  overall_agreement <- sum(
    matched_by_component
  ) / nrow(
    mapping_input
  )
  
  mapping_table <- data.frame(
    Internal_component = 1:3,
    Published_cluster = best_mapping,
    Matched_training_samples = matched_by_component,
    Total_training_samples = nrow(mapping_input),
    Overall_mapping_agreement = overall_agreement,
    stringsAsFactors = FALSE
  )
  
  list(
    mapping_table = mapping_table,
    contingency = contingency
  )
}


log_dirichlet_multinomial_kernel <- function(
    count_vector,
    alpha_vector,
    log_mixture_weight
) {
  
  total_count <- sum(
    count_vector
  )
  
  theta <- sum(
    alpha_vector
  )
  
  log_mixture_weight +
    lgamma(theta) -
    lgamma(
      total_count +
        theta
    ) +
    sum(
      lgamma(
        count_vector +
          alpha_vector
      ) -
        lgamma(
          alpha_vector
        )
    )
}


calculate_fixed_dmm_posterior <- function(
    count_matrix,
    alpha_matrix,
    mixture_weights
) {
  
  count_matrix <- as.matrix(
    count_matrix
  )
  
  alpha_matrix <- as.matrix(
    alpha_matrix
  )
  
  if (
    ncol(count_matrix) !=
    nrow(alpha_matrix)
  ) {
    stop(
      "Feature dimension mismatch between paired counts and model alpha.",
      call. = FALSE
    )
  }
  
  if (
    length(mixture_weights) !=
    ncol(alpha_matrix)
  ) {
    stop(
      "Mixture-weight count does not match model component count.",
      call. = FALSE
    )
  }
  
  if (
    any(!is.finite(alpha_matrix)) ||
    any(alpha_matrix <= 0)
  ) {
    stop(
      "Invalid alpha parameters in fixed DMM model.",
      call. = FALSE
    )
  }
  
  mixture_weights <- mixture_weights /
    sum(
      mixture_weights
    )
  
  log_weights <- log(
    mixture_weights
  )
  
  log_unnormalized <- matrix(
    NA_real_,
    nrow = nrow(count_matrix),
    ncol = ncol(alpha_matrix),
    dimnames = list(
      rownames(count_matrix),
      paste0(
        "Internal_component_",
        seq_len(ncol(alpha_matrix))
      )
    )
  )
  
  for (
    sample_index in seq_len(
      nrow(count_matrix)
    )
  ) {
    
    count_vector <- as.numeric(
      count_matrix[sample_index, ]
    )
    
    for (
      component_index in seq_len(
        ncol(alpha_matrix)
      )
    ) {
      
      alpha_vector <- as.numeric(
        alpha_matrix[, component_index]
      )
      
      log_unnormalized[
        sample_index,
        component_index
      ] <- log_dirichlet_multinomial_kernel(
        count_vector = count_vector,
        alpha_vector = alpha_vector,
        log_mixture_weight = log_weights[component_index]
      )
    }
  }
  
  posterior <- matrix(
    NA_real_,
    nrow = nrow(log_unnormalized),
    ncol = ncol(log_unnormalized),
    dimnames = dimnames(log_unnormalized)
  )
  
  for (
    sample_index in seq_len(
      nrow(log_unnormalized)
    )
  ) {
    
    current_logs <- log_unnormalized[
      sample_index,
      ,
      drop = TRUE
    ]
    
    max_log <- max(
      current_logs
    )
    
    stabilized <- exp(
      current_logs -
        max_log
    )
    
    posterior[
      sample_index,
    ] <- stabilized /
      sum(
        stabilized
      )
  }
  
  posterior
}


# ==============================================================================
# 4. Check files
# ==============================================================================

print_section(
  "FILE CHECK"
)

required_files <- c(
  model_file,
  primary_assignment_file,
  asv_count_file,
  taxonomy_file,
  metadata_file
)

file_check <- data.frame(
  File = required_files,
  Exists = file.exists(required_files),
  stringsAsFactors = FALSE
)

print(
  file_check,
  row.names = FALSE
)

if (any(!file_check$Exists)) {
  
  stop(
    paste0(
      "Missing required file(s):\n",
      paste(
        file_check$File[
          !file_check$Exists
        ],
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 5. Read fixed K3 model
# ==============================================================================

print_section(
  "READ FIXED K=3 MODEL"
)

k3_model <- readRDS(
  model_file
)

cat("Model class:\n")

print(
  class(k3_model)
)

if (!isS4(k3_model)) {
  stop(
    "Saved K=3 model is not an S4 object.",
    call. = FALSE
  )
}


model_slots <- methods::slotNames(
  k3_model
)

cat(
  "Model slots:\n"
)

print(
  model_slots
)


required_model_slots <- c(
  "fit",
  "mixture",
  "group"
)

missing_model_slots <- setdiff(
  required_model_slots,
  model_slots
)

if (length(missing_model_slots) > 0) {
  
  stop(
    paste0(
      "Model missing required slot(s): ",
      paste(
        missing_model_slots,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 6. Extract fixed alpha parameters
# ==============================================================================

print_section(
  "EXTRACT FIXED MODEL PARAMETERS"
)

fit_slot <- methods::slot(
  k3_model,
  "fit"
)

if (
  !is.list(fit_slot) ||
  !"Estimate" %in%
  names(fit_slot)
) {
  stop(
    "Model fit slot does not contain Estimate.",
    call. = FALSE
  )
}

alpha_matrix <- fit_slot[["Estimate"]]

alpha_matrix <- as.matrix(
  alpha_matrix
)

model_feature_names <- rownames(
  alpha_matrix
)

if (is.null(model_feature_names)) {
  stop(
    "Alpha matrix does not have genus row names.",
    call. = FALSE
  )
}

if (anyDuplicated(model_feature_names) > 0) {
  stop(
    "Duplicated genus names found in fixed K=3 model.",
    call. = FALSE
  )
}

number_model_features <- nrow(
  alpha_matrix
)

number_components <- ncol(
  alpha_matrix
)

if (number_components != 3) {
  stop(
    paste0(
      "Expected K=3 model but found ",
      number_components,
      " components."
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 7. Extract mixture weights
# ==============================================================================

mixture_slot <- methods::slot(
  k3_model,
  "mixture"
)

if (
  !is.list(mixture_slot) ||
  !"Weight" %in%
  names(mixture_slot)
) {
  stop(
    "Model mixture slot does not contain Weight.",
    call. = FALSE
  )
}

mixture_weights <- as.numeric(
  mixture_slot[["Weight"]]
)

if (length(mixture_weights) != 3) {
  stop(
    "Expected exactly three K=3 mixture weights.",
    call. = FALSE
  )
}

mixture_weights <- mixture_weights /
  sum(
    mixture_weights
  )

cat(
  "Fixed-model genera: ",
  number_model_features,
  "\n",
  sep = ""
)

cat(
  "Mixture weights:\n"
)

print(
  mixture_weights
)


# ==============================================================================
# 8. Recover model training internal components
# ==============================================================================

group_slot <- methods::slot(
  k3_model,
  "group"
)

group_matrix <- as.matrix(
  group_slot
)

training_sample_ids <- rownames(
  group_matrix
)

if (is.null(training_sample_ids)) {
  stop(
    "Model group slot does not have training sample row names.",
    call. = FALSE
  )
}

if (ncol(group_matrix) != 3) {
  stop(
    "Model group slot does not contain three components.",
    call. = FALSE
  )
}

training_internal_component <- max.col(
  group_matrix,
  ties.method = "first"
)

names(
  training_internal_component
) <- training_sample_ids

cat(
  "Training samples in model: ",
  length(training_sample_ids),
  "\n",
  sep = ""
)


# ==============================================================================
# 9. Read primary assignment table
# ==============================================================================

primary_assignments <- data.table::fread(
  primary_assignment_file,
  data.table = FALSE,
  check.names = FALSE
)

if (
  !all(
    c(
      "SampleID",
      "DMM_cluster"
    ) %in%
    colnames(primary_assignments)
  )
) {
  stop(
    "Primary assignment table must contain SampleID and DMM_cluster.",
    call. = FALSE
  )
}

primary_assignments$SampleID <- as.character(
  primary_assignments$SampleID
)

primary_assignments$DMM_cluster <- standardize_cluster_label(
  primary_assignments$DMM_cluster
)


# ==============================================================================
# 10. Recover internal component -> C1/C2/C3 mapping
# ==============================================================================

print_section(
  "RECOVER INTERNAL COMPONENT MAPPING"
)

primary_index <- match(
  training_sample_ids,
  primary_assignments$SampleID
)

if (anyNA(primary_index)) {
  stop(
    "At least one model training sample is absent from the primary assignment table.",
    call. = FALSE
  )
}

published_training_cluster <- primary_assignments$DMM_cluster[
  primary_index
]

mapping_result <- derive_component_mapping(
  internal_component = training_internal_component,
  published_cluster = published_training_cluster
)

component_mapping <- mapping_result[["mapping_table"]]
mapping_contingency <- mapping_result[["contingency"]]

print(
  component_mapping,
  row.names = FALSE
)

cat(
  "\nMapping contingency:\n"
)

print(
  mapping_contingency
)

mapping_agreement <- unique(
  component_mapping$Overall_mapping_agreement
)

if (
  length(mapping_agreement) != 1 ||
  is.na(mapping_agreement) ||
  mapping_agreement < 0.95
) {
  stop(
    paste0(
      "Internal component mapping agreement is only ",
      signif(
        mapping_agreement,
        4
      ),
      "."
    ),
    call. = FALSE
  )
}

data.table::fwrite(
  component_mapping,
  component_mapping_file,
  sep = "\t",
  quote = FALSE
)

mapping_contingency_df <- as.data.frame.matrix(
  mapping_contingency
)

mapping_contingency_df$Internal_component <- rownames(
  mapping_contingency_df
)

mapping_contingency_df <- mapping_contingency_df[
  ,
  c(
    "Internal_component",
    "C1",
    "C2",
    "C3"
  ),
  drop = FALSE
]

data.table::fwrite(
  mapping_contingency_df,
  component_mapping_contingency_file,
  sep = "\t",
  quote = FALSE
)

cluster_to_internal <- stats::setNames(
  component_mapping$Internal_component,
  component_mapping$Published_cluster
)


# ==============================================================================
# 11. Read paired dataset
# ==============================================================================

print_section(
  "READ FIXED PAIRED 23 + 23 DATA"
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


required_metadata_columns <- c(
  "PairID",
  "SampleID",
  "CA_vs_nonCA_explicit"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    paste0(
      "Paired metadata missing: ",
      paste(
        missing_metadata_columns,
        collapse = ", "
      )
    ),
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

if (!"ASV" %in% colnames(asv_count)) {
  stop(
    "ASV count table must contain ASV.",
    call. = FALSE
  )
}


# ==============================================================================
# 12. Standardize paired sample metadata
# ==============================================================================

sample_ids <- as.character(
  metadata$SampleID
)

sample_groups <- as.character(
  metadata$CA_vs_nonCA_explicit
)

sample_groups[
  sample_groups %in% c(
    "nonCAC",
    "nonCA"
  )
] <- "nonCA"

sample_groups[
  sample_groups %in% c(
    "CAC",
    "CA"
  )
] <- "CA"


if (anyDuplicated(sample_ids) > 0) {
  stop(
    "Duplicated SampleID values detected in paired metadata.",
    call. = FALSE
  )
}

if (
  !all(
    sample_groups %in%
    c(
      "nonCA",
      "CA"
    )
  )
) {
  stop(
    paste0(
      "Unexpected CA_vs_nonCA_explicit values: ",
      paste(
        unique(sample_groups),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

if (sum(sample_groups == "nonCA") != 23) {
  stop(
    "Expected exactly 23 nonCAC samples.",
    call. = FALSE
  )
}

if (sum(sample_groups == "CA") != 23) {
  stop(
    "Expected exactly 23 CAC samples.",
    call. = FALSE
  )
}


# ==============================================================================
# 13. Strict PairID audit
# ==============================================================================

pair_metadata <- data.frame(
  PairID = as.character(metadata$PairID),
  SampleID = sample_ids,
  Group = sample_groups,
  stringsAsFactors = FALSE
)

pair_counts <- table(
  pair_metadata$PairID
)

if (
  length(pair_counts) != 23 ||
  any(pair_counts != 2)
) {
  stop(
    "Expected exactly 23 PairIDs, each containing two samples.",
    call. = FALSE
  )
}

pair_group_table <- table(
  pair_metadata$PairID,
  pair_metadata$Group
)

if (
  !"CA" %in%
  colnames(pair_group_table) ||
  !"nonCA" %in%
  colnames(pair_group_table)
) {
  stop(
    "Pair audit failed because CA/nonCA columns are incomplete.",
    call. = FALSE
  )
}

if (
  any(
    pair_group_table[, "CA"] != 1
  ) ||
  any(
    pair_group_table[, "nonCA"] != 1
  )
) {
  stop(
    "Each PairID must contain exactly one nonCAC and one CAC sample.",
    call. = FALSE
  )
}


# ==============================================================================
# 14. Check paired samples exist in ASV table
# ==============================================================================

missing_samples <- setdiff(
  sample_ids,
  colnames(asv_count)
)

if (length(missing_samples) > 0) {
  stop(
    paste0(
      "Paired samples absent from ASV count table: ",
      paste(
        missing_samples,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 15. Align taxonomy to ASV count table
# ==============================================================================

print_section(
  "ALIGN ASV TAXONOMY"
)

if (anyDuplicated(asv_count$ASV) > 0) {
  stop(
    "Duplicated ASV IDs in paired count table.",
    call. = FALSE
  )
}

if (anyDuplicated(taxonomy$ASV) > 0) {
  stop(
    "Duplicated ASV IDs in paired taxonomy.",
    call. = FALSE
  )
}

taxonomy_match <- match(
  as.character(asv_count$ASV),
  as.character(taxonomy$ASV)
)

if (anyNA(taxonomy_match)) {
  stop(
    paste0(
      sum(is.na(taxonomy_match)),
      " ASVs lack taxonomy."
    ),
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

genus_vector[
  is_bad_genus(genus_vector)
] <- NA_character_

cat(
  "ASVs total: ",
  length(genus_vector),
  "\n",
  sep = ""
)

cat(
  "ASVs with valid strict genus: ",
  sum(!is.na(genus_vector)),
  "\n",
  sep = ""
)


# ==============================================================================
# 16. Build raw ASV count matrix
# ==============================================================================

print_section(
  "BUILD RAW INTEGER ASV MATRIX"
)

asv_count_matrix <- as.matrix(
  asv_count[
    ,
    sample_ids,
    drop = FALSE
  ]
)

storage.mode(
  asv_count_matrix
) <- "numeric"

rownames(
  asv_count_matrix
) <- as.character(
  asv_count$ASV
)

if (any(!is.finite(asv_count_matrix))) {
  stop(
    "ASV count matrix contains non-finite values.",
    call. = FALSE
  )
}

if (any(asv_count_matrix < 0)) {
  stop(
    "Negative ASV counts detected.",
    call. = FALSE
  )
}

if (
  any(
    abs(
      asv_count_matrix -
      round(asv_count_matrix)
    ) > 1e-8
  )
) {
  stop(
    "Non-integer ASV counts detected.",
    call. = FALSE
  )
}

asv_count_matrix <- round(
  asv_count_matrix
)

original_sample_depth <- colSums(
  asv_count_matrix
)

if (any(original_sample_depth <= 0)) {
  stop(
    "At least one paired sample has zero ASV count.",
    call. = FALSE
  )
}


# ==============================================================================
# 17. Collapse paired ASVs to strict genus
# ==============================================================================

print_section(
  "COLLAPSE ASVS TO STRICT GENUS"
)

valid_genus_rows <- !is.na(
  genus_vector
)

if (sum(valid_genus_rows) == 0) {
  stop(
    "No valid genus-level ASVs are available.",
    call. = FALSE
  )
}

genus_count_matrix <- rowsum(
  asv_count_matrix[
    valid_genus_rows,
    ,
    drop = FALSE
  ],
  group = genus_vector[
    valid_genus_rows
  ],
  reorder = FALSE
)

cat(
  "Paired-data genera after strict collapse: ",
  nrow(genus_count_matrix),
  "\n",
  sep = ""
)


# ==============================================================================
# 18. Align paired genus counts to fixed-model feature space
# ==============================================================================

print_section(
  "ALIGN TO FIXED MODEL FEATURE SPACE"
)

shared_features <- intersect(
  model_feature_names,
  rownames(genus_count_matrix)
)

model_only_features <- setdiff(
  model_feature_names,
  rownames(genus_count_matrix)
)

paired_only_features <- setdiff(
  rownames(genus_count_matrix),
  model_feature_names
)

aligned_genus_matrix <- matrix(
  0,
  nrow = number_model_features,
  ncol = length(sample_ids),
  dimnames = list(
    model_feature_names,
    sample_ids
  )
)

aligned_genus_matrix[
  shared_features,
  sample_ids
] <- genus_count_matrix[
  shared_features,
  sample_ids,
  drop = FALSE
]

feature_alignment_summary <- data.frame(
  Metric = c(
    "Fixed_model_features",
    "Shared_features",
    "Model_features_absent_in_paired_data",
    "Paired_data_genera_outside_fixed_model",
    "Paired_data_genera_total"
  ),
  Value = c(
    number_model_features,
    length(shared_features),
    length(model_only_features),
    length(paired_only_features),
    nrow(genus_count_matrix)
  ),
  stringsAsFactors = FALSE
)

print(
  feature_alignment_summary,
  row.names = FALSE
)

data.table::fwrite(
  feature_alignment_summary,
  feature_alignment_file,
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 19. Count-retention QC
# ==============================================================================

print_section(
  "COUNT RETENTION QC"
)

valid_genus_counts <- colSums(
  genus_count_matrix
)

model_feature_counts <- colSums(
  aligned_genus_matrix
)

sample_retention <- data.frame(
  PairID = as.character(metadata$PairID),
  SampleID = sample_ids,
  Group = sample_groups,
  Original_ASV_count = as.numeric(
    original_sample_depth[sample_ids]
  ),
  Count_assigned_to_valid_genus = as.numeric(
    valid_genus_counts[sample_ids]
  ),
  Count_in_fixed_model_features = as.numeric(
    model_feature_counts[sample_ids]
  ),
  stringsAsFactors = FALSE
)

sample_retention$Valid_genus_fraction <-
  sample_retention$Count_assigned_to_valid_genus /
  sample_retention$Original_ASV_count

sample_retention$Fixed_model_feature_fraction <-
  sample_retention$Count_in_fixed_model_features /
  sample_retention$Original_ASV_count

sample_retention$Fixed_model_fraction_of_valid_genus <-
  sample_retention$Count_in_fixed_model_features /
  sample_retention$Count_assigned_to_valid_genus


if (
  any(
    !is.finite(
      sample_retention[["Fixed_model_fraction_of_valid_genus"]]
    )
  )
) {
  stop(
    "At least one sample has zero valid-genus counts.",
    call. = FALSE
  )
}

if (
  any(
    sample_retention[["Count_in_fixed_model_features"]] <= 0
  )
) {
  stop(
    "At least one sample has zero counts in the fixed model feature space.",
    call. = FALSE
  )
}

data.table::fwrite(
  sample_retention,
  sample_retention_file,
  sep = "\t",
  quote = FALSE
)

cat(
  "Fixed-model feature fraction summary:\n"
)

print(
  summary(
    sample_retention[["Fixed_model_feature_fraction"]]
  )
)


# ==============================================================================
# 20. Build sample x model-genus projection matrix
# ==============================================================================

prediction_count_matrix <- t(
  aligned_genus_matrix
)

prediction_count_matrix <- prediction_count_matrix[
  sample_ids,
  model_feature_names,
  drop = FALSE
]

storage.mode(
  prediction_count_matrix
) <- "numeric"

if (
  !identical(
    rownames(prediction_count_matrix),
    sample_ids
  )
) {
  stop(
    "Projection sample order is incorrect.",
    call. = FALSE
  )
}

if (
  !identical(
    colnames(prediction_count_matrix),
    model_feature_names
  )
) {
  stop(
    "Projection genus feature order is incorrect.",
    call. = FALSE
  )
}

cat(
  "Projection matrix dimensions: ",
  nrow(prediction_count_matrix),
  " samples x ",
  ncol(prediction_count_matrix),
  " genera\n",
  sep = ""
)


# ==============================================================================
# 21. Calculate fixed-model posterior probabilities
# ==============================================================================

print_section(
  "CALCULATE FIXED K=3 POSTERIORS"
)

posterior_internal <- calculate_fixed_dmm_posterior(
  count_matrix = prediction_count_matrix,
  alpha_matrix = alpha_matrix,
  mixture_weights = mixture_weights
)

if (
  any(
    abs(
      rowSums(posterior_internal) -
      1
    ) > 1e-8
  )
) {
  stop(
    "Fixed-model posterior probabilities do not sum to one.",
    call. = FALSE
  )
}


# ==============================================================================
# 22. Reorder internal posteriors into C1/C2/C3
# ==============================================================================

posterior_cluster <- matrix(
  NA_real_,
  nrow = length(sample_ids),
  ncol = 3,
  dimnames = list(
    sample_ids,
    c(
      "Posterior_C1",
      "Posterior_C2",
      "Posterior_C3"
    )
  )
)

for (
  cluster_name in
  c(
    "C1",
    "C2",
    "C3"
  )
) {
  
  internal_component <- cluster_to_internal[[cluster_name]]
  
  posterior_cluster[
    ,
    paste0(
      "Posterior_",
      cluster_name
    )
  ] <- posterior_internal[
    ,
    paste0(
      "Internal_component_",
      internal_component
    )
  ]
}


# ==============================================================================
# 23. Hard assignments + confidence
# ==============================================================================

assigned_index <- max.col(
  posterior_cluster,
  ties.method = "first"
)

assigned_cluster <- c(
  "C1",
  "C2",
  "C3"
)[assigned_index]

maximum_posterior <- apply(
  posterior_cluster,
  1,
  max
)

sorted_posterior <- t(
  apply(
    posterior_cluster,
    1,
    sort,
    decreasing = TRUE
  )
)

posterior_margin <-
  sorted_posterior[, 1] -
  sorted_posterior[, 2]

posterior_entropy <- apply(
  posterior_cluster,
  1,
  function(probability_vector) {
    
    nonzero_probability <- probability_vector[
      probability_vector > 0
    ]
    
    -sum(
      nonzero_probability *
        log(
          nonzero_probability
        )
    )
  }
)

normalized_entropy <-
  posterior_entropy /
  log(3)

confidence_category <- cut(
  maximum_posterior,
  breaks = c(
    -Inf,
    0.60,
    0.80,
    0.95,
    Inf
  ),
  labels = c(
    "Low_<0.60",
    "Moderate_0.60-0.79",
    "High_0.80-0.94",
    "Very_high_>=0.95"
  ),
  right = FALSE
)


# ==============================================================================
# 24. Build full projection table
# ==============================================================================

print_section(
  "BUILD PROJECTION TABLE"
)

projection_all <- data.frame(
  PairID = as.character(metadata$PairID),
  SampleID = sample_ids,
  Group = sample_groups,
  Assigned_cluster = assigned_cluster,
  Maximum_posterior = maximum_posterior,
  Posterior_margin = posterior_margin,
  Posterior_entropy = posterior_entropy,
  Normalized_entropy = normalized_entropy,
  Confidence_category = as.character(confidence_category),
  Posterior_C1 = posterior_cluster[, "Posterior_C1"],
  Posterior_C2 = posterior_cluster[, "Posterior_C2"],
  Posterior_C3 = posterior_cluster[, "Posterior_C3"],
  Original_ASV_count = sample_retention$Original_ASV_count,
  Count_in_fixed_model_features =
    sample_retention$Count_in_fixed_model_features,
  Fixed_model_feature_fraction =
    sample_retention$Fixed_model_feature_fraction,
  stringsAsFactors = FALSE
)

projection_all <- projection_all[
  order(
    projection_all$PairID,
    factor(
      projection_all$Group,
      levels = c(
        "nonCA",
        "CA"
      )
    )
  ),
  ,
  drop = FALSE
]

projection_nonca <- projection_all[
  projection_all$Group == "nonCA",
  ,
  drop = FALSE
]

projection_ca <- projection_all[
  projection_all$Group == "CA",
  ,
  drop = FALSE
]


if (nrow(projection_nonca) != 23) {
  stop(
    "Projected nonCAC sample count is not 23.",
    call. = FALSE
  )
}

if (nrow(projection_ca) != 23) {
  stop(
    "Projected CAC sample count is not 23.",
    call. = FALSE
  )
}


data.table::fwrite(
  projection_all,
  projection_all_file,
  sep = "\t",
  quote = FALSE
)

data.table::fwrite(
  projection_all,
  projection_compatibility_file,
  sep = "\t",
  quote = FALSE
)

data.table::fwrite(
  projection_nonca,
  projection_nonca_file,
  sep = "\t",
  quote = FALSE
)

data.table::fwrite(
  projection_ca,
  projection_ca_file,
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 25. Projected state distribution by group
# ==============================================================================

print_section(
  "PROJECTED STATE DISTRIBUTION"
)

group_levels <- c(
  "nonCA",
  "CA"
)

cluster_levels <- c(
  "C1",
  "C2",
  "C3"
)

cluster_table <- table(
  factor(
    projection_all$Group,
    levels = group_levels
  ),
  factor(
    projection_all$Assigned_cluster,
    levels = cluster_levels
  )
)

cluster_distribution <- as.data.frame(
  cluster_table,
  stringsAsFactors = FALSE
)

colnames(
  cluster_distribution
) <- c(
  "Group",
  "Assigned_cluster",
  "N"
)

cluster_distribution$Group_total <- ave(
  cluster_distribution$N,
  cluster_distribution$Group,
  FUN = sum
)

cluster_distribution$Proportion <-
  cluster_distribution$N /
  cluster_distribution$Group_total

cluster_distribution$Percent <-
  100 *
  cluster_distribution$Proportion

data.table::fwrite(
  cluster_distribution,
  cluster_distribution_file,
  sep = "\t",
  quote = FALSE
)

print(
  cluster_distribution,
  row.names = FALSE
)


# ==============================================================================
# 26. Assignment confidence summary
# ==============================================================================

print_section(
  "ASSIGNMENT CONFIDENCE"
)

confidence_summary_list <- lapply(
  group_levels,
  function(group_name) {
    
    current_data <- projection_all[
      projection_all$Group == group_name,
      ,
      drop = FALSE
    ]
    
    data.frame(
      Group = group_name,
      N = nrow(current_data),
      Median_maximum_posterior = median(
        current_data$Maximum_posterior
      ),
      IQR_maximum_posterior = IQR(
        current_data$Maximum_posterior
      ),
      Minimum_maximum_posterior = min(
        current_data$Maximum_posterior
      ),
      N_below_0_60 = sum(
        current_data$Maximum_posterior < 0.60
      ),
      N_below_0_80 = sum(
        current_data$Maximum_posterior < 0.80
      ),
      Median_posterior_margin = median(
        current_data$Posterior_margin
      ),
      stringsAsFactors = FALSE
    )
  }
)

confidence_summary <- do.call(
  rbind,
  confidence_summary_list
)

data.table::fwrite(
  confidence_summary,
  confidence_summary_file,
  sep = "\t",
  quote = FALSE
)

print(
  confidence_summary,
  row.names = FALSE
)


# ==============================================================================
# 27. Validate projected CAC against primary K3 assignments
# ==============================================================================

print_section(
  "CAC VALIDATION AGAINST PRIMARY ASSIGNMENT"
)

ca_primary_index <- match(
  projection_ca$SampleID,
  primary_assignments$SampleID
)

ca_validation <- data.frame(
  PairID = projection_ca$PairID,
  SampleID = projection_ca$SampleID,
  Projected_cluster = projection_ca$Assigned_cluster,
  Maximum_posterior = projection_ca$Maximum_posterior,
  Posterior_C1 = projection_ca$Posterior_C1,
  Posterior_C2 = projection_ca$Posterior_C2,
  Posterior_C3 = projection_ca$Posterior_C3,
  Present_in_primary_model = !is.na(ca_primary_index),
  stringsAsFactors = FALSE
)

ca_validation$Primary_cluster <- NA_character_

valid_ca_primary <- !is.na(
  ca_primary_index
)

ca_validation$Primary_cluster[
  valid_ca_primary
] <- primary_assignments$DMM_cluster[
  ca_primary_index[
    valid_ca_primary
  ]
]

ca_validation$Exact_agreement <- NA

comparable_ca <-
  ca_validation$Present_in_primary_model &
  !is.na(
    ca_validation$Primary_cluster
  )

ca_validation$Exact_agreement[
  comparable_ca
] <-
  ca_validation$Projected_cluster[
    comparable_ca
  ] ==
  ca_validation$Primary_cluster[
    comparable_ca
  ]

number_comparable_ca <- sum(
  comparable_ca
)

exact_matches <- sum(
  ca_validation$Exact_agreement,
  na.rm = TRUE
)

exact_agreement_rate <- if (
  number_comparable_ca > 0
) {
  exact_matches /
    number_comparable_ca
} else {
  NA_real_
}

ca_validation_summary <- data.frame(
  CA_samples_projected = nrow(
    ca_validation
  ),
  CA_samples_found_in_primary_assignment = sum(
    ca_validation$Present_in_primary_model
  ),
  CA_samples_comparable = number_comparable_ca,
  Exact_matches = exact_matches,
  Exact_agreement_rate = exact_agreement_rate,
  stringsAsFactors = FALSE
)

data.table::fwrite(
  ca_validation,
  ca_validation_file,
  sep = "\t",
  quote = FALSE
)

data.table::fwrite(
  ca_validation_summary,
  ca_validation_summary_file,
  sep = "\t",
  quote = FALSE
)

print(
  ca_validation_summary,
  row.names = FALSE
)


# ==============================================================================
# 28. Build paired nonCAC -> CAC transition table
# ==============================================================================

print_section(
  "BUILD PAIRED TRANSITION TABLE"
)

nonca_pair <- projection_nonca[
  ,
  c(
    "PairID",
    "SampleID",
    "Assigned_cluster",
    "Maximum_posterior"
  ),
  drop = FALSE
]

colnames(
  nonca_pair
) <- c(
  "PairID",
  "nonCA_SampleID",
  "nonCA_state",
  "nonCA_maximum_posterior"
)

ca_pair <- projection_ca[
  ,
  c(
    "PairID",
    "SampleID",
    "Assigned_cluster",
    "Maximum_posterior"
  ),
  drop = FALSE
]

colnames(
  ca_pair
) <- c(
  "PairID",
  "CA_SampleID",
  "CA_state",
  "CA_maximum_posterior"
)

pair_transition <- merge(
  nonca_pair,
  ca_pair,
  by = "PairID",
  all = FALSE,
  sort = TRUE
)

if (nrow(pair_transition) != 23) {
  stop(
    paste0(
      "Expected 23 paired transitions but found ",
      nrow(pair_transition),
      "."
    ),
    call. = FALSE
  )
}

data.table::fwrite(
  pair_transition,
  pair_transition_file,
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 29. Transition matrix
# ==============================================================================

transition_matrix <- table(
  factor(
    pair_transition$nonCA_state,
    levels = cluster_levels
  ),
  factor(
    pair_transition$CA_state,
    levels = cluster_levels
  )
)

cat(
  "Rows = paired nonCAC state; columns = paired CAC state:\n"
)

print(
  transition_matrix
)

transition_matrix_df <- as.data.frame.matrix(
  transition_matrix
)

transition_matrix_df$nonCA_state <- rownames(
  transition_matrix_df
)

transition_matrix_df <- transition_matrix_df[
  ,
  c(
    "nonCA_state",
    "C1",
    "C2",
    "C3"
  ),
  drop = FALSE
]

data.table::fwrite(
  transition_matrix_df,
  pair_transition_matrix_file,
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 30. Final summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("7KB fixed-K3 paired projection completed successfully.\n")
cat("============================================================\n")

cat("\n")
cat(
  "Fixed-model genus features: ",
  number_model_features,
  "\n",
  sep = ""
)

cat("\nFeature alignment:\n")
print(
  feature_alignment_summary,
  row.names = FALSE
)

cat("\nInternal component mapping:\n")
print(
  component_mapping,
  row.names = FALSE
)

cat("\nProjected state distribution:\n")
print(
  cluster_distribution,
  row.names = FALSE
)

cat("\nAssignment confidence:\n")
print(
  confidence_summary,
  row.names = FALSE
)

cat("\nCAC validation against primary assignment:\n")
print(
  ca_validation_summary,
  row.names = FALSE
)

cat("\nPaired transition matrix:\n")
print(
  transition_matrix
)

cat("\nProjection compatibility file:\n")
cat(
  projection_compatibility_file,
  "\n"
)

cat("\nOutput directory:\n")
cat(
  output_dir,
  "\n"
)