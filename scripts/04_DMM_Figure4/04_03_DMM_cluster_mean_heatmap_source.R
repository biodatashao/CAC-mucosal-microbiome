#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_03_DMM_cluster_mean_heatmap_source.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Purpose:
## Generate the 7KB K=3 DMM state-defining genus selection and
## row-normalized cluster-mean heatmap source.
##
## Original selection rule:
## 1. Use all genera in frozen primary DMM prevalence10 input matrix.
## 2. Calculate mean relative abundance in C1, C2, C3.
## 3. Assign each genus to the cluster with highest mean abundance.
## 4. Leading margin = highest mean - second-highest mean.
## 5. Require:
##       Highest mean >= 0.5%
##       Leading margin >= 0.5 percentage points
## 6. Rank qualified genera within assigned cluster by leading margin.
## 7. Retain at most 13 genera per cluster.
## 8. Do NOT fill a cluster to 13 if fewer qualify.
## 9. Highest-to-second ratio is audit only, NOT a filter.
## 10. Row-normalize C1/C2/C3 cluster means.
##
## IMPORTANT:
## This is the upstream source-generation step.
## Final Figure 4C curation is performed later.
############################################################


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "readr",
  "tibble",
  "stringr",
  "ggplot2",
  "scales"
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
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(scales)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

dmm_primary_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "DMM_progression127_primary_prevalence10_7KB"
)

clean_data_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "00_clean_data",
  "progression127"
)

genus_count_file <- file.path(
  dmm_primary_dir,
  "tables",
  "genus_raw_count_progression127_prevalence10_frozen_7KB.tsv"
)

assignment_file <- file.path(
  dmm_primary_dir,
  "tables",
  "DMM_sample_assignments_optimalK3_native_7KB.tsv"
)

taxonomy_file <- file.path(
  clean_data_dir,
  "taxonomy_7KB_progression127.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "Figure4C_Lavelle_cluster_mean_heatmap_7KB"
)

table_dir <- file.path(
  output_dir,
  "source_tables"
)

figure_dir <- file.path(
  output_dir,
  "figures"
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Original parameters
# ==============================================================================

cluster_order <- c(
  "C1",
  "C2",
  "C3"
)

maximum_genera_per_cluster <- 13L

# 0.005 proportion = 0.5%
minimum_highest_mean <- 0.005

# 0.005 proportion = 0.5 percentage points
minimum_leading_margin <- 0.005

# Original script clips row z-scores at +/-1.5
row_z_limit <- 1.5


# ==============================================================================
# 3. Helper functions
# ==============================================================================

stop_if_missing <- function(paths) {
  
  missing_paths <- paths[
    !file.exists(paths)
  ]
  
  if (length(missing_paths) > 0) {
    stop(
      paste(
        "The following required files do not exist:",
        paste(missing_paths, collapse = "\n"),
        sep = "\n"
      ),
      call. = FALSE
    )
  }
}


# Same strict genus parsing principle already used in the 7KB Figure 3/DMM run.
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
    
    candidate[
      looks_higher_tax
    ] <- NA_character_
    
    candidate[
      candidate == ""
    ] <- NA_character_
    
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
  
  x <- stringr::str_trim(x)
  
  x <- stringr::str_replace(
    x,
    "^[A-Za-z]__",
    ""
  )
  
  x <- stringr::str_replace_all(
    x,
    "\\[|\\]",
    ""
  )
  
  x <- stringr::str_replace_all(
    x,
    "[^A-Za-z0-9]+",
    "_"
  )
  
  x <- stringr::str_replace_all(
    x,
    "_+",
    "_"
  )
  
  x <- stringr::str_replace_all(
    x,
    "^_|_$",
    ""
  )
  
  stringr::str_to_lower(x)
}


extract_phylum <- function(taxonomy_string) {
  
  taxonomy_string <- as.character(
    taxonomy_string
  )
  
  phylum_value <- stringr::str_extract(
    taxonomy_string,
    "(?<=p__)[^;]+"
  )
  
  phylum_value[
    is.na(phylum_value) |
      phylum_value == ""
  ] <- "Unclassified"
  
  phylum_value
}


standardize_phylum_group <- function(phylum_name) {
  
  dplyr::case_when(
    
    phylum_name %in% c(
      "Bacillota",
      "Firmicutes"
    ) ~ "Firmicutes/Bacillota",
    
    phylum_name %in% c(
      "Pseudomonadota",
      "Proteobacteria"
    ) ~ "Proteobacteria/Pseudomonadota",
    
    phylum_name %in% c(
      "Bacteroidota",
      "Bacteroidetes"
    ) ~ "Bacteroidetes/Bacteroidota",
    
    phylum_name %in% c(
      "Actinomycetota",
      "Actinobacteria"
    ) ~ "Actinobacteria/Actinomycetota",
    
    phylum_name %in% c(
      "Verrucomicrobiota",
      "Verrucomicrobia"
    ) ~ "Verrucomicrobiota",
    
    phylum_name %in% c(
      "Fusobacteriota",
      "Fusobacteria"
    ) ~ "Fusobacteriota",
    
    TRUE ~ "Other"
  )
}


row_normalize_three_values <- function(x) {
  
  x <- as.numeric(x)
  
  current_sd <- stats::sd(
    x,
    na.rm = TRUE
  )
  
  if (
    is.na(current_sd) ||
    current_sd == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }
  
  normalized_values <- (
    x -
      mean(
        x,
        na.rm = TRUE
      )
  ) / current_sd
  
  normalized_values[
    normalized_values >
      row_z_limit
  ] <- row_z_limit
  
  normalized_values[
    normalized_values <
      -row_z_limit
  ] <- -row_z_limit
  
  normalized_values
}


# ==============================================================================
# 4. Validate inputs
# ==============================================================================

stop_if_missing(
  c(
    genus_count_file,
    assignment_file,
    taxonomy_file
  )
)


# ==============================================================================
# 5. Read inputs
# ==============================================================================

genus_count_raw <- data.table::fread(
  genus_count_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

assignment_raw <- data.table::fread(
  assignment_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

taxonomy_raw <- data.table::fread(
  taxonomy_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)


# ==============================================================================
# 6. Validate columns
# ==============================================================================

if (
  !"SampleID" %in%
  colnames(genus_count_raw)
) {
  stop(
    "Genus count matrix does not contain SampleID.",
    call. = FALSE
  )
}

required_assignment_columns <- c(
  "SampleID",
  "DMM_cluster"
)

missing_assignment_columns <- setdiff(
  required_assignment_columns,
  colnames(assignment_raw)
)

if (
  length(
    missing_assignment_columns
  ) > 0
) {
  stop(
    "Missing assignment columns: ",
    paste(
      missing_assignment_columns,
      collapse = ", "
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
    colnames(
      taxonomy_raw
    )
  )
) {
  stop(
    "7KB taxonomy file must contain ASV and Taxonomy.",
    call. = FALSE
  )
}


# ==============================================================================
# 7. Build genus count and relative-abundance matrices
# ==============================================================================

sample_ids <- as.character(
  genus_count_raw$SampleID
)

if (
  anyDuplicated(
    sample_ids
  ) > 0
) {
  stop(
    "Duplicated sample IDs detected in genus count matrix.",
    call. = FALSE
  )
}

genus_count_matrix <- genus_count_raw %>%
  dplyr::select(
    -SampleID
  ) %>%
  as.matrix()

storage.mode(
  genus_count_matrix
) <- "numeric"

rownames(
  genus_count_matrix
) <- sample_ids

if (
  anyNA(
    genus_count_matrix
  )
) {
  stop(
    "NA values detected in genus count matrix.",
    call. = FALSE
  )
}

if (
  any(
    genus_count_matrix < 0
  )
) {
  stop(
    "Negative genus counts detected.",
    call. = FALSE
  )
}

if (
  any(
    rowSums(
      genus_count_matrix
    ) <= 0
  )
) {
  stop(
    "At least one sample has zero genus-level counts.",
    call. = FALSE
  )
}

relative_abundance_matrix <- sweep(
  genus_count_matrix,
  1,
  rowSums(
    genus_count_matrix
  ),
  "/"
)


# ==============================================================================
# 8. Prepare K=3 assignment metadata
# ==============================================================================

sample_metadata <- assignment_raw %>%
  transmute(
    SampleID = as.character(
      SampleID
    ),
    
    DMM_cluster = factor(
      as.character(
        DMM_cluster
      ),
      levels = cluster_order
    )
  )

if (
  anyDuplicated(
    sample_metadata$SampleID
  ) > 0
) {
  stop(
    "Duplicated sample IDs detected in assignment table.",
    call. = FALSE
  )
}

if (
  anyNA(
    sample_metadata$DMM_cluster
  )
) {
  stop(
    "Unexpected or missing DMM state labels.",
    call. = FALSE
  )
}

if (
  !setequal(
    sample_metadata$SampleID,
    rownames(
      relative_abundance_matrix
    )
  )
) {
  stop(
    paste0(
      "Sample IDs do not match between assignment table ",
      "and frozen genus count matrix."
    ),
    call. = FALSE
  )
}

# Strict audit from the completed 7KB DMM run
state_counts <- sample_metadata %>%
  count(
    DMM_cluster,
    name = "N"
  )

print(
  state_counts
)

expected_state_counts <- c(
  C1 = 55,
  C2 = 46,
  C3 = 26
)

observed_state_counts <- setNames(
  state_counts$N,
  as.character(
    state_counts$DMM_cluster
  )
)

if (
  !identical(
    as.integer(
      observed_state_counts[
        names(
          expected_state_counts
        )
      ]
    ),
    as.integer(
      expected_state_counts
    )
  )
) {
  stop(
    "Expected state counts C1=55, C2=46, C3=26 were not recovered.",
    call. = FALSE
  )
}


# ==============================================================================
# 9. Calculate cluster-level summary statistics
# ==============================================================================

cluster_summary_long <- dplyr::bind_rows(
  lapply(
    cluster_order,
    function(current_cluster) {
      
      current_samples <- sample_metadata %>%
        filter(
          DMM_cluster ==
            current_cluster
        ) %>%
        pull(
          SampleID
        )
      
      current_matrix <- relative_abundance_matrix[
        current_samples,
        ,
        drop = FALSE
      ]
      
      tibble(
        Genus = colnames(
          current_matrix
        ),
        
        DMM_cluster =
          current_cluster,
        
        Cluster_sample_n =
          length(
            current_samples
          ),
        
        Mean_relative_abundance =
          colMeans(
            current_matrix,
            na.rm = TRUE
          ),
        
        Median_relative_abundance =
          apply(
            current_matrix,
            2,
            stats::median,
            na.rm = TRUE
          ),
        
        Prevalence_fraction =
          colMeans(
            current_matrix > 0,
            na.rm = TRUE
          )
      )
    }
  )
)

write_tsv(
  cluster_summary_long,
  file.path(
    table_dir,
    "all_genera_cluster_summary_long_7KB.tsv"
  )
)


# ==============================================================================
# 10. Convert cluster means to wide form
# ==============================================================================

cluster_mean_wide <- cluster_summary_long %>%
  select(
    Genus,
    DMM_cluster,
    Mean_relative_abundance
  ) %>%
  pivot_wider(
    names_from = DMM_cluster,
    values_from = Mean_relative_abundance,
    names_prefix = "Mean_"
  )

required_mean_columns <- c(
  "Mean_C1",
  "Mean_C2",
  "Mean_C3"
)

missing_mean_columns <- setdiff(
  required_mean_columns,
  colnames(
    cluster_mean_wide
  )
)

if (
  length(
    missing_mean_columns
  ) > 0
) {
  stop(
    "Missing cluster mean columns: ",
    paste(
      missing_mean_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 11. Assign genus to its highest-mean state
# ==============================================================================

mean_value_matrix <- cluster_mean_wide %>%
  select(
    Mean_C1,
    Mean_C2,
    Mean_C3
  ) %>%
  as.matrix()

highest_mean_values <- apply(
  mean_value_matrix,
  1,
  max,
  na.rm = TRUE
)

second_highest_mean_values <- apply(
  mean_value_matrix,
  1,
  function(x) {
    
    sorted_values <- sort(
      as.numeric(x),
      decreasing = TRUE,
      na.last = TRUE
    )
    
    sorted_values[2]
  }
)

lowest_mean_values <- apply(
  mean_value_matrix,
  1,
  min,
  na.rm = TRUE
)

highest_cluster_index <- max.col(
  mean_value_matrix,
  ties.method = "first"
)

assigned_cluster_values <- cluster_order[
  highest_cluster_index
]


genus_assignment <- cluster_mean_wide %>%
  mutate(
    Highest_mean =
      highest_mean_values,
    
    Second_highest_mean =
      second_highest_mean_values,
    
    Lowest_mean =
      lowest_mean_values,
    
    Assigned_cluster = factor(
      assigned_cluster_values,
      levels = cluster_order
    ),
    
    Leading_margin =
      Highest_mean -
      Second_highest_mean,
    
    Highest_to_second_ratio = if_else(
      Second_highest_mean > 0,
      Highest_mean /
        Second_highest_mean,
      NA_real_
    ),
    
    Overall_cluster_range =
      Highest_mean -
      Lowest_mean,
    
    Highest_mean_percent =
      Highest_mean * 100,
    
    Second_highest_mean_percent =
      Second_highest_mean * 100,
    
    Leading_margin_percentage_points =
      Leading_margin * 100,
    
    Pass_highest_mean_threshold =
      Highest_mean >=
      minimum_highest_mean,
    
    Pass_leading_margin_threshold =
      Leading_margin >=
      minimum_leading_margin,
    
    Qualified_for_display =
      Pass_highest_mean_threshold &
      Pass_leading_margin_threshold
  ) %>%
  arrange(
    Assigned_cluster,
    desc(
      Leading_margin
    ),
    desc(
      Highest_mean
    )
  )

write_tsv(
  genus_assignment,
  file.path(
    table_dir,
    "all_genera_assignment_and_qualification_audit_7KB.tsv"
  )
)


# ==============================================================================
# 12. Qualification audit by state
# ==============================================================================

assignment_summary <- genus_assignment %>%
  group_by(
    Assigned_cluster
  ) %>%
  summarise(
    All_assigned_genera_n =
      n(),
    
    Pass_highest_mean_n =
      sum(
        Pass_highest_mean_threshold,
        na.rm = TRUE
      ),
    
    Pass_leading_margin_n =
      sum(
        Pass_leading_margin_threshold,
        na.rm = TRUE
      ),
    
    Qualified_genera_n =
      sum(
        Qualified_for_display,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  arrange(
    Assigned_cluster
  )

write_tsv(
  assignment_summary,
  file.path(
    table_dir,
    "genera_qualification_counts_by_cluster_7KB.tsv"
  )
)


# ==============================================================================
# 13. Select maximum 13 qualified genera per state
# ==============================================================================

qualified_genera <- genus_assignment %>%
  filter(
    Qualified_for_display
  ) %>%
  arrange(
    Assigned_cluster,
    desc(
      Leading_margin
    ),
    desc(
      Highest_mean
    )
  )

selected_genera <- qualified_genera %>%
  group_by(
    Assigned_cluster
  ) %>%
  arrange(
    desc(
      Leading_margin
    ),
    desc(
      Highest_mean
    ),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = maximum_genera_per_cluster
  ) %>%
  ungroup() %>%
  arrange(
    Assigned_cluster,
    desc(
      Leading_margin
    ),
    desc(
      Highest_mean
    )
  ) %>%
  mutate(
    Rank_within_assigned_cluster =
      ave(
        seq_len(
          n()
        ),
        Assigned_cluster,
        FUN = seq_along
      ),
    
    Display_order =
      row_number()
  )

number_selected_genera <- nrow(
  selected_genera
)

if (
  number_selected_genera == 0
) {
  stop(
    paste0(
      "No genus passed the original display thresholds. ",
      "Do not alter thresholds before auditing results."
    ),
    call. = FALSE
  )
}

selected_counts <- selected_genera %>%
  count(
    Assigned_cluster,
    name = "Selected_genera_n"
  ) %>%
  complete(
    Assigned_cluster = factor(
      cluster_order,
      levels = cluster_order
    ),
    fill = list(
      Selected_genera_n = 0L
    )
  ) %>%
  arrange(
    Assigned_cluster
  )

write_tsv(
  selected_counts,
  file.path(
    table_dir,
    "selected_genera_counts_by_cluster_7KB.tsv"
  )
)


# ==============================================================================
# 14. Build genus-to-phylum mapping from 7KB taxonomy
# ==============================================================================

taxonomy_with_genus <- taxonomy_raw %>%
  transmute(
    ASV = as.character(
      ASV
    ),
    
    Taxonomy = as.character(
      Taxonomy
    ),
    
    Genus_strict =
      extract_genus_strict(
        Taxonomy
      )
  ) %>%
  filter(
    !is_bad_genus(
      Genus_strict
    )
  )

genus_phylum_candidates <- taxonomy_with_genus %>%
  transmute(
    Genus_key =
      canonicalize_genus(
        Genus_strict
      ),
    
    Phylum_raw =
      extract_phylum(
        Taxonomy
      )
  ) %>%
  filter(
    !is.na(
      Genus_key
    ),
    Genus_key != ""
  ) %>%
  count(
    Genus_key,
    Phylum_raw,
    name = "ASV_n"
  ) %>%
  group_by(
    Genus_key
  ) %>%
  arrange(
    desc(
      ASV_n
    ),
    .by_group = TRUE
  ) %>%
  slice_head(
    n = 1
  ) %>%
  ungroup() %>%
  mutate(
    Phylum_group =
      standardize_phylum_group(
        Phylum_raw
      )
  )

selected_genera <- selected_genera %>%
  mutate(
    Genus_key =
      canonicalize_genus(
        Genus
      )
  ) %>%
  left_join(
    genus_phylum_candidates %>%
      select(
        Genus_key,
        Phylum_raw,
        Phylum_group
      ),
    by = "Genus_key"
  ) %>%
  mutate(
    Phylum_raw = coalesce(
      Phylum_raw,
      "Unclassified"
    ),
    
    Phylum_group = coalesce(
      Phylum_group,
      "Other"
    )
  )

write_tsv(
  selected_genera,
  file.path(
    table_dir,
    "Figure4C_selected_qualified_genera_7KB.tsv"
  )
)

# Also save the old-style filename so that the final plotting script
# can be migrated with minimal structural changes.
write_tsv(
  selected_genera,
  file.path(
    table_dir,
    "Figure3B_selected_top13_per_cluster.tsv"
  )
)


# ==============================================================================
# 15. Row-normalize C1/C2/C3 mean abundance
# ==============================================================================

mean_matrix <- selected_genera %>%
  select(
    Mean_C1,
    Mean_C2,
    Mean_C3
  ) %>%
  as.matrix()

rownames(
  mean_matrix
) <- selected_genera$Genus

row_normalized_matrix <- t(
  apply(
    mean_matrix,
    1,
    row_normalize_three_values
  )
)

colnames(
  row_normalized_matrix
) <- cluster_order

rownames(
  row_normalized_matrix
) <- selected_genera$Genus


# ==============================================================================
# 16. Build long row-normalized heatmap source
# ==============================================================================

heatmap_long <- as.data.frame(
  row_normalized_matrix,
  check.names = FALSE
) %>%
  rownames_to_column(
    "Genus"
  ) %>%
  pivot_longer(
    cols = all_of(
      cluster_order
    ),
    names_to = "DMM_cluster",
    values_to = "Row_normalized_mean"
  ) %>%
  left_join(
    selected_genera %>%
      select(
        Genus,
        Assigned_cluster,
        Rank_within_assigned_cluster,
        Leading_margin,
        Leading_margin_percentage_points,
        Highest_mean,
        Highest_mean_percent,
        Highest_to_second_ratio,
        Display_order,
        Phylum_raw,
        Phylum_group
      ),
    by = "Genus"
  ) %>%
  mutate(
    DMM_cluster = factor(
      DMM_cluster,
      levels = cluster_order
    ),
    
    Row_position =
      number_selected_genera -
      Display_order +
      1L
  )

write_tsv(
  heatmap_long,
  file.path(
    table_dir,
    "Figure4C_row_normalized_cluster_mean_source_7KB.tsv"
  )
)

# Old-style filename for later final-script compatibility.
write_tsv(
  heatmap_long,
  file.path(
    table_dir,
    "Figure3B_row_normalized_cluster_mean_source.tsv"
  )
)


# ==============================================================================
# 17. Export unscaled selected cluster means
# ==============================================================================

selected_cluster_means <- selected_genera %>%
  select(
    Assigned_cluster,
    Rank_within_assigned_cluster,
    Genus,
    Mean_C1,
    Mean_C2,
    Mean_C3,
    Highest_mean,
    Highest_mean_percent,
    Leading_margin,
    Leading_margin_percentage_points,
    Highest_to_second_ratio,
    Phylum_raw,
    Phylum_group
  )

write_tsv(
  selected_cluster_means,
  file.path(
    table_dir,
    "Figure4C_selected_cluster_mean_abundance_7KB.tsv"
  )
)


# ==============================================================================
# 18. Audit heatmap only
#
# This is NOT the final manuscript Figure 4C.
# ==============================================================================

plot_data <- heatmap_long %>%
  mutate(
    Genus = factor(
      Genus,
      levels = rev(
        selected_genera$Genus
      )
    )
  )

audit_heatmap <- ggplot(
  plot_data,
  aes(
    x = DMM_cluster,
    y = Genus,
    fill = Row_normalized_mean
  )
) +
  geom_tile(
    linewidth = 0.3
  ) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(
      -row_z_limit,
      row_z_limit
    ),
    oob = scales::squish,
    name = "Row-scaled\nmean abundance"
  ) +
  labs(
    title = "7KB DMM state-defining genera: selection audit",
    subtitle = paste0(
      "Highest mean >= 0.5%; leading margin >= 0.5 percentage points; ",
      "maximum 13 genera/state"
    ),
    x = NULL,
    y = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    axis.text.x = element_text(
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 7,
      color = "black"
    ),
    legend.title = element_text(
      size = 8
    ),
    legend.text = element_text(
      size = 7
    )
  )

grDevices::pdf(
  file.path(
    figure_dir,
    "Figure4C_selection_audit_heatmap_7KB.pdf"
  ),
  width = 5,
  height = 8,
  useDingbats = FALSE
)

print(
  audit_heatmap
)

grDevices::dev.off()

ggsave(
  filename = file.path(
    figure_dir,
    "Figure4C_selection_audit_heatmap_7KB.png"
  ),
  plot = audit_heatmap,
  width = 5,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    figure_dir,
    "Figure4C_selection_audit_heatmap_7KB.tif"
  ),
  plot = audit_heatmap,
  width = 5,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)


# ==============================================================================
# 19. Console audit
# ==============================================================================

message("")
message("============================================================")
message("7KB qualified cluster-mean heatmap source completed.")
message("============================================================")

message("")
message(
  "Samples: ",
  nrow(
    sample_metadata
  )
)

message(
  "Genera in primary DMM prevalence10 matrix: ",
  ncol(
    relative_abundance_matrix
  )
)

message("")
message("Selection thresholds:")

message(
  "Highest cluster mean >= ",
  minimum_highest_mean * 100,
  "%"
)

message(
  "Leading margin >= ",
  minimum_leading_margin * 100,
  " percentage points"
)

message(
  "Maximum genera per state: ",
  maximum_genera_per_cluster
)

message("")
message("Qualification counts:")

print(
  tibble::as_tibble(
    assignment_summary
  )
)

message("")
message("Selected genera per cluster:")

print(
  tibble::as_tibble(
    selected_counts
  )
)

message("")
message(
  "Total selected genera: ",
  number_selected_genera
)

message("")
message("Selected genera:")

print(
  tibble::as_tibble(
    selected_genera %>%
      select(
        Assigned_cluster,
        Rank_within_assigned_cluster,
        Genus,
        Highest_mean_percent,
        Leading_margin_percentage_points,
        Highest_to_second_ratio,
        Phylum_group
      )
  ),
  n = Inf,
  width = Inf
)

message("")
message("Output directory:")
message(output_dir)