#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_02_DMM_characterize_K3_states.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Purpose:
## 1. Characterize C1/C2/C3 in fixed progression127
## 2. Summarize cluster composition and dominant genera
## 3. Generate row-scaled genus heatmap source
## 4. Compare genus-level ecological characteristics
## 5. Compare representative genera
## 6. Test association with progression stage
## 7. Audit technical variables
##
## IMPORTANT:
## - Uses the prevalence-10 genus RAW COUNT matrix used for DMM fitting.
## - Genus-level diversity here characterizes DMM states only.
## - It does NOT replace ASV-level alpha-diversity analyses.
## - C1/C2/C3 are retained exactly as native K=3 components because
##   the 7KB native states reproduce the original state identities.
############################################################


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "stringr",
  "tibble",
  "scales",
  "vegan",
  "pheatmap",
  "RColorBrewer"
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
      "Missing required package(s):\n",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(scales)
  library(vegan)
  library(pheatmap)
  library(RColorBrewer)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

dmm_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "DMM_progression127_primary_prevalence10_7KB"
)

input_table_dir <- file.path(
  dmm_dir,
  "tables"
)

genus_count_file <- file.path(
  input_table_dir,
  "genus_raw_count_progression127_prevalence10_frozen_7KB.tsv"
)

assignment_file <- file.path(
  input_table_dir,
  "DMM_sample_assignments_optimalK3_native_7KB.tsv"
)

posterior_file <- file.path(
  input_table_dir,
  "DMM_posterior_probabilities_optimalK3_7KB.tsv"
)

metadata_file <- file.path(
  input_table_dir,
  "metadata_progression127_frozen_7KB.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "DMM_progression127_K3_characterization_7KB"
)

table_dir <- file.path(
  output_dir,
  "tables"
)

figure_dir <- file.path(
  output_dir,
  "figures"
)

log_dir <- file.path(
  output_dir,
  "logs"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
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

dir.create(
  log_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Analysis parameters
#
# Retained from original 09 script.
# ==============================================================================

sample_id_col <- "SampleID"

group_col <- "Progression5"

cluster_col <- "DMM_cluster"

group_order <- c(
  "Polyp",
  "UC_remission",
  "UC_active",
  "Dysplasia",
  "CA"
)

cluster_order <- c(
  "C1",
  "C2",
  "C3"
)

n_heatmap_genera <- 30L

minimum_mean_relative_abundance <- 0.001

candidate_genera <- c(
  "UCG-005",
  "Streptococcus",
  "Desulfovibrio",
  "Mediterraneibacter",
  "Dorea",
  "Monoglobus",
  "[Eubacterium]_eligens_group",
  "Enterococcus",
  "Lactobacillus",
  "Ligilactobacillus",
  "Bacteroides"
)

set.seed(20260715)


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
        paste(
          missing_paths,
          collapse = "\n"
        ),
        sep = "\n"
      ),
      call. = FALSE
    )
  }
}


save_plot_all_formats <- function(
    plot_object,
    filename_base,
    width = 7,
    height = 5,
    dpi = 600
) {
  
  grDevices::pdf(
    paste0(
      filename_base,
      ".pdf"
    ),
    width = width,
    height = height,
    useDingbats = FALSE
  )
  
  print(plot_object)
  
  grDevices::dev.off()
  
  if (
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ) {
    
    ragg::agg_png(
      paste0(
        filename_base,
        ".png"
      ),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      background = "white"
    )
    
    print(plot_object)
    
    grDevices::dev.off()
    
    ragg::agg_tiff(
      paste0(
        filename_base,
        ".tif"
      ),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      background = "white"
    )
    
    print(plot_object)
    
    grDevices::dev.off()
    
  } else {
    
    if (Sys.info()[["sysname"]] != "Darwin") {
      stop(
        "Package 'ragg' is required for PNG/TIFF export on non-macOS systems.",
        call. = FALSE
      )
    }
    
    grDevices::png(
      paste0(
        filename_base,
        ".png"
      ),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      type = "quartz",
      bg = "white"
    )
    
    print(plot_object)
    
    grDevices::dev.off()
    
    grDevices::tiff(
      paste0(
        filename_base,
        ".tif"
      ),
      width = width,
      height = height,
      units = "in",
      res = dpi,
      type = "quartz",
      compression = "lzw",
      bg = "white"
    )
    
    print(plot_object)
    
    grDevices::dev.off()
  }
}


safe_kruskal <- function(
    data,
    outcome,
    group
) {
  
  formula_object <- reformulate(
    termlabels = group,
    response = outcome
  )
  
  fit <- kruskal.test(
    formula_object,
    data = data
  )
  
  tibble(
    Variable = outcome,
    Test = "Kruskal-Wallis",
    Statistic = unname(
      fit$statistic
    ),
    Degrees_of_freedom = unname(
      fit$parameter
    ),
    P_value = fit$p.value
  )
}


pairwise_wilcoxon_bh <- function(
    data,
    outcome,
    group
) {
  
  group_values <- unique(
    as.character(
      data[[group]]
    )
  )
  
  group_values <- group_values[
    !is.na(group_values)
  ]
  
  group_pairs <- combn(
    group_values,
    2,
    simplify = FALSE
  )
  
  results <- lapply(
    group_pairs,
    function(current_pair) {
      
      subset_data <- data %>%
        filter(
          .data[[group]] %in%
            current_pair
        )
      
      x <- subset_data %>%
        filter(
          .data[[group]] ==
            current_pair[1]
        ) %>%
        pull(
          .data[[outcome]]
        )
      
      y <- subset_data %>%
        filter(
          .data[[group]] ==
            current_pair[2]
        ) %>%
        pull(
          .data[[outcome]]
        )
      
      test_result <- wilcox.test(
        x,
        y,
        exact = FALSE
      )
      
      tibble(
        Variable = outcome,
        Group1 = current_pair[1],
        Group2 = current_pair[2],
        N1 = length(x),
        N2 = length(y),
        Median1 = median(
          x,
          na.rm = TRUE
        ),
        Median2 = median(
          y,
          na.rm = TRUE
        ),
        P_value = test_result$p.value
      )
    }
  )
  
  bind_rows(
    results
  ) %>%
    mutate(
      Q_value_BH = p.adjust(
        P_value,
        method = "BH"
      )
    )
}


scale_rows_safely <- function(
    matrix_object
) {
  
  row_means <- rowMeans(
    matrix_object,
    na.rm = TRUE
  )
  
  row_sds <- apply(
    matrix_object,
    1,
    sd,
    na.rm = TRUE
  )
  
  row_sds[
    is.na(row_sds) |
      row_sds == 0
  ] <- 1
  
  scaled_matrix <- sweep(
    matrix_object,
    1,
    row_means,
    "-"
  )
  
  scaled_matrix <- sweep(
    scaled_matrix,
    1,
    row_sds,
    "/"
  )
  
  scaled_matrix
}


# ==============================================================================
# 4. Validate and read input
# ==============================================================================

stop_if_missing(
  c(
    genus_count_file,
    assignment_file,
    posterior_file,
    metadata_file
  )
)

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

posterior_raw <- data.table::fread(
  posterior_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

metadata_raw <- data.table::fread(
  metadata_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)


# ==============================================================================
# 5. Standardize 7KB technical-variable names
#
# Original 09 expected:
# Reads_after_step2
# Retention_fraction
# High_contamination_risk
#
# 7KB frozen metadata contains the corresponding *_7KB_audit fields.
# We create aliases only; values are unchanged.
# ==============================================================================

if (
  !"Reads_after_step2" %in% names(metadata_raw) &&
  "Reads_after_step2_7KB_audit" %in% names(metadata_raw)
) {
  metadata_raw$Reads_after_step2 <-
    metadata_raw$Reads_after_step2_7KB_audit
}

if (
  !"Retention_fraction" %in% names(metadata_raw) &&
  "Retention_fraction_7KB_audit" %in% names(metadata_raw)
) {
  metadata_raw$Retention_fraction <-
    metadata_raw$Retention_fraction_7KB_audit
}

if (
  !"High_contamination_risk" %in% names(metadata_raw) &&
  "High_contamination_risk_7KB_audit" %in% names(metadata_raw)
) {
  metadata_raw$High_contamination_risk <-
    metadata_raw$High_contamination_risk_7KB_audit
}


# ==============================================================================
# 6. Validate required columns
# ==============================================================================

required_assignment_cols <- c(
  sample_id_col,
  cluster_col,
  "Maximum_posterior_probability"
)

missing_assignment_cols <- setdiff(
  required_assignment_cols,
  colnames(assignment_raw)
)

if (length(missing_assignment_cols) > 0) {
  stop(
    "Missing columns in assignment file: ",
    paste(
      missing_assignment_cols,
      collapse = ", "
    ),
    call. = FALSE
  )
}

required_metadata_cols <- c(
  sample_id_col,
  group_col,
  "Reads_after_step2",
  "Retention_fraction",
  "High_contamination_risk"
)

missing_metadata_cols <- setdiff(
  required_metadata_cols,
  colnames(metadata_raw)
)

if (length(missing_metadata_cols) > 0) {
  stop(
    "Missing columns in metadata file: ",
    paste(
      missing_metadata_cols,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  !sample_id_col %in%
  colnames(genus_count_raw)
) {
  stop(
    "The frozen genus count matrix does not contain column: ",
    sample_id_col,
    call. = FALSE
  )
}


# ==============================================================================
# 7. Build sample-by-genus count matrix
# ==============================================================================

sample_ids <- as.character(
  genus_count_raw[[sample_id_col]]
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
  select(
    -all_of(
      sample_id_col
    )
  ) %>%
  as.matrix()

storage.mode(
  genus_count_matrix
) <- "numeric"

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
    "Negative values detected in genus count matrix.",
    call. = FALSE
  )
}

rownames(
  genus_count_matrix
) <- sample_ids

if (
  any(
    rowSums(
      genus_count_matrix
    ) == 0
  )
) {
  stop(
    "At least one sample has zero genus-level counts.",
    call. = FALSE
  )
}


# ==============================================================================
# 8. Merge assignment and metadata
# ==============================================================================

metadata_selected <- metadata_raw %>%
  mutate(
    SampleID = as.character(
      .data[[sample_id_col]]
    )
  )

assignment_selected <- assignment_raw %>%
  mutate(
    SampleID = as.character(
      .data[[sample_id_col]]
    ),
    DMM_cluster = factor(
      .data[[cluster_col]],
      levels = cluster_order
    )
  ) %>%
  select(
    SampleID,
    DMM_cluster,
    DMM_cluster_number,
    Maximum_posterior_probability,
    Posterior_entropy,
    High_confidence_assignment
  )

sample_data <- metadata_selected %>%
  inner_join(
    assignment_selected,
    by = "SampleID"
  )

if (
  nrow(
    sample_data
  ) !=
  nrow(
    genus_count_matrix
  )
) {
  stop(
    paste0(
      "Merged metadata/assignment sample count does not match ",
      "genus count matrix sample count."
    ),
    call. = FALSE
  )
}

missing_from_counts <- setdiff(
  sample_data$SampleID,
  rownames(
    genus_count_matrix
  )
)

missing_from_metadata <- setdiff(
  rownames(
    genus_count_matrix
  ),
  sample_data$SampleID
)

if (
  length(
    missing_from_counts
  ) > 0 ||
  length(
    missing_from_metadata
  ) > 0
) {
  stop(
    "Sample IDs do not match perfectly across input files.",
    call. = FALSE
  )
}

sample_data <- sample_data %>%
  arrange(
    match(
      SampleID,
      rownames(
        genus_count_matrix
      )
    )
  )

genus_count_matrix <- genus_count_matrix[
  sample_data$SampleID,
  ,
  drop = FALSE
]


# ==============================================================================
# 9. Strict K=3 / progression127 audit
# ==============================================================================

cluster_size_audit <- sample_data %>%
  count(
    DMM_cluster,
    name = "N"
  )

expected_cluster_sizes <- c(
  C1 = 55,
  C2 = 46,
  C3 = 26
)

observed_cluster_sizes <- setNames(
  cluster_size_audit$N,
  as.character(
    cluster_size_audit$DMM_cluster
  )
)

if (
  !identical(
    as.integer(
      observed_cluster_sizes[
        names(
          expected_cluster_sizes
        )
      ]
    ),
    as.integer(
      expected_cluster_sizes
    )
  )
) {
  stop(
    paste0(
      "K=3 cluster-size audit failed. ",
      "Expected C1=55, C2=46, C3=26."
    ),
    call. = FALSE
  )
}

progression_audit <- sample_data %>%
  count(
    Progression5,
    name = "N"
  )

expected_progression_counts <- c(
  Polyp = 26,
  UC_remission = 36,
  UC_active = 25,
  Dysplasia = 17,
  CA = 23
)

observed_progression_counts <- setNames(
  progression_audit$N,
  progression_audit$Progression5
)

if (
  !identical(
    as.integer(
      observed_progression_counts[
        names(
          expected_progression_counts
        )
      ]
    ),
    as.integer(
      expected_progression_counts
    )
  )
) {
  stop(
    "progression127 group-count audit failed.",
    call. = FALSE
  )
}


# ==============================================================================
# 10. Set disease-group factor order
# ==============================================================================

observed_groups <- unique(
  as.character(
    sample_data[[group_col]]
  )
)

if (
  all(
    group_order %in%
    observed_groups
  )
) {
  
  sample_data[[group_col]] <- factor(
    sample_data[[group_col]],
    levels = group_order
  )
  
} else {
  
  message(
    "The predefined group labels did not exactly match metadata."
  )
  
  message(
    "Observed group labels: ",
    paste(
      sort(
        observed_groups
      ),
      collapse = ", "
    )
  )
  
  sample_data[[group_col]] <- factor(
    sample_data[[group_col]],
    levels = sort(
      observed_groups
    )
  )
}


# ==============================================================================
# 11. Calculate relative abundance
# ==============================================================================

relative_abundance_matrix <- sweep(
  genus_count_matrix,
  1,
  rowSums(
    genus_count_matrix
  ),
  "/"
)

if (
  anyNA(
    relative_abundance_matrix
  )
) {
  stop(
    "NA values generated during relative-abundance conversion.",
    call. = FALSE
  )
}

write_tsv(
  as.data.frame(
    relative_abundance_matrix,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "SampleID"
    ),
  file.path(
    table_dir,
    "DMM_K3_genus_relative_abundance_7KB.tsv"
  )
)


# ==============================================================================
# 12. Calculate genus-level diversity
#
# Exact original 09 metrics.
# ==============================================================================

genus_shannon <- vegan::diversity(
  genus_count_matrix,
  index = "shannon"
)

genus_simpson <- vegan::diversity(
  genus_count_matrix,
  index = "simpson"
)

observed_genera <- rowSums(
  genus_count_matrix > 0
)

dominant_genus_relative_abundance <- apply(
  relative_abundance_matrix,
  1,
  max
)

dominant_genus_name <- apply(
  relative_abundance_matrix,
  1,
  function(x) {
    names(x)[
      which.max(
        x
      )
    ]
  }
)

diversity_table <- sample_data %>%
  transmute(
    SampleID,
    DMM_cluster,
    Progression_group = .data[[group_col]],
    Genus_Shannon = genus_shannon[SampleID],
    Genus_Simpson = genus_simpson[SampleID],
    Observed_genera = observed_genera[SampleID],
    Dominant_genus = dominant_genus_name[SampleID],
    Dominant_genus_relative_abundance =
      dominant_genus_relative_abundance[SampleID],
    Reads_after_step2,
    Retention_fraction,
    High_contamination_risk
  )

write_tsv(
  diversity_table,
  file.path(
    table_dir,
    "DMM_K3_genus_diversity_and_dominance_7KB.tsv"
  )
)


# ==============================================================================
# 13. Diversity summary by DMM state
# ==============================================================================

diversity_summary <- diversity_table %>%
  group_by(
    DMM_cluster
  ) %>%
  summarise(
    N = n(),
    
    Genus_Shannon_median = median(
      Genus_Shannon
    ),
    
    Genus_Shannon_IQR = IQR(
      Genus_Shannon
    ),
    
    Genus_Simpson_median = median(
      Genus_Simpson
    ),
    
    Genus_Simpson_IQR = IQR(
      Genus_Simpson
    ),
    
    Observed_genera_median = median(
      Observed_genera
    ),
    
    Observed_genera_IQR = IQR(
      Observed_genera
    ),
    
    Dominant_abundance_median = median(
      Dominant_genus_relative_abundance
    ),
    
    Dominant_abundance_IQR = IQR(
      Dominant_genus_relative_abundance
    ),
    
    .groups = "drop"
  )

write_tsv(
  diversity_summary,
  file.path(
    table_dir,
    "DMM_K3_diversity_summary_by_cluster_7KB.tsv"
  )
)


# ==============================================================================
# 14. Overall and pairwise diversity tests
# ==============================================================================

diversity_variables <- c(
  "Genus_Shannon",
  "Genus_Simpson",
  "Observed_genera",
  "Dominant_genus_relative_abundance"
)

diversity_overall_tests <- bind_rows(
  lapply(
    diversity_variables,
    function(current_variable) {
      safe_kruskal(
        data = diversity_table,
        outcome = current_variable,
        group = "DMM_cluster"
      )
    }
  )
)

write_tsv(
  diversity_overall_tests,
  file.path(
    table_dir,
    "DMM_K3_diversity_Kruskal_Wallis_tests_7KB.tsv"
  )
)

diversity_pairwise_tests <- bind_rows(
  lapply(
    diversity_variables,
    function(current_variable) {
      pairwise_wilcoxon_bh(
        data = diversity_table,
        outcome = current_variable,
        group = "DMM_cluster"
      )
    }
  )
)

write_tsv(
  diversity_pairwise_tests,
  file.path(
    table_dir,
    "DMM_K3_diversity_pairwise_Wilcoxon_BH_7KB.tsv"
  )
)


# ==============================================================================
# 15. Audit diversity plot
#
# This is NOT yet the final Figure 4D.
# ==============================================================================

diversity_long <- diversity_table %>%
  select(
    SampleID,
    DMM_cluster,
    Genus_Shannon,
    Genus_Simpson,
    Observed_genera,
    Dominant_genus_relative_abundance
  ) %>%
  pivot_longer(
    cols = c(
      Genus_Shannon,
      Genus_Simpson,
      Observed_genera,
      Dominant_genus_relative_abundance
    ),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c(
        "Genus_Shannon",
        "Genus_Simpson",
        "Observed_genera",
        "Dominant_genus_relative_abundance"
      ),
      labels = c(
        "Genus-level Shannon",
        "Genus-level Simpson",
        "Observed genera",
        "Dominant-genus proportion"
      )
    )
  )

diversity_plot <- ggplot(
  diversity_long,
  aes(
    x = DMM_cluster,
    y = Value
  )
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.4,
    alpha = 0.65
  ) +
  facet_wrap(
    ~ Metric,
    scales = "free_y",
    ncol = 2
  ) +
  labs(
    title = "Ecological characteristics of DMM community states",
    subtitle = "Genus-level metrics calculated from the frozen DMM input matrix",
    x = "DMM community state",
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

save_plot_all_formats(
  plot_object = diversity_plot,
  filename_base = file.path(
    figure_dir,
    "DMM_K3_genus_diversity_and_dominance_7KB"
  ),
  width = 8,
  height = 7
)


# ==============================================================================
# 16. Calculate cluster mean relative abundance
# ==============================================================================

relative_abundance_df <- as.data.frame(
  relative_abundance_matrix,
  check.names = FALSE
) %>%
  rownames_to_column(
    "SampleID"
  ) %>%
  left_join(
    sample_data %>%
      select(
        SampleID,
        DMM_cluster
      ),
    by = "SampleID"
  )

cluster_mean_relative_abundance <- relative_abundance_df %>%
  group_by(
    DMM_cluster
  ) %>%
  summarise(
    across(
      -SampleID,
      \(x) mean(
        x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

write_tsv(
  cluster_mean_relative_abundance,
  file.path(
    table_dir,
    "DMM_K3_cluster_mean_genus_relative_abundance_7KB.tsv"
  )
)


# ==============================================================================
# 17. Convert cluster means to genus-by-cluster matrix
# ==============================================================================

cluster_mean_matrix <- cluster_mean_relative_abundance %>%
  column_to_rownames(
    "DMM_cluster"
  ) %>%
  as.matrix()

cluster_mean_matrix <- t(
  cluster_mean_matrix
)

cluster_mean_matrix <- cluster_mean_matrix[
  ,
  cluster_order[
    cluster_order %in%
      colnames(
        cluster_mean_matrix
      )
  ],
  drop = FALSE
]


# ==============================================================================
# 18. Rank genera by between-cluster contrast
# ==============================================================================

global_mean_abundance <- colMeans(
  relative_abundance_matrix
)

between_cluster_range <- apply(
  cluster_mean_matrix,
  1,
  function(x) {
    max(x) - min(x)
  }
)

between_cluster_sd <- apply(
  cluster_mean_matrix,
  1,
  sd
)

dominant_cluster <- apply(
  cluster_mean_matrix,
  1,
  function(x) {
    colnames(
      cluster_mean_matrix
    )[
      which.max(
        x
      )
    ]
  }
)

genus_cluster_ranking <- tibble(
  Genus = rownames(
    cluster_mean_matrix
  ),
  
  Global_mean_relative_abundance =
    global_mean_abundance[
      rownames(
        cluster_mean_matrix
      )
    ],
  
  Between_cluster_range =
    between_cluster_range,
  
  Between_cluster_SD =
    between_cluster_sd,
  
  Highest_mean_cluster =
    dominant_cluster
) %>%
  arrange(
    desc(
      Between_cluster_range
    ),
    desc(
      Global_mean_relative_abundance
    )
  )

write_tsv(
  genus_cluster_ranking,
  file.path(
    table_dir,
    "DMM_K3_genus_ranking_by_between_cluster_difference_7KB.tsv"
  )
)


# ==============================================================================
# 19. Export Top20 genera for each cluster
# ==============================================================================

cluster_mean_long <- as.data.frame(
  cluster_mean_matrix,
  check.names = FALSE
) %>%
  rownames_to_column(
    "Genus"
  ) %>%
  pivot_longer(
    cols = -Genus,
    names_to = "DMM_cluster",
    values_to = "Mean_relative_abundance"
  )

top_genera_by_cluster <- cluster_mean_long %>%
  group_by(
    DMM_cluster
  ) %>%
  arrange(
    desc(
      Mean_relative_abundance
    ),
    .by_group = TRUE
  ) %>%
  mutate(
    Rank_within_cluster =
      row_number()
  ) %>%
  filter(
    Rank_within_cluster <= 20
  ) %>%
  ungroup()

write_tsv(
  top_genera_by_cluster,
  file.path(
    table_dir,
    "DMM_K3_top20_genera_by_cluster_7KB.tsv"
  )
)


# ==============================================================================
# 20. Original 30-genus row-scaled heatmap source
#
# This reproduces the original 09 characterization step.
# It is NOT yet the curated final Figure 4C heatmap.
# ==============================================================================

heatmap_genera <- genus_cluster_ranking %>%
  slice_head(
    n = min(
      n_heatmap_genera,
      nrow(
        genus_cluster_ranking
      )
    )
  ) %>%
  pull(
    Genus
  )

heatmap_matrix_unscaled <- cluster_mean_matrix[
  heatmap_genera,
  ,
  drop = FALSE
]

heatmap_matrix_scaled <- scale_rows_safely(
  heatmap_matrix_unscaled
)

write_tsv(
  as.data.frame(
    heatmap_matrix_unscaled,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "Genus"
    ),
  file.path(
    table_dir,
    "DMM_K3_heatmap_top30_unscaled_cluster_means_7KB.tsv"
  )
)

write_tsv(
  as.data.frame(
    heatmap_matrix_scaled,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "Genus"
    ),
  file.path(
    table_dir,
    "DMM_K3_heatmap_top30_row_scaled_7KB.tsv"
  )
)


# ==============================================================================
# 21. Audit heatmap
# ==============================================================================

heatmap_annotation <- data.frame(
  Community_state = factor(
    colnames(
      heatmap_matrix_scaled
    ),
    levels = cluster_order
  )
)

rownames(
  heatmap_annotation
) <- colnames(
  heatmap_matrix_scaled
)

heatmap_breaks <- seq(
  -2,
  2,
  length.out = 101
)

heatmap_palette <- colorRampPalette(
  rev(
    RColorBrewer::brewer.pal(
      11,
      "RdBu"
    )
  )
)(100)

grDevices::pdf(
  file.path(
    figure_dir,
    "DMM_K3_top30_row_scaled_heatmap_7KB.pdf"
  ),
  width = 6.5,
  height = 10,
  useDingbats = FALSE
)

pheatmap(
  heatmap_matrix_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  scale = "none",
  color = heatmap_palette,
  breaks = heatmap_breaks,
  border_color = NA,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 11,
  angle_col = 0,
  annotation_col = heatmap_annotation,
  main = "DMM community-state defining genera"
)

grDevices::dev.off()


# ==============================================================================
# 22. Candidate representative genera
# ==============================================================================

available_candidate_genera <- intersect(
  candidate_genera,
  colnames(
    relative_abundance_matrix
  )
)

missing_candidate_genera <- setdiff(
  candidate_genera,
  colnames(
    relative_abundance_matrix
  )
)

write_lines(
  missing_candidate_genera,
  file.path(
    table_dir,
    "DMM_K3_candidate_genera_not_found_7KB.txt"
  )
)

if (
  length(
    available_candidate_genera
  ) < 8
) {
  
  additional_genera <- genus_cluster_ranking %>%
    filter(
      Global_mean_relative_abundance >=
        minimum_mean_relative_abundance
    ) %>%
    pull(
      Genus
    )
  
  available_candidate_genera <- unique(
    c(
      available_candidate_genera,
      additional_genera
    )
  )
  
  available_candidate_genera <-
    available_candidate_genera[
      seq_len(
        min(
          12,
          length(
            available_candidate_genera
          )
        )
      )
    ]
}

representative_genus_long <- as.data.frame(
  relative_abundance_matrix[
    ,
    available_candidate_genera,
    drop = FALSE
  ],
  check.names = FALSE
) %>%
  rownames_to_column(
    "SampleID"
  ) %>%
  left_join(
    sample_data %>%
      select(
        SampleID,
        DMM_cluster,
        all_of(
          group_col
        )
      ),
    by = "SampleID"
  ) %>%
  pivot_longer(
    cols = all_of(
      available_candidate_genera
    ),
    names_to = "Genus",
    values_to = "Relative_abundance"
  )

write_tsv(
  representative_genus_long,
  file.path(
    table_dir,
    "DMM_K3_representative_genus_sample_values_7KB.tsv"
  )
)


# ==============================================================================
# 23. Representative-genus statistical tests
# ==============================================================================

representative_genus_overall_tests <- representative_genus_long %>%
  group_by(
    Genus
  ) %>%
  group_modify(
    ~ {
      
      test_result <- kruskal.test(
        Relative_abundance ~ DMM_cluster,
        data = .x
      )
      
      tibble(
        Statistic = unname(
          test_result$statistic
        ),
        
        Degrees_of_freedom = unname(
          test_result$parameter
        ),
        
        P_value =
          test_result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    Q_value_BH = p.adjust(
      P_value,
      method = "BH"
    )
  )

write_tsv(
  representative_genus_overall_tests,
  file.path(
    table_dir,
    "DMM_K3_representative_genus_Kruskal_Wallis_BH_7KB.tsv"
  )
)

representative_genus_pairwise_tests <- bind_rows(
  lapply(
    available_candidate_genera,
    function(current_genus) {
      
      current_data <- representative_genus_long %>%
        filter(
          Genus ==
            current_genus
        )
      
      pairwise_wilcoxon_bh(
        data = current_data,
        outcome = "Relative_abundance",
        group = "DMM_cluster"
      ) %>%
        mutate(
          Genus =
            current_genus,
          .before = 1
        )
    }
  )
)

write_tsv(
  representative_genus_pairwise_tests,
  file.path(
    table_dir,
    "DMM_K3_representative_genus_pairwise_Wilcoxon_BH_7KB.tsv"
  )
)


# ==============================================================================
# 24. Representative-genus audit plot
# ==============================================================================

representative_genus_plot <- ggplot(
  representative_genus_long,
  aes(
    x = DMM_cluster,
    y = Relative_abundance
  )
) +
  geom_boxplot(
    width = 0.62,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.2,
    alpha = 0.6
  ) +
  facet_wrap(
    ~ Genus,
    scales = "free_y",
    ncol = 3
  ) +
  scale_y_continuous(
    labels = percent_format(
      accuracy = 0.1
    )
  ) +
  labs(
    title = "Representative genera across DMM community states",
    x = "DMM community state",
    y = "Relative abundance"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "italic"
    )
  )

save_plot_all_formats(
  plot_object =
    representative_genus_plot,
  
  filename_base =
    file.path(
      figure_dir,
      "DMM_K3_representative_genera_by_cluster_7KB"
    ),
  
  width = 10,
  height = 9
)


# ==============================================================================
# 25. Cluster distribution across progression groups
# ==============================================================================

cluster_by_group_counts <- sample_data %>%
  count(
    .data[[group_col]],
    DMM_cluster,
    name = "N"
  ) %>%
  group_by(
    .data[[group_col]]
  ) %>%
  mutate(
    Total_group =
      sum(
        N
      ),
    
    Proportion_within_group =
      N /
      Total_group
  ) %>%
  ungroup()

write_tsv(
  cluster_by_group_counts,
  file.path(
    table_dir,
    "DMM_K3_cluster_distribution_within_disease_group_7KB.tsv"
  )
)

group_by_cluster_counts <- sample_data %>%
  count(
    DMM_cluster,
    .data[[group_col]],
    name = "N"
  ) %>%
  group_by(
    DMM_cluster
  ) %>%
  mutate(
    Total_cluster =
      sum(
        N
      ),
    
    Proportion_within_cluster =
      N /
      Total_cluster
  ) %>%
  ungroup()

write_tsv(
  group_by_cluster_counts,
  file.path(
    table_dir,
    "DMM_K3_disease_group_distribution_within_cluster_7KB.tsv"
  )
)


# ==============================================================================
# 26. Progression distribution audit plot
# ==============================================================================

cluster_group_plot <- ggplot(
  cluster_by_group_counts,
  aes(
    x = .data[[group_col]],
    y = Proportion_within_group,
    fill = DMM_cluster
  )
) +
  geom_col(
    width = 0.75
  ) +
  scale_y_continuous(
    labels = percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0,
        0.04
      )
    )
  ) +
  labs(
    title = "DMM community states across the progression axis",
    x = NULL,
    y = "Proportion within disease group",
    fill = "Community state"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

save_plot_all_formats(
  plot_object =
    cluster_group_plot,
  
  filename_base =
    file.path(
      figure_dir,
      "DMM_K3_cluster_distribution_across_disease_groups_7KB"
    ),
  
  width = 8,
  height = 5.5
)


# ==============================================================================
# 27. Overall Fisher exact test
# ==============================================================================

cluster_group_contingency <- table(
  sample_data[[group_col]],
  sample_data$DMM_cluster
)

set.seed(
  20260715
)

overall_fisher <- fisher.test(
  cluster_group_contingency,
  simulate.p.value = TRUE,
  B = 100000
)

overall_fisher_summary <- tibble(
  Test =
    "Fisher exact test with Monte Carlo simulation",
  
  P_value =
    overall_fisher$p.value,
  
  Simulation_replicates =
    100000
)

write_tsv(
  overall_fisher_summary,
  file.path(
    table_dir,
    "DMM_K3_cluster_disease_group_overall_Fisher_7KB.tsv"
  )
)

write.table(
  cluster_group_contingency,
  file = file.path(
    table_dir,
    "DMM_K3_cluster_disease_group_contingency_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ==============================================================================
# 28. Cluster-specific enrichment tests
# ==============================================================================

group_levels <- levels(
  sample_data[[group_col]]
)

cluster_enrichment_results <- list()

result_index <- 1L

for (
  current_cluster in
  cluster_order
) {
  
  for (
    current_group in
    group_levels
  ) {
    
    cluster_binary <-
      sample_data$DMM_cluster ==
      current_cluster
    
    group_binary <-
      sample_data[[group_col]] ==
      current_group
    
    current_table <- table(
      Cluster_membership =
        cluster_binary,
      Disease_group =
        group_binary
    )
    
    current_fisher <- fisher.test(
      current_table
    )
    
    cluster_enrichment_results[[result_index]] <- tibble(
      DMM_cluster =
        current_cluster,
      
      Disease_group =
        current_group,
      
      N_cluster_and_group =
        sum(
          cluster_binary &
            group_binary
        ),
      
      N_cluster_total =
        sum(
          cluster_binary
        ),
      
      N_group_total =
        sum(
          group_binary
        ),
      
      Odds_ratio =
        unname(
          current_fisher$estimate
        ),
      
      CI_lower =
        current_fisher$conf.int[1],
      
      CI_upper =
        current_fisher$conf.int[2],
      
      P_value =
        current_fisher$p.value
    )
    
    result_index <-
      result_index + 1L
  }
}

cluster_enrichment_table <- bind_rows(
  cluster_enrichment_results
) %>%
  mutate(
    Q_value_BH = p.adjust(
      P_value,
      method = "BH"
    )
  )

write_tsv(
  cluster_enrichment_table,
  file.path(
    table_dir,
    "DMM_K3_cluster_specific_disease_group_enrichment_7KB.tsv"
  )
)


# ==============================================================================
# 29. Technical-variable summaries
#
# Same original characterization analysis, using 7KB audit columns.
# ==============================================================================

technical_summary <- sample_data %>%
  group_by(
    DMM_cluster
  ) %>%
  summarise(
    N = n(),
    
    Reads_after_step2_median =
      median(
        Reads_after_step2,
        na.rm = TRUE
      ),
    
    Reads_after_step2_IQR =
      IQR(
        Reads_after_step2,
        na.rm = TRUE
      ),
    
    Retention_fraction_median =
      median(
        Retention_fraction,
        na.rm = TRUE
      ),
    
    Retention_fraction_IQR =
      IQR(
        Retention_fraction,
        na.rm = TRUE
      ),
    
    High_contamination_risk_n =
      sum(
        High_contamination_risk %in%
          TRUE,
        na.rm = TRUE
      ),
    
    High_contamination_risk_proportion =
      mean(
        High_contamination_risk %in%
          TRUE,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

write_tsv(
  technical_summary,
  file.path(
    table_dir,
    "DMM_K3_technical_variable_summary_by_cluster_7KB.tsv"
  )
)


# ==============================================================================
# 30. Technical-variable tests
# ==============================================================================

technical_continuous_tests <- bind_rows(
  
  safe_kruskal(
    data = sample_data,
    outcome = "Reads_after_step2",
    group = "DMM_cluster"
  ),
  
  safe_kruskal(
    data = sample_data,
    outcome = "Retention_fraction",
    group = "DMM_cluster"
  )
)

write_tsv(
  technical_continuous_tests,
  file.path(
    table_dir,
    "DMM_K3_technical_continuous_Kruskal_Wallis_7KB.tsv"
  )
)

contamination_table <- table(
  sample_data$DMM_cluster,
  sample_data$High_contamination_risk
)

set.seed(
  20260715
)

contamination_fisher <- fisher.test(
  contamination_table,
  simulate.p.value = TRUE,
  B = 100000
)

contamination_fisher_summary <- tibble(
  Test =
    "Fisher exact test with Monte Carlo simulation",
  
  P_value =
    contamination_fisher$p.value,
  
  Simulation_replicates =
    100000
)

write_tsv(
  contamination_fisher_summary,
  file.path(
    table_dir,
    "DMM_K3_high_contamination_risk_Fisher_7KB.tsv"
  )
)

write.table(
  contamination_table,
  file = file.path(
    table_dir,
    "DMM_K3_high_contamination_risk_contingency_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ==============================================================================
# 31. Dominant genus frequencies
# ==============================================================================

dominant_genus_frequency <- diversity_table %>%
  count(
    DMM_cluster,
    Dominant_genus,
    name = "N"
  ) %>%
  group_by(
    DMM_cluster
  ) %>%
  mutate(
    Cluster_total =
      sum(
        N
      ),
    
    Proportion_within_cluster =
      N /
      Cluster_total
  ) %>%
  arrange(
    DMM_cluster,
    desc(
      N
    )
  ) %>%
  ungroup()

write_tsv(
  dominant_genus_frequency,
  file.path(
    table_dir,
    "DMM_K3_dominant_genus_frequency_by_cluster_7KB.tsv"
  )
)


# ==============================================================================
# 32. UCG-005 dominance audit
# ==============================================================================

if (
  "UCG-005" %in%
  colnames(
    relative_abundance_matrix
  )
) {
  
  ucg005_audit <- sample_data %>%
    transmute(
      SampleID,
      DMM_cluster,
      Progression_group =
        .data[[group_col]],
      
      UCG_005_relative_abundance =
        relative_abundance_matrix[
          SampleID,
          "UCG-005"
        ],
      
      Genus_Shannon =
        genus_shannon[
          SampleID
        ],
      
      Genus_Simpson =
        genus_simpson[
          SampleID
        ],
      
      Dominant_genus =
        dominant_genus_name[
          SampleID
        ],
      
      UCG_005_is_dominant =
        dominant_genus_name[
          SampleID
        ] ==
        "UCG-005"
    )
  
  write_tsv(
    ucg005_audit,
    file.path(
      table_dir,
      "DMM_K3_UCG005_dominance_audit_7KB.tsv"
    )
  )
  
  ucg005_summary <- ucg005_audit %>%
    group_by(
      DMM_cluster
    ) %>%
    summarise(
      N = n(),
      
      Median_UCG005 =
        median(
          UCG_005_relative_abundance
        ),
      
      IQR_UCG005 =
        IQR(
          UCG_005_relative_abundance
        ),
      
      UCG005_dominant_n =
        sum(
          UCG_005_is_dominant
        ),
      
      UCG005_dominant_proportion =
        mean(
          UCG_005_is_dominant
        ),
      
      .groups = "drop"
    )
  
  write_tsv(
    ucg005_summary,
    file.path(
      table_dir,
      "DMM_K3_UCG005_summary_by_cluster_7KB.tsv"
    )
  )
  
  ucg005_correlations <- tibble(
    Comparison = c(
      "UCG-005 vs genus-level Shannon",
      "UCG-005 vs genus-level Simpson",
      "UCG-005 vs dominant-genus proportion"
    ),
    
    Spearman_rho = c(
      cor(
        ucg005_audit$UCG_005_relative_abundance,
        ucg005_audit$Genus_Shannon,
        method = "spearman"
      ),
      
      cor(
        ucg005_audit$UCG_005_relative_abundance,
        ucg005_audit$Genus_Simpson,
        method = "spearman"
      ),
      
      cor(
        ucg005_audit$UCG_005_relative_abundance,
        diversity_table$Dominant_genus_relative_abundance,
        method = "spearman"
      )
    ),
    
    P_value = c(
      cor.test(
        ucg005_audit$UCG_005_relative_abundance,
        ucg005_audit$Genus_Shannon,
        method = "spearman",
        exact = FALSE
      )$p.value,
      
      cor.test(
        ucg005_audit$UCG_005_relative_abundance,
        ucg005_audit$Genus_Simpson,
        method = "spearman",
        exact = FALSE
      )$p.value,
      
      cor.test(
        ucg005_audit$UCG_005_relative_abundance,
        diversity_table$Dominant_genus_relative_abundance,
        method = "spearman",
        exact = FALSE
      )$p.value
    )
  ) %>%
    mutate(
      Q_value_BH = p.adjust(
        P_value,
        method = "BH"
      )
    )
  
  write_tsv(
    ucg005_correlations,
    file.path(
      table_dir,
      "DMM_K3_UCG005_ecological_correlations_7KB.tsv"
    )
  )
}


# ==============================================================================
# 33. Cluster characterization summary
# ==============================================================================

cluster_characterization_summary <- cluster_mean_long %>%
  group_by(
    DMM_cluster
  ) %>%
  slice_max(
    order_by =
      Mean_relative_abundance,
    n = 5,
    with_ties = FALSE
  ) %>%
  summarise(
    Top_five_genera = paste(
      Genus,
      collapse = "; "
    ),
    
    Top_five_mean_abundances = paste(
      scales::percent(
        Mean_relative_abundance,
        accuracy = 0.1
      ),
      collapse = "; "
    ),
    
    .groups = "drop"
  ) %>%
  left_join(
    diversity_summary,
    by = "DMM_cluster"
  ) %>%
  left_join(
    technical_summary,
    by = c(
      "DMM_cluster",
      "N"
    )
  )

write_tsv(
  cluster_characterization_summary,
  file.path(
    table_dir,
    "DMM_K3_cluster_characterization_summary_7KB.tsv"
  )
)


# ==============================================================================
# 34. Run summary
# ==============================================================================

run_summary <- tibble(
  Item = c(
    "Input frozen genus count matrix",
    "Input DMM assignment file",
    "Input posterior file",
    "Input metadata file",
    "Number of samples",
    "Number of genera",
    "Number of DMM states",
    "C1 samples",
    "C2 samples",
    "C3 samples",
    "Initial heatmap genera",
    "Overall cluster-disease Fisher P value",
    "Output directory"
  ),
  
  Value = c(
    genus_count_file,
    assignment_file,
    posterior_file,
    metadata_file,
    
    nrow(
      genus_count_matrix
    ),
    
    ncol(
      genus_count_matrix
    ),
    
    length(
      unique(
        sample_data$DMM_cluster
      )
    ),
    
    observed_cluster_sizes["C1"],
    observed_cluster_sizes["C2"],
    observed_cluster_sizes["C3"],
    
    length(
      heatmap_genera
    ),
    
    overall_fisher$p.value,
    
    output_dir
  )
)

write_tsv(
  run_summary,
  file.path(
    output_dir,
    "DMM_K3_characterization_run_summary_7KB.tsv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    log_dir,
    "sessionInfo_7KB.txt"
  )
)


# ==============================================================================
# 35. Console output
# ==============================================================================

message("")
message("============================================================")
message("7KB DMM K=3 characterization completed successfully.")
message("============================================================")

message("")
message("Cluster sizes:")
print(
  tibble::as_tibble(
    cluster_size_audit
  )
)

message("")
message("Genus-level diversity summary:")
print(
  tibble::as_tibble(
    diversity_summary
  )
)

message("")
message("Overall diversity / dominance tests:")
print(
  tibble::as_tibble(
    diversity_overall_tests
  )
)

message("")
message("Pairwise diversity / dominance tests:")
print(
  tibble::as_tibble(
    diversity_pairwise_tests
  )
)

message("")
message("Technical-variable tests:")
print(
  tibble::as_tibble(
    technical_continuous_tests
  )
)

message("")
message(
  "Overall cluster-disease Fisher P value: ",
  signif(
    overall_fisher$p.value,
    4
  )
)

message("")
message("Top five genera by cluster:")
print(
  tibble::as_tibble(
    top_genera_by_cluster %>%
      filter(
        Rank_within_cluster <= 5
      )
  )
)

message("")
message("Top 30 genera by between-cluster contrast:")
print(
  tibble::as_tibble(
    genus_cluster_ranking %>%
      slice_head(
        n = 30
      )
  )
)

message("")
message("UCG-005 summary by cluster:")

if (
  exists(
    "ucg005_summary"
  )
) {
  
  print(
    tibble::as_tibble(
      ucg005_summary
    )
  )
  
} else {
  
  message(
    "UCG-005 was not present in the prevalence-filtered matrix."
  )
}

message("")
message("Results saved to:")
message(output_dir)