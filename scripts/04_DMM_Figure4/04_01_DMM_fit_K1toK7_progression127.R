#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_01_DMM_fit_K1toK7_progression127.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Script:
##
## Purpose:
## 1. Read fixed 7KB progression127 ASV raw-count table
## 2. Aggregate ASVs to strict genus level
## 3. Freeze genus-level raw-count matrix
## 4. Apply prevalence filtering: >=10% samples
## 5. Fit DMM models for K = 1:7
## 6. Repeat each K 20 times
## 7. Retain best converged model by minimum Laplace
## 8. Export optimal-K assignments/posteriors/model objects
##
## IMPORTANT:
## - DMM input = RAW INTEGER COUNTS
## - No relative abundance is used for model fitting
## - No manual C1/C2/C3 biological relabeling is performed here
############################################################


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "stringr",
  "tibble",
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
      "Missing required package(s):\n",
      paste(
        missing_packages,
        collapse = ", "
      )
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
  library(DirichletMultinomial)
})


# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_root <- PROJECT_ROOT

input_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "00_clean_data",
  "progression127"
)

count_file <- file.path(
  input_dir,
  "asv_count_7KB_progression127.tsv"
)

taxonomy_file <- file.path(
  input_dir,
  "taxonomy_7KB_progression127.tsv"
)

metadata_file <- file.path(
  input_dir,
  "metadata_7KB_progression127.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "DMM_progression127_primary_prevalence10_7KB"
)

model_dir <- file.path(
  output_dir,
  "models"
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
  model_dir,
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


# ------------------------------------------------------------
# 2. Primary DMM settings
#
# EXACT original settings
# ------------------------------------------------------------

prevalence_threshold <- 0.10

minimum_count_for_presence <- 1L

k_values <- 1:7

n_repeats <- 20L

master_seed <- 20260715L

sample_id_col <- "SampleID"

group_col <- "Progression5"


# ------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------

stop_if_missing <- function(paths) {
  
  missing_paths <- paths[
    !file.exists(paths)
  ]
  
  if (length(missing_paths) > 0) {
    
    stop(
      paste(
        "The following required input files do not exist:",
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


extract_genus_strict <- function(x) {
  
  x <- as.character(x)
  
  x <- stringr::str_trim(
    x
  )
  
  genus <- rep(
    NA_character_,
    length(x)
  )
  
  has_g <- stringr::str_detect(
    x,
    "g__"
  )
  
  if (
    any(
      has_g,
      na.rm = TRUE
    )
  ) {
    
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
      g_part == "" |
        is.na(g_part)
    ] <- NA_character_
    
    genus[has_g] <- g_part
  }
  
  no_g_idx <- which(
    !has_g |
      is.na(has_g)
  )
  
  if (
    length(no_g_idx) > 0
  ) {
    
    candidate <- x[
      no_g_idx
    ]
    
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
    
    genus[
      no_g_idx
    ] <- candidate
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
    stringr::str_detect(
      x0,
      "uncultured"
    ) |
    stringr::str_detect(
      x0,
      "unclassified"
    ) |
    stringr::str_detect(
      x0,
      "unknown"
    ) |
    stringr::str_detect(
      x0,
      "metagenome"
    ) |
    stringr::str_detect(
      x0,
      "norank"
    ) |
    stringr::str_detect(
      x0,
      "no_rank"
    ) |
    stringr::str_detect(
      x0,
      "ambiguous"
    ) |
    stringr::str_detect(
      x0,
      "^[dkpcofs]__"
    ) |
    stringr::str_detect(
      x0,
      "k__|p__|c__|o__|f__|s__"
    )
}


extract_model_metrics <- function(
    fit_object
) {
  
  tibble(
    Laplace = as.numeric(
      laplace(
        fit_object
      )
    ),
    
    AIC = as.numeric(
      AIC(
        fit_object
      )
    ),
    
    BIC = as.numeric(
      BIC(
        fit_object
      )
    )
  )
}


fit_single_dmm <- function(
    count_matrix,
    k,
    seed_value
) {
  
  set.seed(
    seed_value
  )
  
  fit <- tryCatch(
    {
      dmn(
        count = count_matrix,
        k = k,
        verbose = FALSE
      )
    },
    
    error = function(e) {
      
      message(
        "DMM failed for K = ",
        k,
        ", seed = ",
        seed_value,
        ": ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  fit
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
  
  print(
    plot_object
  )
  
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
    
    print(
      plot_object
    )
    
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
    
    print(
      plot_object
    )
    
    grDevices::dev.off()
    
  } else {
    
    if (
      Sys.info()[["sysname"]] != "Darwin"
    ) {
      stop(
        "Package 'ragg' is required for PNG/TIFF export on non-macOS systems."
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
    
    print(
      plot_object
    )
    
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
      compression = "lzw",
      type = "quartz",
      bg = "white"
    )
    
    print(
      plot_object
    )
    
    grDevices::dev.off()
  }
}


# ------------------------------------------------------------
# 4. Validate input files
# ------------------------------------------------------------

stop_if_missing(
  c(
    count_file,
    taxonomy_file,
    metadata_file
  )
)

message("Input directory:")
message(input_dir)

message("")
message("Output directory:")
message(output_dir)


# ------------------------------------------------------------
# 5. Read data
# ------------------------------------------------------------

asv_count_raw <- data.table::fread(
  count_file,
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

metadata_raw <- data.table::fread(
  metadata_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)


# ------------------------------------------------------------
# 6. Validate required columns
# ------------------------------------------------------------

if (
  !sample_id_col %in%
  colnames(metadata_raw)
) {
  
  stop(
    "Metadata does not contain SampleID.",
    call. = FALSE
  )
}

if (
  !group_col %in%
  colnames(metadata_raw)
) {
  
  stop(
    "Metadata does not contain Progression5.",
    call. = FALSE
  )
}

if (
  !"ASV" %in%
  colnames(asv_count_raw)
) {
  
  stop(
    "ASV count table does not contain ASV.",
    call. = FALSE
  )
}

if (
  !"ASV" %in%
  colnames(taxonomy_raw)
) {
  
  stop(
    "Taxonomy table does not contain ASV.",
    call. = FALSE
  )
}

if (
  !"Taxonomy" %in%
  colnames(taxonomy_raw)
) {
  
  stop(
    "Taxonomy table does not contain Taxonomy.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 7. Strict progression127 audit
# ------------------------------------------------------------

metadata_raw <- metadata_raw %>%
  mutate(
    SampleID = as.character(
      SampleID
    ),
    
    Progression5 = as.character(
      Progression5
    )
  )

group_order <- c(
  "Polyp",
  "UC_remission",
  "UC_active",
  "Dysplasia",
  "CA"
)

metadata_raw <- metadata_raw %>%
  filter(
    Progression5 %in%
      group_order
  )

group_counts <- metadata_raw %>%
  count(
    Progression5,
    name = "N"
  )

expected_counts <- tibble(
  Progression5 = group_order,
  Expected_N = c(
    26,
    36,
    25,
    17,
    23
  )
)

group_count_audit <- expected_counts %>%
  left_join(
    group_counts,
    by = "Progression5"
  )

print(
  group_count_audit,
  n = Inf
)

if (
  nrow(metadata_raw) != 127 ||
  any(is.na(group_count_audit$N)) ||
  any(
    group_count_audit$N !=
    group_count_audit$Expected_N
  )
) {
  
  stop(
    "progression127 audit failed.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 8. Check sample matching
# ------------------------------------------------------------

count_sample_ids <- setdiff(
  colnames(asv_count_raw),
  "ASV"
)

metadata_sample_ids <- metadata_raw$SampleID

if (
  anyDuplicated(
    count_sample_ids
  ) > 0
) {
  
  stop(
    "Duplicated sample IDs detected in count table.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    metadata_sample_ids
  ) > 0
) {
  
  stop(
    "Duplicated sample IDs detected in metadata.",
    call. = FALSE
  )
}

samples_only_in_count <- setdiff(
  count_sample_ids,
  metadata_sample_ids
)

samples_only_in_metadata <- setdiff(
  metadata_sample_ids,
  count_sample_ids
)

sample_match_audit <- tibble(
  Metric = c(
    "Count table samples",
    "Metadata samples",
    "Samples only in count table",
    "Samples only in metadata"
  ),
  
  Value = c(
    length(count_sample_ids),
    length(metadata_sample_ids),
    length(samples_only_in_count),
    length(samples_only_in_metadata)
  )
)

write_tsv(
  sample_match_audit,
  file.path(
    table_dir,
    "sample_match_audit_7KB.tsv"
  )
)

if (
  length(samples_only_in_count) > 0 ||
  length(samples_only_in_metadata) > 0
) {
  
  stop(
    "Sample IDs do not match perfectly between count table and metadata.",
    call. = FALSE
  )
}


# ------------------------------------------------------------
# 9. Validate ASV raw-count matrix
# ------------------------------------------------------------

asv_ids <- as.character(
  asv_count_raw$ASV
)

if (
  anyDuplicated(
    asv_ids
  ) > 0
) {
  
  stop(
    "Duplicated ASV IDs detected.",
    call. = FALSE
  )
}

count_matrix_numeric <- asv_count_raw %>%
  select(
    all_of(
      count_sample_ids
    )
  ) %>%
  as.matrix()

storage.mode(
  count_matrix_numeric
) <- "numeric"

if (
  anyNA(
    count_matrix_numeric
  )
) {
  
  stop(
    "NA values detected in ASV count matrix.",
    call. = FALSE
  )
}

if (
  any(
    count_matrix_numeric < 0
  )
) {
  
  stop(
    "Negative counts detected.",
    call. = FALSE
  )
}

non_integer_values <- abs(
  count_matrix_numeric -
    round(
      count_matrix_numeric
    )
) > 1e-8

if (
  any(
    non_integer_values
  )
) {
  
  stop(
    "Non-integer values detected. DMM requires raw counts.",
    call. = FALSE
  )
}

count_matrix_numeric <- round(
  count_matrix_numeric
)

storage.mode(
  count_matrix_numeric
) <- "integer"

rownames(
  count_matrix_numeric
) <- asv_ids


# ------------------------------------------------------------
# 10. Taxonomy -> strict genus
# ------------------------------------------------------------

taxonomy_clean <- taxonomy_raw %>%
  transmute(
    ASV = as.character(
      ASV
    ),
    
    Taxonomy = as.character(
      Taxonomy
    ),
    
    Genus = extract_genus_strict(
      Taxonomy
    )
  ) %>%
  mutate(
    Genus_valid =
      !is_bad_genus(
        Genus
      )
  ) %>%
  filter(
    Genus_valid
  ) %>%
  distinct(
    ASV,
    .keep_all = TRUE
  )

taxonomy_coverage_summary <- tibble(
  Metric = c(
    "ASVs in count table",
    "ASVs with valid strict genus",
    "ASVs excluded from genus aggregation"
  ),
  
  Value = c(
    nrow(
      count_matrix_numeric
    ),
    
    sum(
      rownames(
        count_matrix_numeric
      ) %in%
        taxonomy_clean$ASV
    ),
    
    sum(
      !rownames(
        count_matrix_numeric
      ) %in%
        taxonomy_clean$ASV
    )
  )
)

write_tsv(
  taxonomy_coverage_summary,
  file.path(
    table_dir,
    "taxonomy_coverage_summary_7KB.tsv"
  )
)

message("")
message("Taxonomy coverage:")

print(
  taxonomy_coverage_summary,
  n = Inf
)


# ------------------------------------------------------------
# 11. Aggregate ASVs -> genus raw counts
# ------------------------------------------------------------

asv_count_long <- as.data.frame(
  count_matrix_numeric,
  check.names = FALSE
) %>%
  rownames_to_column(
    "ASV"
  ) %>%
  inner_join(
    taxonomy_clean %>%
      select(
        ASV,
        Genus
      ),
    by = "ASV"
  )

genus_count_feature_by_sample <- asv_count_long %>%
  select(
    Genus,
    all_of(
      count_sample_ids
    )
  ) %>%
  group_by(
    Genus
  ) %>%
  summarise(
    across(
      all_of(
        count_sample_ids
      ),
      ~ sum(
        .x,
        na.rm = TRUE
      )
    ),
    .groups = "drop"
  )

genus_ids <- genus_count_feature_by_sample$Genus

genus_count_matrix_feature_by_sample <-
  genus_count_feature_by_sample %>%
  select(
    all_of(
      count_sample_ids
    )
  ) %>%
  as.matrix()

storage.mode(
  genus_count_matrix_feature_by_sample
) <- "integer"

rownames(
  genus_count_matrix_feature_by_sample
) <- genus_ids


# ------------------------------------------------------------
# 12. Reorder samples to metadata order
# ------------------------------------------------------------

metadata_ordered <- metadata_raw %>%
  mutate(
    SampleID_internal =
      as.character(
        SampleID
      )
  )

sample_order <- metadata_ordered$SampleID_internal

genus_count_matrix_feature_by_sample <-
  genus_count_matrix_feature_by_sample[
    ,
    sample_order,
    drop = FALSE
  ]


# ------------------------------------------------------------
# 13. Remove zero-sum genera/samples
# ------------------------------------------------------------

genus_total_counts <- rowSums(
  genus_count_matrix_feature_by_sample
)

sample_total_counts <- colSums(
  genus_count_matrix_feature_by_sample
)

zero_sum_samples <- names(
  sample_total_counts[
    sample_total_counts == 0
  ]
)

if (
  length(
    zero_sum_samples
  ) > 0
) {
  
  stop(
    "Samples with zero genus-level counts detected: ",
    paste(
      zero_sum_samples,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  any(
    genus_total_counts == 0
  )
) {
  
  genus_count_matrix_feature_by_sample <-
    genus_count_matrix_feature_by_sample[
      genus_total_counts > 0,
      ,
      drop = FALSE
    ]
}


# ------------------------------------------------------------
# 14. Freeze unfiltered genus matrix
# ------------------------------------------------------------

genus_count_unfiltered_sample_by_genus <-
  t(
    genus_count_matrix_feature_by_sample
  )

storage.mode(
  genus_count_unfiltered_sample_by_genus
) <- "integer"

write_tsv(
  as.data.frame(
    genus_count_unfiltered_sample_by_genus,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "SampleID"
    ),
  file.path(
    table_dir,
    "genus_raw_count_progression127_unfiltered_7KB.tsv"
  )
)


# ------------------------------------------------------------
# 15. Primary prevalence filter: >=10%
# ------------------------------------------------------------

presence_matrix <-
  genus_count_unfiltered_sample_by_genus >=
  minimum_count_for_presence

genus_prevalence_n <- colSums(
  presence_matrix
)

genus_prevalence_fraction <-
  genus_prevalence_n /
  nrow(
    genus_count_unfiltered_sample_by_genus
  )

minimum_samples_required <- ceiling(
  prevalence_threshold *
    nrow(
      genus_count_unfiltered_sample_by_genus
    )
)

keep_genera <- names(
  genus_prevalence_n[
    genus_prevalence_n >=
      minimum_samples_required
  ]
)

if (
  length(
    keep_genera
  ) < 2
) {
  
  stop(
    "Fewer than 2 genera passed prevalence filtering.",
    call. = FALSE
  )
}

genus_count_filtered_sample_by_genus <-
  genus_count_unfiltered_sample_by_genus[
    ,
    keep_genera,
    drop = FALSE
  ]

storage.mode(
  genus_count_filtered_sample_by_genus
) <- "integer"


# ------------------------------------------------------------
# 16. Export filtering audit
# ------------------------------------------------------------

genus_filter_audit <- tibble(
  Genus = colnames(
    genus_count_unfiltered_sample_by_genus
  ),
  
  Prevalence_n = as.integer(
    genus_prevalence_n
  ),
  
  Prevalence_fraction = as.numeric(
    genus_prevalence_fraction
  ),
  
  Total_count = as.numeric(
    colSums(
      genus_count_unfiltered_sample_by_genus
    )
  ),
  
  Mean_count = as.numeric(
    colMeans(
      genus_count_unfiltered_sample_by_genus
    )
  ),
  
  Retained_primary_prevalence10 =
    Genus %in%
    keep_genera
) %>%
  arrange(
    desc(
      Prevalence_fraction
    ),
    desc(
      Total_count
    )
  )

write_tsv(
  genus_filter_audit,
  file.path(
    table_dir,
    "genus_prevalence_filter_audit_7KB.tsv"
  )
)

filter_summary <- tibble(
  Parameter = c(
    "Number of samples",
    "Genera before prevalence filtering",
    "Primary prevalence threshold",
    "Minimum count defining presence",
    "Minimum samples required",
    "Genera retained after filtering"
  ),
  
  Value = c(
    nrow(
      genus_count_unfiltered_sample_by_genus
    ),
    
    ncol(
      genus_count_unfiltered_sample_by_genus
    ),
    
    prevalence_threshold,
    
    minimum_count_for_presence,
    
    minimum_samples_required,
    
    ncol(
      genus_count_filtered_sample_by_genus
    )
  )
)

write_tsv(
  filter_summary,
  file.path(
    table_dir,
    "primary_filter_summary_7KB.tsv"
  )
)

message("")
message("Primary prevalence filter summary:")

print(
  filter_summary,
  n = Inf
)


# ------------------------------------------------------------
# 17. Freeze filtered genus matrix + metadata
# ------------------------------------------------------------

write_tsv(
  as.data.frame(
    genus_count_filtered_sample_by_genus,
    check.names = FALSE
  ) %>%
    rownames_to_column(
      "SampleID"
    ),
  file.path(
    table_dir,
    "genus_raw_count_progression127_prevalence10_frozen_7KB.tsv"
  )
)

write_tsv(
  metadata_ordered,
  file.path(
    table_dir,
    "metadata_progression127_frozen_7KB.tsv"
  )
)


# ------------------------------------------------------------
# 18. DMM fitting: K=1:7 x 20 repeats
# ------------------------------------------------------------

set.seed(
  master_seed
)

seed_grid <- tidyr::expand_grid(
  K = k_values,
  Repeat = seq_len(
    n_repeats
  )
) %>%
  mutate(
    Seed =
      master_seed +
      K * 10000L +
      Repeat
  )

all_fit_records <- vector(
  mode = "list",
  length = nrow(
    seed_grid
  )
)

all_fit_objects <- vector(
  mode = "list",
  length = nrow(
    seed_grid
  )
)

message("")
message("Starting DMM fitting...")

message(
  "Samples: ",
  nrow(
    genus_count_filtered_sample_by_genus
  )
)

message(
  "Genera: ",
  ncol(
    genus_count_filtered_sample_by_genus
  )
)

message(
  "K values: ",
  paste(
    k_values,
    collapse = ", "
  )
)

message(
  "Repeats per K: ",
  n_repeats
)

for (
  i in seq_len(
    nrow(
      seed_grid
    )
  )
) {
  
  current_k <- seed_grid$K[i]
  
  current_repeat <- seed_grid$Repeat[i]
  
  current_seed <- seed_grid$Seed[i]
  
  message(
    "Fitting K = ",
    current_k,
    ", repeat = ",
    current_repeat,
    "/",
    n_repeats
  )
  
  fit_object <- fit_single_dmm(
    count_matrix =
      genus_count_filtered_sample_by_genus,
    
    k = current_k,
    
    seed_value =
      current_seed
  )
  
  all_fit_objects[[i]] <-
    fit_object
  
  if (
    is.null(
      fit_object
    )
  ) {
    
    all_fit_records[[i]] <- tibble(
      K = current_k,
      Repeat = current_repeat,
      Seed = current_seed,
      Converged = FALSE,
      Laplace = NA_real_,
      AIC = NA_real_,
      BIC = NA_real_
    )
    
  } else {
    
    metrics <- extract_model_metrics(
      fit_object
    )
    
    all_fit_records[[i]] <- tibble(
      K = current_k,
      Repeat = current_repeat,
      Seed = current_seed,
      Converged = TRUE,
      Laplace = metrics$Laplace,
      AIC = metrics$AIC,
      BIC = metrics$BIC
    )
  }
}

fit_results <- bind_rows(
  all_fit_records
)

write_tsv(
  fit_results,
  file.path(
    table_dir,
    "DMM_all_repeats_K1_to_K7_7KB.tsv"
  )
)


# ------------------------------------------------------------
# 19. Select best repeat for each K
# ------------------------------------------------------------

best_repeat_by_k <- fit_results %>%
  filter(
    Converged,
    is.finite(
      Laplace
    )
  ) %>%
  group_by(
    K
  ) %>%
  slice_min(
    order_by = Laplace,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

if (
  !all(
    k_values %in%
    best_repeat_by_k$K
  )
) {
  
  missing_k <- setdiff(
    k_values,
    best_repeat_by_k$K
  )
  
  stop(
    "No converged DMM model available for K = ",
    paste(
      missing_k,
      collapse = ", "
    ),
    call. = FALSE
  )
}

write_tsv(
  best_repeat_by_k,
  file.path(
    table_dir,
    "DMM_best_repeat_per_K_7KB.tsv"
  )
)


# ------------------------------------------------------------
# 20. Retrieve/save best model for each K
# ------------------------------------------------------------

best_model_list <- list()

for (
  current_k in
  k_values
) {
  
  best_row <- best_repeat_by_k %>%
    filter(
      K ==
        current_k
    )
  
  object_index <- which(
    seed_grid$K ==
      current_k &
      seed_grid$Repeat ==
      best_row$Repeat
  )
  
  if (
    length(
      object_index
    ) != 1L
  ) {
    
    stop(
      "Could not uniquely identify best model object for K = ",
      current_k,
      call. = FALSE
    )
  }
  
  best_model <- all_fit_objects[[object_index]]
  
  best_model_list[[paste0("K", current_k)]] <- best_model
  
  saveRDS(
    best_model,
    file.path(
      model_dir,
      paste0(
        "DMM_best_model_K",
        current_k,
        "_7KB.rds"
      )
    )
  )
}

saveRDS(
  best_model_list,
  file.path(
    model_dir,
    "DMM_best_models_K1_to_K7_7KB.rds"
  )
)


# ------------------------------------------------------------
# 21. Select optimal K by minimum Laplace
# ------------------------------------------------------------

optimal_row <- best_repeat_by_k %>%
  slice_min(
    order_by = Laplace,
    n = 1,
    with_ties = FALSE
  )

optimal_k <- optimal_row$K

message("")
message(
  "Optimal K by minimum Laplace = ",
  optimal_k
)

write_tsv(
  optimal_row,
  file.path(
    table_dir,
    "DMM_optimal_K_by_minimum_Laplace_7KB.tsv"
  )
)

optimal_model <- best_model_list[[paste0("K", optimal_k)]]

saveRDS(
  optimal_model,
  file.path(
    model_dir,
    paste0(
      "DMM_primary_optimal_model_K",
      optimal_k,
      "_7KB.rds"
    )
  )
)


# ------------------------------------------------------------
# 22. Laplace curve
# ------------------------------------------------------------

laplace_plot_data <- best_repeat_by_k %>%
  arrange(
    K
  )

laplace_plot <- ggplot(
  laplace_plot_data,
  aes(
    x = K,
    y = Laplace
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_vline(
    xintercept =
      optimal_k,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_x_continuous(
    breaks =
      k_values
  ) +
  labs(
    title =
      "Dirichlet-multinomial mixture model selection",
    
    subtitle = paste0(
      "7KB progression127; prevalence >= ",
      prevalence_threshold * 100,
      "%; ",
      n_repeats,
      " repeated fits per K"
    ),
    
    x =
      "Number of DMM components (K)",
    
    y =
      "Minimum Laplace approximation"
  ) +
  theme_classic(
    base_size = 12
  )

save_plot_all_formats(
  plot_object =
    laplace_plot,
  
  filename_base =
    file.path(
      figure_dir,
      "DMM_Laplace_curve_K1_to_K7_7KB"
    ),
  
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 23. Repeat-level diagnostics
# ------------------------------------------------------------

repeat_plot <- ggplot(
  fit_results %>%
    filter(
      Converged,
      is.finite(
        Laplace
      )
    ),
  
  aes(
    x = factor(
      K
    ),
    y = Laplace
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.12,
    height = 0,
    size = 1.5,
    alpha = 0.65
  ) +
  labs(
    title =
      "Repeat-level DMM fitting diagnostics",
    
    subtitle = paste0(
      n_repeats,
      " repeated fits for each K"
    ),
    
    x =
      "Number of DMM components (K)",
    
    y =
      "Laplace approximation"
  ) +
  theme_classic(
    base_size = 12
  )

save_plot_all_formats(
  plot_object =
    repeat_plot,
  
  filename_base =
    file.path(
      figure_dir,
      "DMM_repeat_level_Laplace_diagnostics_7KB"
    ),
  
  width = 7,
  height = 5
)


# ------------------------------------------------------------
# 24. Posterior probabilities + native component assignments
# ------------------------------------------------------------

posterior_matrix <- mixture(
  optimal_model
)

if (
  nrow(
    posterior_matrix
  ) !=
  nrow(
    genus_count_filtered_sample_by_genus
  )
) {
  
  stop(
    "Posterior row count does not match sample count.",
    call. = FALSE
  )
}

rownames(
  posterior_matrix
) <- rownames(
  genus_count_filtered_sample_by_genus
)

colnames(
  posterior_matrix
) <- paste0(
  "Cluster",
  seq_len(
    ncol(
      posterior_matrix
    )
  )
)

hard_assignment <- apply(
  posterior_matrix,
  1,
  which.max
)

maximum_posterior <- apply(
  posterior_matrix,
  1,
  max
)

posterior_entropy <- apply(
  posterior_matrix,
  1,
  function(probabilities) {
    
    probabilities <- probabilities[
      probabilities > 0
    ]
    
    -sum(
      probabilities *
        log(
          probabilities
        )
    )
  }
)

assignment_table <- tibble(
  SampleID = rownames(
    posterior_matrix
  ),
  
  DMM_cluster_number =
    hard_assignment,
  
  DMM_cluster_native =
    paste0(
      "C",
      hard_assignment
    ),
  
  DMM_cluster =
    paste0(
      "C",
      hard_assignment
    ),
  
  Maximum_posterior_probability =
    maximum_posterior,
  
  Posterior_entropy =
    posterior_entropy,
  
  High_confidence_assignment =
    maximum_posterior >=
    0.70
) %>%
  left_join(
    metadata_ordered,
    by = c(
      "SampleID" =
        "SampleID_internal"
    )
  )

write_tsv(
  assignment_table,
  file.path(
    table_dir,
    paste0(
      "DMM_sample_assignments_optimalK",
      optimal_k,
      "_native_7KB.tsv"
    )
  )
)

posterior_table <- as.data.frame(
  posterior_matrix,
  check.names = FALSE
) %>%
  rownames_to_column(
    "SampleID"
  )

write_tsv(
  posterior_table,
  file.path(
    table_dir,
    paste0(
      "DMM_posterior_probabilities_optimalK",
      optimal_k,
      "_7KB.tsv"
    )
  )
)


# ------------------------------------------------------------
# 25. Cluster sizes
# ------------------------------------------------------------

cluster_size_table <- assignment_table %>%
  count(
    DMM_cluster_number,
    DMM_cluster,
    name = "N"
  ) %>%
  mutate(
    Proportion =
      N /
      sum(
        N
      )
  ) %>%
  arrange(
    DMM_cluster_number
  )

write_tsv(
  cluster_size_table,
  file.path(
    table_dir,
    paste0(
      "DMM_cluster_sizes_optimalK",
      optimal_k,
      "_native_7KB.tsv"
    )
  )
)


# ------------------------------------------------------------
# 26. Assignment confidence
# ------------------------------------------------------------

confidence_summary <- assignment_table %>%
  summarise(
    N = n(),
    
    Median_max_posterior =
      median(
        Maximum_posterior_probability
      ),
    
    IQR_max_posterior =
      IQR(
        Maximum_posterior_probability
      ),
    
    N_high_confidence =
      sum(
        High_confidence_assignment
      ),
    
    Proportion_high_confidence =
      mean(
        High_confidence_assignment
      )
  )

write_tsv(
  confidence_summary,
  file.path(
    table_dir,
    paste0(
      "DMM_assignment_confidence_summary_optimalK",
      optimal_k,
      "_native_7KB.tsv"
    )
  )
)


# ------------------------------------------------------------
# 27. Cluster x progression group
# ------------------------------------------------------------

cluster_group_table <- assignment_table %>%
  count(
    Progression5,
    DMM_cluster,
    name = "N"
  ) %>%
  group_by(
    Progression5
  ) %>%
  mutate(
    Proportion_within_group =
      N /
      sum(
        N
      )
  ) %>%
  ungroup()

write_tsv(
  cluster_group_table,
  file.path(
    table_dir,
    paste0(
      "DMM_cluster_by_progression_group_optimalK",
      optimal_k,
      "_native_7KB.tsv"
    )
  )
)


# ------------------------------------------------------------
# 28. Fisher test
# ------------------------------------------------------------

contingency_matrix <- table(
  assignment_table$Progression5,
  assignment_table$DMM_cluster
)

set.seed(
  master_seed
)

fisher_result <- fisher.test(
  contingency_matrix,
  simulate.p.value = TRUE,
  B = 100000
)

fisher_summary <- tibble(
  Test =
    "Fisher exact test with Monte Carlo simulation",
  
  Optimal_K =
    optimal_k,
  
  P_value =
    fisher_result$p.value,
  
  Simulation_replicates =
    100000
)

write_tsv(
  fisher_summary,
  file.path(
    table_dir,
    paste0(
      "DMM_cluster_progression_Fisher_test_K",
      optimal_k,
      "_native_7KB.tsv"
    )
  )
)

write.table(
  contingency_matrix,
  file = file.path(
    table_dir,
    paste0(
      "DMM_cluster_progression_contingency_K",
      optimal_k,
      "_native_7KB.tsv"
    )
  ),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)


# ------------------------------------------------------------
# 29. Component weights
# ------------------------------------------------------------

component_weights <- fitted(
  optimal_model
)

component_weights_df <- as.data.frame(
  component_weights,
  check.names = FALSE
)

write_tsv(
  component_weights_df,
  file.path(
    table_dir,
    paste0(
      "DMM_component_weights_optimalK",
      optimal_k,
      "_raw_7KB.tsv"
    )
  )
)


# ------------------------------------------------------------
# 30. Final console summary
# ------------------------------------------------------------

message("")
message("============================================================")
message("7KB primary DMM completed.")
message("============================================================")

message("")
message("Primary prevalence filter:")
print(
  filter_summary,
  n = Inf,
  width = Inf
)

message("")
message("Best repeat per K:")
print(
  best_repeat_by_k,
  n = Inf,
  width = Inf
)

message("")
message(
  "Optimal K by minimum Laplace: ",
  optimal_k
)

message("")
message("Native component sizes:")
print(
  cluster_size_table,
  n = Inf,
  width = Inf
)

message("")
message("Assignment confidence:")
print(
  confidence_summary,
  n = Inf,
  width = Inf
)

message("")
message("Native component x progression:")
print(
  cluster_group_table,
  n = Inf,
  width = Inf
)

message("")
message("Fisher test:")
print(
  fisher_summary,
  n = Inf,
  width = Inf
)

message("")
message(
  "IMPORTANT: DMM_cluster is still the native component numbering."
)

message(
  "No biological C1/C2/C3 remapping has been applied."
)

message("")
message("Output directory:")
message(output_dir)