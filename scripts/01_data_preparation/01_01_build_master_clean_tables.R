


############################################################
## 01_01_build_master_clean_tables.R
##
## Module 01 - Data preparation
##
## FINAL 7-KB MASTER CLEAN DATASET
##
## Purpose
## 1. Use the completed 7-KB decontamination result.
##
## 2. Build ONE master clean ASV table containing ALL
##    biological samples.
##
## 3. Construct the progression127 cohort.
##
## 4. Also generate all explicit CA/nonCA samples as a
##    candidate pool. 
############################################################


rm(list = ls())

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})


############################################################
## 0. Project paths
############################################################

PROJECT_DIR <- PROJECT_ROOT

RAW_DIR <- file.path(
  PROJECT_DIR,
  "raw"
)

OLD_127_DIR <- file.path(
  PROJECT_DIR,
  "output",
  "decontamination_8controls",
  "clean_tables_previous127_v0.5_step2_retention_flagged"
)

SENS_DIR <- file.path(
  PROJECT_DIR,
  "output",
  "decontamination_7controls_final",
  "decontamination"
)

FINAL_ROOT <- file.path(
  PROJECT_DIR,
  "output",
  "analysis"
)

CLEAN_ROOT <- file.path(
  FINAL_ROOT,
  "00_clean_data"
)

MASTER_DIR <- file.path(
  CLEAN_ROOT,
  "master"
)

PROGRESSION_DIR <- file.path(
  CLEAN_ROOT,
  "progression127"
)

EXPLICIT_DIR <- file.path(
  CLEAN_ROOT,
  "explicit_CA_nonCA_candidate"
)

AUDIT_DIR <- file.path(
  CLEAN_ROOT,
  "audit"
)


for (d in c(
  MASTER_DIR,
  PROGRESSION_DIR,
  EXPLICIT_DIR,
  AUDIT_DIR
)) {
  
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


############################################################
## 1. Input files
############################################################

RAW_COUNT_FILE <- file.path(
  RAW_DIR,
  "asv_count_raw.tsv"
)

RAW_META_FILE <- file.path(
  RAW_DIR,
  "metadata_main.tsv"
)

RAW_TAX_FILE <- file.path(
  RAW_DIR,
  "taxonomy_from_featuretable.tsv"
)


OLD_127_COUNT_FILE <- file.path(
  OLD_127_DIR,
  "asv_count_clean_v0.5_step2_previous127.tsv"
)

OLD_127_META_FILE <- file.path(
  OLD_127_DIR,
  "metadata_clean_v0.5_step2_previous127.tsv"
)


REMOVED_STEP0_FILE <- file.path(
  SENS_DIR,
  "removed_ASVs_step0_non_target.tsv"
)

REMOVED_STEP1_FILE <- file.path(
  SENS_DIR,
  "removed_ASVs_step1_decontam.tsv"
)

REMOVED_STEP2_FILE <- file.path(
  SENS_DIR,
  "removed_ASVs_step2_genus_blacklist.tsv"
)

DECONTAM_7KB_FILE <- file.path(
  SENS_DIR,
  "decontam_ASV_results.tsv"
)

RETENTION_7KB_FILE <- file.path(
  SENS_DIR,
  "sample_retention_fraction.tsv"
)


required_files <- c(
  RAW_COUNT_FILE,
  RAW_META_FILE,
  RAW_TAX_FILE,
  OLD_127_COUNT_FILE,
  OLD_127_META_FILE,
  REMOVED_STEP0_FILE,
  REMOVED_STEP1_FILE,
  REMOVED_STEP2_FILE
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "Missing required files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


############################################################
## 2. Helper functions
############################################################

detect_asv_column <- function(df) {
  
  candidates <- c(
    "ASV",
    "ASV_ID",
    "Feature.ID",
    "FeatureID",
    "#OTU ID"
  )
  
  hit <- candidates[
    candidates %in% names(df)
  ]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  names(df)[1]
}


logical_from_text <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  tolower(
    trimws(
      as.character(x)
    )
  ) %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y",
    "blank"
  )
}


clean_character <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x[
    x %in% c(
      "",
      "NA",
      "NaN",
      "NULL"
    )
  ] <- NA_character_
  
  x
}


write_count_table <- function(mat, path) {
  
  mat %>%
    as.data.frame(
      check.names = FALSE
    ) %>%
    rownames_to_column(
      "ASV"
    ) %>%
    write_tsv(
      path
    )
}


############################################################
## 3. Read raw count table
############################################################

count_raw <- read_tsv(
  RAW_COUNT_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


names(count_raw) <- trimws(
  names(count_raw)
)


ASV_COL_COUNT <- names(count_raw)[1]


count_mat_raw <- count_raw %>%
  column_to_rownames(
    ASV_COL_COUNT
  ) %>%
  as.matrix()


storage.mode(
  count_mat_raw
) <- "numeric"


if (anyNA(count_mat_raw)) {
  
  stop(
    "NA values detected in raw count matrix."
  )
}


non_integer_n <- sum(
  abs(
    count_mat_raw -
      round(count_mat_raw)
  ) > 1e-8
)


if (non_integer_n > 0) {
  
  stop(
    paste0(
      "Raw count matrix contains ",
      non_integer_n,
      " non-integer values."
    )
  )
}


count_mat_raw <- round(
  count_mat_raw
)


storage.mode(
  count_mat_raw
) <- "integer"


############################################################
## 4. Read raw metadata
############################################################

meta_raw <- read_tsv(
  RAW_META_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


names(meta_raw) <- trimws(
  names(meta_raw)
)


if (!"SampleID" %in% names(meta_raw)) {
  
  stop(
    "SampleID column not found in raw metadata."
  )
}


meta_raw$SampleID <- clean_character(
  meta_raw$SampleID
)


if (anyDuplicated(meta_raw$SampleID) > 0) {
  
  dup <- unique(
    meta_raw$SampleID[
      duplicated(meta_raw$SampleID)
    ]
  )
  
  stop(
    paste0(
      "Duplicated SampleIDs in raw metadata:\n",
      paste(
        dup,
        collapse = ", "
      )
    )
  )
}


############################################################
## 5. Identify blanks
############################################################

if ("IsBlank" %in% names(meta_raw)) {
  
  meta_raw <- meta_raw %>%
    mutate(
      IsBlank_master = logical_from_text(
        IsBlank
      )
    )
  
} else {
  
  meta_raw <- meta_raw %>%
    mutate(
      IsBlank_master = str_detect(
        tolower(SampleID),
        "^kb[0-9]+$"
      )
    )
}


cat("\n============================================\n")
cat("RAW SAMPLE INVENTORY\n")
cat("============================================\n\n")


cat(
  "Raw metadata samples: ",
  nrow(meta_raw),
  "\n",
  sep = ""
)


cat(
  "Blank controls detected: ",
  sum(
    meta_raw$IsBlank_master,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Biological samples detected: ",
  sum(
    !meta_raw$IsBlank_master,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat("\nBlank IDs:\n")

print(
  meta_raw %>%
    filter(
      IsBlank_master
    ) %>%
    select(
      SampleID
    )
)


############################################################
## 6. Match metadata to count table
############################################################

missing_meta_in_count <- setdiff(
  meta_raw$SampleID,
  colnames(count_mat_raw)
)


if (length(missing_meta_in_count) > 0) {
  
  warning(
    paste0(
      "Metadata samples missing from raw count table:\n",
      paste(
        missing_meta_in_count,
        collapse = ", "
      )
    )
  )
}


common_samples <- intersect(
  meta_raw$SampleID,
  colnames(count_mat_raw)
)


meta_raw <- meta_raw %>%
  filter(
    SampleID %in% common_samples
  ) %>%
  arrange(
    match(
      SampleID,
      common_samples
    )
  )


count_mat_raw <- count_mat_raw[
  ,
  common_samples,
  drop = FALSE
]


stopifnot(
  identical(
    meta_raw$SampleID,
    colnames(count_mat_raw)
  )
)


############################################################
## 7. Read taxonomy
############################################################

tax_raw <- read_tsv(
  RAW_TAX_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


names(tax_raw) <- trimws(
  names(tax_raw)
)


ASV_COL_TAX <- names(tax_raw)[1]


tax_raw <- tax_raw %>%
  rename(
    ASV = all_of(
      ASV_COL_TAX
    )
  )


tax_raw$ASV <- clean_character(
  tax_raw$ASV
)


if (anyDuplicated(tax_raw$ASV) > 0) {
  
  stop(
    "Duplicated ASVs detected in raw taxonomy."
  )
}


############################################################
## 8. Read 7-KB removal sets
############################################################

removed_step0 <- read_tsv(
  REMOVED_STEP0_FILE,
  show_col_types = FALSE,
  progress = FALSE
)


removed_step1 <- read_tsv(
  REMOVED_STEP1_FILE,
  show_col_types = FALSE,
  progress = FALSE
)


removed_step2 <- read_tsv(
  REMOVED_STEP2_FILE,
  show_col_types = FALSE,
  progress = FALSE
)


step0_col <- detect_asv_column(
  removed_step0
)


step1_col <- detect_asv_column(
  removed_step1
)


step2_col <- detect_asv_column(
  removed_step2
)


removed_step0_asvs <- unique(
  clean_character(
    removed_step0[[step0_col]]
  )
)


removed_step1_asvs <- unique(
  clean_character(
    removed_step1[[step1_col]]
  )
)


removed_step2_asvs <- unique(
  clean_character(
    removed_step2[[step2_col]]
  )
)


removed_step0_asvs <- na.omit(
  removed_step0_asvs
)


removed_step1_asvs <- na.omit(
  removed_step1_asvs
)


removed_step2_asvs <- na.omit(
  removed_step2_asvs
)


removed_all_asvs <- unique(
  c(
    removed_step0_asvs,
    removed_step1_asvs,
    removed_step2_asvs
  )
)


cat("\n============================================\n")
cat("7-KB ASV CLEANING\n")
cat("============================================\n\n")


cat(
  "Step 0 removed ASVs: ",
  length(removed_step0_asvs),
  "\n",
  sep = ""
)


cat(
  "Step 1 removed ASVs: ",
  length(removed_step1_asvs),
  "\n",
  sep = ""
)


cat(
  "Step 2 removed ASVs: ",
  length(removed_step2_asvs),
  "\n",
  sep = ""
)


cat(
  "Unique ASVs removed overall: ",
  length(removed_all_asvs),
  "\n",
  sep = ""
)


############################################################
## 9. Build MASTER biological sample set
##
## IMPORTANT:
## - remove ALL blanks from downstream master
## - therefore kb1-kb8 are not biological samples
## - kb2 has already been excluded from the decontam model
## - NO retention-based biological sample removal here
############################################################

bio_samples <- meta_raw %>%
  filter(
    !IsBlank_master
  ) %>%
  pull(
    SampleID
  )


count_master <- count_mat_raw[
  !rownames(count_mat_raw) %in%
    removed_all_asvs,
  bio_samples,
  drop = FALSE
]


############################################################
## 10. Remove ASVs absent from ALL biological samples
############################################################

zero_master_asvs <- rownames(
  count_master
)[
  rowSums(count_master) == 0
]


count_master <- count_master[
  rowSums(count_master) > 0,
  ,
  drop = FALSE
]


############################################################
## 11. Master taxonomy
############################################################

missing_tax_asvs <- setdiff(
  rownames(count_master),
  tax_raw$ASV
)


if (length(missing_tax_asvs) > 0) {
  
  stop(
    paste0(
      "Taxonomy missing for ",
      length(missing_tax_asvs),
      " clean ASVs."
    )
  )
}


tax_master <- tax_raw %>%
  filter(
    ASV %in% rownames(count_master)
  ) %>%
  arrange(
    match(
      ASV,
      rownames(count_master)
    )
  )


stopifnot(
  identical(
    tax_master$ASV,
    rownames(count_master)
  )
)


############################################################
## 12. Master metadata
############################################################

meta_master <- meta_raw %>%
  filter(
    SampleID %in% bio_samples
  ) %>%
  arrange(
    match(
      SampleID,
      bio_samples
    )
  )


############################################################
## 13. Attach 7-KB retention as AUDIT information
##
## This does NOT determine sample inclusion.
############################################################

if (file.exists(
  RETENTION_7KB_FILE
)) {
  
  retention7 <- read_tsv(
    RETENTION_7KB_FILE,
    show_col_types = FALSE,
    progress = FALSE
  )
  
  
  retention7 <- retention7 %>%
    select(
      any_of(
        c(
          "SampleID",
          "Raw_reads",
          "Reads_after_step2",
          "Retention_fraction",
          "High_contamination_risk",
          "Removal_reason"
        )
      )
    )
  
  
  old_names <- names(
    retention7
  )
  
  
  names(retention7) <- ifelse(
    old_names == "SampleID",
    "SampleID",
    paste0(
      old_names,
      "_7KB_audit"
    )
  )
  
  
  meta_master <- meta_master %>%
    left_join(
      retention7,
      by = "SampleID"
    )
}


meta_master <- meta_master %>%
  mutate(
    Master_7KB_biological_sample = TRUE
  )


stopifnot(
  identical(
    colnames(count_master),
    meta_master$SampleID
  )
)


############################################################
## 14. Write MASTER clean files
############################################################

MASTER_COUNT_OUT <- file.path(
  MASTER_DIR,
  "asv_count_clean_7KB_master.tsv"
)


MASTER_TAX_OUT <- file.path(
  MASTER_DIR,
  "taxonomy_clean_7KB_master.tsv"
)


MASTER_META_OUT <- file.path(
  MASTER_DIR,
  "metadata_clean_7KB_master.tsv"
)


write_count_table(
  count_master,
  MASTER_COUNT_OUT
)


write_tsv(
  tax_master,
  MASTER_TAX_OUT
)


write_tsv(
  meta_master,
  MASTER_META_OUT
)


############################################################
## 15. Reconstruct EXACT original progression127 cohort
############################################################

meta_old127 <- read_tsv(
  OLD_127_META_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


if (!"SampleID" %in% names(meta_old127)) {
  
  stop(
    "SampleID missing from old previous127 metadata."
  )
}


meta_old127$SampleID <- clean_character(
  meta_old127$SampleID
)


progression_samples <- meta_old127$SampleID


if (length(progression_samples) != 127) {
  
  stop(
    paste0(
      "Old previous127 metadata contains ",
      length(progression_samples),
      " samples instead of 127."
    )
  )
}


missing_progression <- setdiff(
  progression_samples,
  colnames(count_master)
)


if (length(missing_progression) > 0) {
  
  stop(
    paste0(
      "The following progression127 samples are missing ",
      "from the 7-KB master:\n",
      paste(
        missing_progression,
        collapse = ", "
      )
    )
  )
}


count_progression <- count_master[
  ,
  progression_samples,
  drop = FALSE
]


count_progression <- count_progression[
  rowSums(count_progression) > 0,
  ,
  drop = FALSE
]


############################################################
## 16. Preserve original previous127 metadata definitions
##
## We use the OLD cohort metadata because it represents the
## exact cohort/grouping used in the manuscript.
##
## Then attach selected fresh 7KB audit variables.
############################################################

meta_progression <- meta_old127 %>%
  arrange(
    match(
      SampleID,
      progression_samples
    )
  )


if ("Retention_fraction_7KB_audit" %in%
    names(meta_master)) {
  
  fresh_audit <- meta_master %>%
    select(
      SampleID,
      ends_with(
        "_7KB_audit"
      )
    )
  
  
  meta_progression <- meta_progression %>%
    select(
      -any_of(
        names(fresh_audit)[
          names(fresh_audit) != "SampleID"
        ]
      )
    ) %>%
    left_join(
      fresh_audit,
      by = "SampleID"
    )
}


meta_progression <- meta_progression %>%
  mutate(
    Final_7KB_progression127 = TRUE
  )


stopifnot(
  identical(
    colnames(count_progression),
    meta_progression$SampleID
  )
)


############################################################
## 17. Progression taxonomy
############################################################

tax_progression <- tax_master %>%
  filter(
    ASV %in% rownames(count_progression)
  ) %>%
  arrange(
    match(
      ASV,
      rownames(count_progression)
    )
  )


############################################################
## 18. Write progression127 files
############################################################

write_count_table(
  count_progression,
  file.path(
    PROGRESSION_DIR,
    "asv_count_7KB_progression127.tsv"
  )
)


write_tsv(
  tax_progression,
  file.path(
    PROGRESSION_DIR,
    "taxonomy_7KB_progression127.tsv"
  )
)


write_tsv(
  meta_progression,
  file.path(
    PROGRESSION_DIR,
    "metadata_7KB_progression127.tsv"
  )
)


############################################################
## 19. Verify Progression5 counts
############################################################

if (!"Progression5" %in%
    names(meta_progression)) {
  
  stop(
    "Progression5 missing from previous127 metadata."
  )
}


progression_group_counts <- meta_progression %>%
  count(
    Progression5,
    name = "N"
  )


write_csv(
  progression_group_counts,
  file.path(
    AUDIT_DIR,
    "AUDIT_progression127_group_counts.csv"
  )
)


############################################################
## 20. Identify explicit CA/nonCA candidate pool
##
## IMPORTANT:
## This is NOT yet necessarily the final paired 23+23 set.
############################################################

explicit_col_candidates <- c(
  "CA_vs_nonCA_explicit",
  "CA_vs_nonCA",
  "CA_nonCA",
  "CAC_vs_nonCAC"
)


explicit_col <- explicit_col_candidates[
  explicit_col_candidates %in%
    names(meta_master)
]


if (length(explicit_col) == 0) {
  
  warning(
    paste0(
      "No explicit CA/nonCA grouping column detected.\n",
      "Candidate CA/nonCA files will not be generated."
    )
  )
  
  
  explicit_meta <- NULL
  
} else {
  
  explicit_col <- explicit_col[1]
  
  
  cat(
    "\nExplicit CA/nonCA grouping column detected: ",
    explicit_col,
    "\n",
    sep = ""
  )
  
  
  explicit_values <- clean_character(
    meta_master[[explicit_col]]
  )
  
  
  cat("\nRaw values in explicit grouping column:\n")
  
  print(
    table(
      explicit_values,
      useNA = "ifany"
    )
  )
  
  
  explicit_meta <- meta_master %>%
    mutate(
      Explicit_group_raw =
        clean_character(
          .data[[explicit_col]]
        )
    ) %>%
    filter(
      !is.na(
        Explicit_group_raw
      )
    )
  
  
  ##########################################################
  ## 20A. Standardize CA/nonCA labels cautiously
  ##########################################################
  
  explicit_meta <- explicit_meta %>%
    mutate(
      Explicit_CA_nonCA = case_when(
        
        tolower(
          Explicit_group_raw
        ) %in%
          c(
            "ca",
            "cac",
            "cancer"
          ) ~
          "CAC",
        
        tolower(
          Explicit_group_raw
        ) %in%
          c(
            "nonca",
            "noncac",
            "non-cac",
            "non-ca"
          ) ~
          "nonCAC",
        
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(
        Explicit_CA_nonCA
      )
    )
  
  
  explicit_counts <- explicit_meta %>%
    count(
      Explicit_CA_nonCA,
      name = "N"
    )
  
  
  write_csv(
    explicit_counts,
    file.path(
      AUDIT_DIR,
      "AUDIT_explicit_CA_nonCA_candidate_counts.csv"
    )
  )
  
  
  explicit_samples <- explicit_meta$SampleID
  
  
  count_explicit <- count_master[
    ,
    explicit_samples,
    drop = FALSE
  ]
  
  
  count_explicit <- count_explicit[
    rowSums(
      count_explicit
    ) > 0,
    ,
    drop = FALSE
  ]
  
  
  tax_explicit <- tax_master %>%
    filter(
      ASV %in%
        rownames(
          count_explicit
        )
    ) %>%
    arrange(
      match(
        ASV,
        rownames(
          count_explicit
        )
      )
    )
  
  
  explicit_meta <- explicit_meta %>%
    arrange(
      match(
        SampleID,
        explicit_samples
      )
    )
  
  
  stopifnot(
    identical(
      colnames(count_explicit),
      explicit_meta$SampleID
    )
  )
  
  
  write_count_table(
    count_explicit,
    file.path(
      EXPLICIT_DIR,
      "asv_count_7KB_explicit_CA_nonCA_candidate.tsv"
    )
  )
  
  
  write_tsv(
    tax_explicit,
    file.path(
      EXPLICIT_DIR,
      "taxonomy_7KB_explicit_CA_nonCA_candidate.tsv"
    )
  )
  
  
  write_tsv(
    explicit_meta,
    file.path(
      EXPLICIT_DIR,
      "metadata_7KB_explicit_CA_nonCA_candidate.tsv"
    )
  )
}


############################################################
## 21. ASV1 audit in master
############################################################

ASV1_ID <- "ASV1"


asv1_master_audit <- tibble(
  Metric = c(
    "ASV1_present_in_7KB_master",
    "ASV1_present_in_progression127"
  ),
  
  Value = c(
    ASV1_ID %in%
      rownames(count_master),
    
    ASV1_ID %in%
      rownames(count_progression)
  )
)


write_csv(
  asv1_master_audit,
  file.path(
    AUDIT_DIR,
    "AUDIT_ASV1_presence.csv"
  )
)


############################################################
## 22. Compare old 8-KB vs new 7-KB progression127
############################################################

old127_count <- read_tsv(
  OLD_127_COUNT_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


OLD_ASV_COL <- names(
  old127_count
)[1]


old127_mat <- old127_count %>%
  column_to_rownames(
    OLD_ASV_COL
  ) %>%
  as.matrix()


storage.mode(
  old127_mat
) <- "numeric"


old_asvs <- rownames(
  old127_mat
)


new_asvs <- rownames(
  count_progression
)


shared_asvs <- intersect(
  old_asvs,
  new_asvs
)


old_only <- setdiff(
  old_asvs,
  new_asvs
)


new_only <- setdiff(
  new_asvs,
  old_asvs
)


asv_set_compare <- tibble(
  Metric = c(
    "Old_8KB_ASVs_progression127",
    "New_7KB_ASVs_progression127",
    "Shared_ASVs",
    "Old_only_ASVs",
    "New_only_ASVs"
  ),
  
  N = c(
    length(old_asvs),
    length(new_asvs),
    length(shared_asvs),
    length(old_only),
    length(new_only)
  )
)


write_csv(
  asv_set_compare,
  file.path(
    AUDIT_DIR,
    "AUDIT_8KB_vs_7KB_progression127_ASV_set.csv"
  )
)


write_csv(
  tibble(
    ASV = old_only
  ),
  file.path(
    AUDIT_DIR,
    "AUDIT_ASVs_old8KB_only_progression127.csv"
  )
)


write_csv(
  tibble(
    ASV = new_only
  ),
  file.path(
    AUDIT_DIR,
    "AUDIT_ASVs_new7KB_only_progression127.csv"
  )
)


############################################################
## 23. Sample-by-sample sequencing depth comparison
############################################################

shared_samples <- intersect(
  colnames(old127_mat),
  colnames(count_progression)
)


old_depth <- colSums(
  old127_mat[
    ,
    shared_samples,
    drop = FALSE
  ]
)


new_depth <- colSums(
  count_progression[
    ,
    shared_samples,
    drop = FALSE
  ]
)


depth_compare <- tibble(
  SampleID = shared_samples,
  
  Reads_old_8KB =
    as.numeric(
      old_depth[
        shared_samples
      ]
    ),
  
  Reads_new_7KB =
    as.numeric(
      new_depth[
        shared_samples
      ]
    )
) %>%
  mutate(
    Difference_reads =
      Reads_new_7KB -
      Reads_old_8KB,
    
    Difference_percent_old =
      100 *
      Difference_reads /
      Reads_old_8KB,
    
    Absolute_difference_reads =
      abs(
        Difference_reads
      )
  ) %>%
  arrange(
    desc(
      Absolute_difference_reads
    )
  )


write_csv(
  depth_compare,
  file.path(
    AUDIT_DIR,
    "AUDIT_sample_depth_8KB_vs_7KB_progression127.csv"
  )
)


############################################################
## 24. Dataset manifest
############################################################

manifest <- tibble(
  Item = c(
    "Analysis_version",
    "Excluded_negative_control",
    "Negative_controls_used_for_decontam",
    "Master_sample_definition",
    "Progression127_definition",
    "Explicit_CA_nonCA_definition",
    "Retention_rule_for_master",
    "Master_count",
    "Master_taxonomy",
    "Master_metadata"
  ),
  
  Value = c(
    "analysis",
    "kb2",
    "kb1,kb3,kb4,kb5,kb6,kb7,kb8",
    "All nonblank biological samples after 7-KB ASV decontamination",
    "Exact SampleIDs from original previous127 metadata",
    "Candidate pool only; exact paired 23+23 cohort to be locked separately",
    "Retention metrics retained as audit variables; not used for automatic master sample exclusion",
    MASTER_COUNT_OUT,
    MASTER_TAX_OUT,
    MASTER_META_OUT
  )
)


write_csv(
  manifest,
  file.path(
    AUDIT_DIR,
    "DATASET_MANIFEST_7KB_master.csv"
  )
)


############################################################
## 25. Console report
############################################################

cat("\n")
cat("============================================================\n")
cat("FINAL 7-KB MASTER CLEAN DATASET COMPLETE\n")
cat("============================================================\n\n")


cat(
  "MASTER biological samples: ",
  ncol(count_master),
  "\n",
  sep = ""
)


cat(
  "MASTER ASVs: ",
  nrow(count_master),
  "\n\n",
  sep = ""
)


cat(
  "Progression cohort samples: ",
  ncol(count_progression),
  "\n",
  sep = ""
)


cat(
  "Progression cohort ASVs: ",
  nrow(count_progression),
  "\n\n",
  sep = ""
)


cat("Progression5 group counts:\n")

print(
  as.data.frame(
    progression_group_counts
  ),
  row.names = FALSE
)


cat("\n8-KB vs 7-KB progression127 ASV comparison:\n")

print(
  as.data.frame(
    asv_set_compare
  ),
  row.names = FALSE
)


cat("\nASV1:\n")

print(
  as.data.frame(
    asv1_master_audit
  ),
  row.names = FALSE
)


if (!is.null(
  explicit_meta
)) {
  
  cat("\nExplicit CA/nonCA CANDIDATE pool:\n")
  
  print(
    as.data.frame(
      explicit_counts
    ),
    row.names = FALSE
  )
  
  
  cat(
    "\nCandidate total: ",
    nrow(explicit_meta),
    "\n",
    sep = ""
  )
  
  
  if (nrow(explicit_meta) != 46) {
    
    cat(
      "\nIMPORTANT:\n",
      "Candidate pool is NOT assumed to be the final paired 23+23 cohort.\n",
      "We will recover the exact paired sample definition from the old 07 script next.\n",
      sep = ""
    )
  }
}


cat("\nLargest sequencing-depth changes in progression127:\n")

print(
  head(
    as.data.frame(
      depth_compare
    ),
    15
  ),
  row.names = FALSE
)


cat("\nMASTER files:\n")

cat(
  MASTER_COUNT_OUT,
  "\n"
)

cat(
  MASTER_TAX_OUT,
  "\n"
)

cat(
  MASTER_META_OUT,
  "\n"
)


cat("\nDone.\n")