


############################################################
## 01_02_build_paired_CAC23_nonCAC23_cohort.R
##
## Module 01 - Data preparation
##
## Purpose:
##   Build the exact paired CA23 vs nonCA23 cohort.
##
## Inputs:
##   7-KB master clean count / taxonomy / metadata
##
## Outputs:
##   - paired count table
##   - paired taxonomy table
##   - paired metadata table
##   - pair manifest
##   - audit tables
##
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
  library(tibble)
})


############################################################
## 0. Paths
############################################################

PROJECT_DIR <- PROJECT_ROOT


FINAL_ROOT <- file.path(
  PROJECT_DIR,
  "output",
  "analysis"
)


MASTER_DIR <- file.path(
  FINAL_ROOT,
  "00_clean_data",
  "master"
)


OUT_DIR <- file.path(
  FINAL_ROOT,
  "00_clean_data",
  "CA23_nonCA23_paired"
)


AUDIT_DIR <- file.path(
  OUT_DIR,
  "audit"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


COUNT_FILE <- file.path(
  MASTER_DIR,
  "asv_count_clean_7KB_master.tsv"
)


TAX_FILE <- file.path(
  MASTER_DIR,
  "taxonomy_clean_7KB_master.tsv"
)


META_FILE <- file.path(
  MASTER_DIR,
  "metadata_clean_7KB_master.tsv"
)


############################################################
## 1. Check files
############################################################

required_files <- c(
  COUNT_FILE,
  TAX_FILE,
  META_FILE
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
## 2. Exact pair map from the previous paired analysis
############################################################

pair_map <- tibble(
  PairID = sprintf(
    "P%02d",
    1:23
  ),
  
  nonCA = c(
    "nonCA5",
    "nonCA2",
    "nonCA3",
    "nonCA4",
    "nonCA6",
    "nonCA7",
    "nonCA9",
    "nonCA10",
    "nonCA11",
    "nonCA12",
    "nonCA13",
    "nonCA16",
    "nonCA17",
    "nonCA18",
    "nonCA19",
    "nonCA20",
    "nonCA21",
    "nonCA22",
    "nonCA23",
    "nonCA8",
    "nonCA25",
    "nonCA15",
    "nonCA29"
  ),
  
  CA = c(
    "CA1",
    "CA2",
    "CA3",
    "CA4",
    "CA6",
    "CA7",
    "CA9",
    "CA10",
    "CA11",
    "CA12",
    "CA14",
    "CA16",
    "CA17",
    "CA18",
    "CA19",
    "CA20",
    "CA21",
    "CA22",
    "CA23",
    "CA24",
    "CA25",
    "CA27",
    "CA28"
  )
)


############################################################
## 3. Save pair manifest immediately
############################################################

write_csv(
  pair_map,
  file.path(
    OUT_DIR,
    "pair_manifest_CA23_nonCA23.csv"
  )
)


############################################################
## 4. Read master count table
############################################################

count_df <- read_tsv(
  COUNT_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


names(count_df) <- trimws(
  names(count_df)
)


ASV_COL <- names(count_df)[1]


count_mat <- count_df %>%
  column_to_rownames(
    ASV_COL
  ) %>%
  as.matrix()


storage.mode(
  count_mat
) <- "numeric"


############################################################
## 5. Read taxonomy
############################################################

tax_df <- read_tsv(
  TAX_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


names(tax_df) <- trimws(
  names(tax_df)
)


if (!"ASV" %in% names(tax_df)) {
  
  names(tax_df)[1] <- "ASV"
}


tax_df$ASV <- as.character(
  tax_df$ASV
)


############################################################
## 6. Read metadata
############################################################

meta_df <- read_tsv(
  META_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


names(meta_df) <- trimws(
  names(meta_df)
)


if (!"SampleID" %in% names(meta_df)) {
  
  stop(
    "SampleID column not found in master metadata."
  )
}


meta_df$SampleID <- as.character(
  meta_df$SampleID
)


############################################################
## 7. Build long pair table
############################################################

pair_long <- bind_rows(
  
  pair_map %>%
    transmute(
      PairID = PairID,
      SampleID = nonCA,
      Group = "nonCAC"
    ),
  
  pair_map %>%
    transmute(
      PairID = PairID,
      SampleID = CA,
      Group = "CAC"
    )
) %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "nonCAC",
        "CAC"
      )
    )
  )


############################################################
## 8. Verify all 46 samples exist in master count table
############################################################

paired_samples <- pair_long$SampleID


missing_in_count <- setdiff(
  paired_samples,
  colnames(count_mat)
)


missing_in_meta <- setdiff(
  paired_samples,
  meta_df$SampleID
)


if (length(missing_in_count) > 0) {
  
  stop(
    paste0(
      "Paired samples missing from master count table:\n",
      paste(
        missing_in_count,
        collapse = ", "
      )
    )
  )
}


if (length(missing_in_meta) > 0) {
  
  stop(
    paste0(
      "Paired samples missing from master metadata:\n",
      paste(
        missing_in_meta,
        collapse = ", "
      )
    )
  )
}


############################################################
## 9. Verify pair structure
############################################################

pair_check <- pair_long %>%
  count(
    PairID,
    Group,
    name = "N"
  )


bad_pairs <- pair_check %>%
  filter(
    N != 1
  )


if (nrow(bad_pairs) > 0) {
  
  stop(
    "Pair structure is invalid."
  )
}


if (n_distinct(pair_long$PairID) != 23) {
  
  stop(
    "Expected 23 PairIDs."
  )
}


if (nrow(pair_long) != 46) {
  
  stop(
    "Expected 46 paired samples."
  )
}


############################################################
## 10. Extract paired count table
############################################################

count_paired <- count_mat[
  ,
  paired_samples,
  drop = FALSE
]


############################################################
## 11. Remove ASVs absent across all 46 paired samples
############################################################

zero_asvs <- rownames(
  count_paired
)[
  rowSums(count_paired) == 0
]


count_paired <- count_paired[
  rowSums(count_paired) > 0,
  ,
  drop = FALSE
]


############################################################
## 12. Build paired taxonomy
############################################################

tax_paired <- tax_df %>%
  filter(
    ASV %in% rownames(count_paired)
  ) %>%
  arrange(
    match(
      ASV,
      rownames(count_paired)
    )
  )


if (!identical(
  tax_paired$ASV,
  rownames(count_paired)
)) {
  
  stop(
    "Taxonomy order does not match paired count table."
  )
}


############################################################
## 13. Build paired metadata
############################################################

meta_paired <- pair_long %>%
  left_join(
    meta_df,
    by = "SampleID"
  ) %>%
  arrange(
    PairID,
    Group
  )


############################################################
## 14. Reorder count table to match metadata
############################################################

count_paired <- count_paired[
  ,
  meta_paired$SampleID,
  drop = FALSE
]


stopifnot(
  identical(
    colnames(count_paired),
    meta_paired$SampleID
  )
)


############################################################
## 15. Check original explicit labels
############################################################

if ("CA_vs_nonCA_explicit" %in%
    names(meta_paired)) {
  
  explicit_check <- meta_paired %>%
    select(
      PairID,
      SampleID,
      Group,
      CA_vs_nonCA_explicit
    )
  
  write_csv(
    explicit_check,
    file.path(
      AUDIT_DIR,
      "AUDIT_pair_group_vs_original_explicit_label.csv"
    )
  )
  
} else {
  
  explicit_check <- NULL
}


############################################################
## 16. Add paired-analysis flag
############################################################

meta_paired <- meta_paired %>%
  mutate(
    Final_7KB_CA23_nonCA23_paired = TRUE
  )


############################################################
## 17. Save outputs
############################################################

COUNT_OUT <- file.path(
  OUT_DIR,
  "asv_count_7KB_CA23_nonCA23_paired.tsv"
)


TAX_OUT <- file.path(
  OUT_DIR,
  "taxonomy_7KB_CA23_nonCA23_paired.tsv"
)


META_OUT <- file.path(
  OUT_DIR,
  "metadata_7KB_CA23_nonCA23_paired.tsv"
)


count_paired %>%
  as.data.frame(
    check.names = FALSE
  ) %>%
  rownames_to_column(
    "ASV"
  ) %>%
  write_tsv(
    COUNT_OUT
  )


write_tsv(
  tax_paired,
  TAX_OUT
)


write_tsv(
  meta_paired,
  META_OUT
)


############################################################
## 18. Sample depth audit
############################################################

sample_depth <- tibble(
  SampleID = colnames(
    count_paired
  ),
  
  Clean_reads_7KB = colSums(
    count_paired
  )
) %>%
  left_join(
    meta_paired %>%
      select(
        PairID,
        SampleID,
        Group
      ),
    by = "SampleID"
  ) %>%
  arrange(
    PairID,
    Group
  )


write_csv(
  sample_depth,
  file.path(
    AUDIT_DIR,
    "AUDIT_paired_sample_clean_depth_7KB.csv"
  )
)


############################################################
## 19. Pair-level sequencing depth audit
############################################################

pair_depth <- sample_depth %>%
  select(
    PairID,
    Group,
    Clean_reads_7KB
  ) %>%
  tidyr::pivot_wider(
    names_from = Group,
    values_from = Clean_reads_7KB
  ) %>%
  mutate(
    CA_minus_nonCA_reads =
      CAC - nonCAC,
    
    CA_to_nonCA_depth_ratio =
      CAC / nonCAC
  )


write_csv(
  pair_depth,
  file.path(
    AUDIT_DIR,
    "AUDIT_pair_depth_comparison_7KB.csv"
  )
)


############################################################
## 20. Group count audit
############################################################

group_counts <- meta_paired %>%
  count(
    Group,
    name = "N"
  )


write_csv(
  group_counts,
  file.path(
    AUDIT_DIR,
    "AUDIT_group_counts_CA23_nonCA23.csv"
  )
)


############################################################
## 21. Pair completeness audit
############################################################

pair_completeness <- meta_paired %>%
  count(
    PairID,
    name = "N_samples"
  ) %>%
  mutate(
    Complete_pair =
      N_samples == 2
  )


write_csv(
  pair_completeness,
  file.path(
    AUDIT_DIR,
    "AUDIT_pair_completeness.csv"
  )
)


############################################################
## 22. ASV1 audit
############################################################

ASV1_ID <- "ASV1"


asv1_present <- ASV1_ID %in%
  rownames(
    count_paired
  )


if (asv1_present) {
  
  asv1_reads <- count_paired[
    ASV1_ID,
    ,
    drop = TRUE
  ]
  
  
  asv1_audit <- tibble(
    SampleID = names(
      asv1_reads
    ),
    
    ASV1_reads = as.numeric(
      asv1_reads
    )
  ) %>%
    left_join(
      meta_paired %>%
        select(
          PairID,
          SampleID,
          Group
        ),
      by = "SampleID"
    ) %>%
    mutate(
      ASV1_relative_percent =
        100 *
        ASV1_reads /
        sample_depth$Clean_reads_7KB[
          match(
            SampleID,
            sample_depth$SampleID
          )
        ]
    ) %>%
    arrange(
      PairID,
      Group
    )
  
  
  write_csv(
    asv1_audit,
    file.path(
      AUDIT_DIR,
      "AUDIT_ASV1_CA23_nonCA23_paired.csv"
    )
  )
}


############################################################
## 23. Dataset manifest
############################################################

manifest <- tibble(
  Item = c(
    "Analysis_version",
    "Cohort",
    "Pairs",
    "Samples",
    "nonCAC_samples",
    "CAC_samples",
    "Source_count_table",
    "Source_taxonomy",
    "Source_metadata",
    "Output_count_table",
    "Output_taxonomy",
    "Output_metadata"
  ),
  
  Value = c(
    "analysis",
    "Exact paired CA23 vs nonCA23 cohort from previous analysis",
    "23",
    "46",
    "23",
    "23",
    COUNT_FILE,
    TAX_FILE,
    META_FILE,
    COUNT_OUT,
    TAX_OUT,
    META_OUT
  )
)


write_csv(
  manifest,
  file.path(
    OUT_DIR,
    "DATASET_MANIFEST_CA23_nonCA23_paired_7KB.csv"
  )
)


############################################################
## 24. Console report
############################################################

cat("\n")
cat("============================================================\n")
cat("7-KB CA23 vs nonCA23 PAIRED COHORT BUILT\n")
cat("============================================================\n\n")


cat(
  "Paired samples: ",
  ncol(count_paired),
  "\n",
  sep = ""
)


cat(
  "Complete pairs: ",
  sum(
    pair_completeness$Complete_pair
  ),
  "\n",
  sep = ""
)


cat(
  "ASVs retained in paired cohort: ",
  nrow(count_paired),
  "\n\n",
  sep = ""
)


cat("Group counts:\n")

print(
  as.data.frame(
    group_counts
  ),
  row.names = FALSE
)


cat("\nPair map:\n")

print(
  as.data.frame(
    pair_map
  ),
  row.names = FALSE
)


cat("\nASV1 present:\n")
cat(
  asv1_present,
  "\n"
)


cat("\nSequencing depth summary:\n")

depth_summary <- sample_depth %>%
  group_by(
    Group
  ) %>%
  summarise(
    n = n(),
    median_reads = median(
      Clean_reads_7KB
    ),
    min_reads = min(
      Clean_reads_7KB
    ),
    max_reads = max(
      Clean_reads_7KB
    ),
    .groups = "drop"
  )


print(
  as.data.frame(
    depth_summary
  ),
  row.names = FALSE
)


cat("\nOutput files:\n")

cat(
  COUNT_OUT,
  "\n"
)

cat(
  TAX_OUT,
  "\n"
)

cat(
  META_OUT,
  "\n"
)


cat("\nDone.\n")