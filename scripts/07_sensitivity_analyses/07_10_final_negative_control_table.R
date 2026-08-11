#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_10_final_negative_control_table.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Assemble the final negative-control supplementary table.
############################################################



options(stringsAsFactors = FALSE)
options(width = 220)

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

RAW_DIR <- file.path(PROJECT_ROOT, "raw")

COUNT_FILE <- file.path(
  RAW_DIR,
  "asv_count_raw.tsv"
)

TAX_FILE <- file.path(
  RAW_DIR,
  "taxonomy_from_featuretable.tsv"
)

OUT_DIR <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "05_manuscript_audit/negative_controls_7KB"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

KB_KEEP <- c(
  "kb1",
  "kb3",
  "kb4",
  "kb5",
  "kb6",
  "kb7",
  "kb8"
)

# ==============================================================================
# 2. Read raw ASV count table
# ==============================================================================

count_df <- fread(
  COUNT_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (!"ASV_ID" %in% colnames(count_df)) {
  stop("Expected ASV_ID column was not found in asv_count_raw.tsv.")
}

missing_kb <- setdiff(
  KB_KEEP,
  colnames(count_df)
)

if (length(missing_kb) > 0) {
  stop(
    paste0(
      "Missing retained negative-control columns: ",
      paste(
        missing_kb,
        collapse = ", "
      )
    )
  )
}

count_7kb <- count_df[
  ,
  c(
    "ASV_ID",
    KB_KEEP
  ),
  drop = FALSE
]

count_mat <- as.matrix(
  count_7kb[
    ,
    KB_KEEP,
    drop = FALSE
  ]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- as.character(
  count_7kb$ASV_ID
)

if (any(!is.finite(count_mat))) {
  stop("Non-finite values detected in raw count matrix.")
}

if (any(count_mat < 0)) {
  stop("Negative values detected in raw count matrix.")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop("Raw count matrix contains non-integer values.")
}

count_mat <- round(
  count_mat
)

# ==============================================================================
# 3. Raw reads and detected ASVs
# ==============================================================================

raw_reads <- colSums(
  count_mat
)

detected_asvs <- colSums(
  count_mat > 0
)

# ==============================================================================
# 4. Read taxonomy
# ==============================================================================

tax_df <- fread(
  TAX_FILE,
  data.table = FALSE,
  check.names = FALSE
)

required_tax_cols <- c(
  "ASV_ID",
  "Taxonomy"
)

missing_tax_cols <- setdiff(
  required_tax_cols,
  colnames(tax_df)
)

if (length(missing_tax_cols) > 0) {
  stop(
    paste0(
      "Taxonomy table missing column(s): ",
      paste(
        missing_tax_cols,
        collapse = ", "
      )
    )
  )
}

# ==============================================================================
# 5. Extract genus
# ==============================================================================

extract_genus <- function(x) {
  
  x <- as.character(x)
  
  genus <- rep(
    NA_character_,
    length(x)
  )
  
  for (i in seq_along(x)) {
    
    if (is.na(x[i]) || x[i] == "") {
      next
    }
    
    pieces <- trimws(
      unlist(
        strsplit(
          x[i],
          ";",
          fixed = TRUE
        )
      )
    )
    
    genus_piece <- pieces[
      grepl(
        "^g__",
        pieces
      )
    ]
    
    if (length(genus_piece) > 0) {
      genus[i] <- sub(
        "^g__",
        "",
        genus_piece[1]
      )
    }
  }
  
  genus
}

normalize_genus <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  x <- gsub(
    "_",
    "-",
    x
  )
  
  x <- gsub(
    "\\s+",
    "",
    x
  )
  
  x
}

tax_df$Genus <- extract_genus(
  tax_df$Taxonomy
)

# ==============================================================================
# 6. Match taxonomy to count table
# ==============================================================================

tax_idx <- match(
  rownames(count_mat),
  tax_df$ASV_ID
)

if (any(is.na(tax_idx))) {
  stop(
    paste0(
      "Taxonomy match failed for ",
      sum(is.na(tax_idx)),
      " raw ASVs."
    )
  )
}

genus_vec <- tax_df$Genus[
  tax_idx
]

# ==============================================================================
# 7. UCG-005 ASVs
# ==============================================================================

ucg_mask <- normalize_genus(
  genus_vec
) == "ucg-005"

ucg_mask[
  is.na(ucg_mask)
] <- FALSE

n_ucg_asvs <- sum(
  ucg_mask
)

if (n_ucg_asvs == 0) {
  stop("No UCG-005 ASVs were identified in raw taxonomy.")
}

ucg_reads <- colSums(
  count_mat[
    ucg_mask,
    ,
    drop = FALSE
  ]
)

ucg_percent <- 100 *
  ucg_reads /
  raw_reads

# ==============================================================================
# 8. Build final numeric table
# ==============================================================================

final_table <- data.frame(
  Negative_control = KB_KEEP,
  Raw_reads = as.numeric(
    raw_reads[
      KB_KEEP
    ]
  ),
  Detected_ASVs = as.numeric(
    detected_asvs[
      KB_KEEP
    ]
  ),
  UCG005_reads = as.numeric(
    ucg_reads[
      KB_KEEP
    ]
  ),
  UCG005_percent = as.numeric(
    ucg_percent[
      KB_KEEP
    ]
  ),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 9. Robust validation
# ==============================================================================

expected_table <- data.frame(
  Negative_control = c(
    "kb1",
    "kb3",
    "kb4",
    "kb5",
    "kb6",
    "kb7",
    "kb8"
  ),
  Expected_UCG005_reads = c(
    0,
    5,
    0,
    0,
    6,
    0,
    0
  ),
  stringsAsFactors = FALSE
)

match_idx <- match(
  expected_table$Negative_control,
  final_table$Negative_control
)

if (any(is.na(match_idx))) {
  stop(
    paste0(
      "Could not match expected control(s): ",
      paste(
        expected_table$Negative_control[
          is.na(match_idx)
        ],
        collapse = ", "
      )
    )
  )
}

observed_ucg_reads <- final_table$UCG005_reads[
  match_idx
]

comparison_ok <- observed_ucg_reads ==
  expected_table$Expected_UCG005_reads

if (any(is.na(comparison_ok))) {
  stop(
    "UCG-005 validation produced NA values."
  )
}

if (!all(comparison_ok)) {
  
  validation_table <- data.frame(
    Negative_control = expected_table$Negative_control,
    Expected_UCG005_reads = expected_table$Expected_UCG005_reads,
    Observed_UCG005_reads = observed_ucg_reads,
    stringsAsFactors = FALSE
  )
  
  print(
    validation_table,
    row.names = FALSE
  )
  
  stop(
    "UCG-005 read counts do not reproduce the previously audited result."
  )
}

if (sum(final_table$UCG005_reads) != 11) {
  stop(
    paste0(
      "Expected total UCG-005 reads = 11, observed = ",
      sum(final_table$UCG005_reads)
    )
  )
}

# ==============================================================================
# 10. Manuscript-formatted table
# ==============================================================================

manuscript_table <- data.frame(
  `Negative control` = final_table$Negative_control,
  `Raw reads` = final_table$Raw_reads,
  `Detected ASVs` = final_table$Detected_ASVs,
  `UCG-005 reads` = final_table$UCG005_reads,
  `UCG-005 (%)` = sprintf(
    "%.4f",
    final_table$UCG005_percent
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ==============================================================================
# 11. Console
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("FINAL 7KB NEGATIVE CONTROL TABLE\n")
cat("============================================================\n")

print(
  manuscript_table,
  row.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("UCG-005 SUMMARY\n")
cat("============================================================\n")

cat(
  "Negative controls = 7\n"
)

cat(
  "UCG-005 positive controls = ",
  sum(final_table$UCG005_reads > 0),
  "/7\n",
  sep = ""
)

cat(
  "UCG-005 negative controls = ",
  sum(final_table$UCG005_reads == 0),
  "/7\n",
  sep = ""
)

cat(
  "Total UCG-005 reads = ",
  sum(final_table$UCG005_reads),
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Median UCG-005 abundance = %.4f%%\n",
    median(final_table$UCG005_percent)
  )
)

cat(
  sprintf(
    "Maximum UCG-005 abundance = %.4f%%\n",
    max(final_table$UCG005_percent)
  )
)

# ==============================================================================
# 12. Save tables
# ==============================================================================

fwrite(
  final_table,
  file.path(
    OUT_DIR,
    "Supplementary_Table_negative_controls_7KB_numeric.tsv"
  ),
  sep = "\t"
)

fwrite(
  manuscript_table,
  file.path(
    OUT_DIR,
    "Supplementary_Table_negative_controls_7KB.csv"
  )
)

# ==============================================================================
# 13. Save title and footnote
# ==============================================================================

table_text <- c(
  "Supplementary Table Sx. Sequencing characteristics and UCG-005 detection in negative controls.",
  "",
  paste0(
    "Raw sequencing data from the seven retained negative controls ",
    "(kb1 and kb3-kb8) are shown before contamination filtering. "
  ),
  paste0(
    "UCG-005 was not detected in five of seven controls and was present only ",
    "at trace levels in kb3 and kb6 (5 and 6 reads, respectively). "
  ),
  paste0(
    "Across all seven controls, only 11 UCG-005 reads were detected, ",
    "with a maximum relative abundance of 0.0059%. "
  ),
  "kb2 was excluded from the final analysis."
)

writeLines(
  table_text,
  file.path(
    OUT_DIR,
    "Supplementary_Table_negative_controls_7KB_title_and_footnote.txt"
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