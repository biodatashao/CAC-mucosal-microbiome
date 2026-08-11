#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 05_01_lock_clinical_metadata.R
##
## Module 05 - Clinical association and recurrence (Figure 7, Supplementary Figure 3)
##
## Purpose:
## Audit and lock clinical metadata for Figure 7 and Supplementary Figure 3.
##
## Targets:
##   CAC23
##     - AJCC Stage I / II / III
##     - recurrence / no recurrence
##     - DFS status
##     - DFS time
##
##   Dysplasia17
##     - LGD / HGD
##
## This script performs NO plotting.
############################################################

options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "readxl",
  "dplyr"
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
  library(readxl)
  library(dplyr)
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

clean_root <- file.path(
  rerun_root,
  "00_clean_data"
)

progression_metadata_file <- file.path(
  clean_root,
  "progression127",
  "metadata_7KB_progression127.tsv"
)

ca_clinical_file <- file.path(
  project_root,
  "data",
  "clinical",
  "CA_0629.xlsx"
)

dysplasia_clinical_file <- file.path(
  project_root,
  "data",
  "clinical",
  "dysplasia_0629.xlsx"
)

output_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "00_clinical_metadata_audit"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Input checks
#
# If your Excel files are stored elsewhere, only change the two paths above.
# ==============================================================================

required_files <- c(
  progression_metadata_file,
  ca_clinical_file,
  dysplasia_clinical_file
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
# 3. Helpers
# ==============================================================================

section <- function(title) {
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
}


clean_sample_id <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}


parse_date_safe <- function(x) {
  
  x <- as.character(x)
  x <- trimws(x)
  
  x[
    is.na(x) |
      x == "" |
      x == "NA"
  ] <- NA_character_
  
  result <- as.Date(
    rep(
      NA_character_,
      length(x)
    )
  )
  
  formats <- c(
    "%Y-%m-%d",
    "%Y/%m/%d",
    "%Y.%m.%d",
    "%m/%d/%Y",
    "%d/%m/%Y"
  )
  
  for (fmt in formats) {
    
    unresolved <- is.na(result) &
      !is.na(x)
    
    if (!any(unresolved)) {
      break
    }
    
    parsed <- as.Date(
      x[unresolved],
      format = fmt
    )
    
    result[unresolved] <- parsed
  }
  
  result
}


normalize_recurrence_status <- function(
    dfs_status,
    recurrence_date
) {
  
  dfs_text <- tolower(
    trimws(
      as.character(dfs_status)
    )
  )
  
  recurrence_date_present <- !is.na(
    recurrence_date
  )
  
  output <- rep(
    NA_character_,
    length(dfs_text)
  )
  
  recurrence_tokens <- c(
    "1",
    "recurrence",
    "recurred",
    "event",
    "yes"
  )
  
  no_recurrence_tokens <- c(
    "0",
    "no recurrence",
    "non-recurrence",
    "nonrecurrence",
    "censored",
    "no"
  )
  
  output[
    dfs_text %in%
      recurrence_tokens
  ] <- "Recurrence"
  
  output[
    dfs_text %in%
      no_recurrence_tokens
  ] <- "No recurrence"
  
  output[
    is.na(output) &
      recurrence_date_present
  ] <- "Recurrence"
  
  output
}


# ==============================================================================
# 4. Read fixed progression127 metadata
# ==============================================================================

section(
  "READ FIXED PROGRESSION127 METADATA"
)

progression <- data.table::fread(
  progression_metadata_file,
  data.table = FALSE,
  check.names = FALSE
)

if (!"SampleID" %in% colnames(progression)) {
  stop(
    "progression127 metadata does not contain SampleID.",
    call. = FALSE
  )
}

progression$SampleID <- clean_sample_id(
  progression$SampleID
)


# Identify progression group column
group_candidates <- c(
  "Progression5",
  "Group",
  "group"
)

group_col <- group_candidates[
  group_candidates %in%
    colnames(progression)
][1]

if (is.na(group_col)) {
  stop(
    "Could not identify progression group column.",
    call. = FALSE
  )
}

progression_group <- as.character(
  progression[[group_col]]
)

cac_mask <- progression_group %in%
  c(
    "CA",
    "CAC"
  )

dysplasia_mask <- progression_group %in%
  c(
    "Dysplasia",
    "dysplasia"
  )

cac_samples <- progression$SampleID[
  cac_mask
]

dysplasia_samples <- progression$SampleID[
  dysplasia_mask
]


if (length(cac_samples) != 23) {
  stop(
    paste0(
      "Expected 23 CAC samples in progression127, found ",
      length(cac_samples),
      "."
    ),
    call. = FALSE
  )
}

if (length(dysplasia_samples) != 17) {
  stop(
    paste0(
      "Expected 17 dysplasia samples in progression127, found ",
      length(dysplasia_samples),
      "."
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
# 5. Read CAC clinical Excel
# ==============================================================================

section(
  "READ CAC CLINICAL METADATA"
)

ca_clinical <- readxl::read_excel(
  ca_clinical_file,
  sheet = 1
)

ca_clinical <- as.data.frame(
  ca_clinical,
  check.names = FALSE
)

required_ca_columns <- c(
  "Sample_ID",
  "Stage_group",
  "DFS_status",
  "DFS_time_months",
  "Surgery_date",
  "Recurrence_date"
)

missing_ca_columns <- setdiff(
  required_ca_columns,
  colnames(ca_clinical)
)

if (length(missing_ca_columns) > 0) {
  stop(
    paste0(
      "CAC clinical table missing required column(s): ",
      paste(missing_ca_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

ca_clinical$Sample_ID <- clean_sample_id(
  ca_clinical$Sample_ID
)

if (anyDuplicated(ca_clinical$Sample_ID) > 0) {
  stop(
    "Duplicated Sample_ID detected in CAC clinical table.",
    call. = FALSE
  )
}


# ==============================================================================
# 6. Match fixed CAC23
# ==============================================================================

section(
  "MATCH CAC23 TO CLINICAL METADATA"
)

ca_match_index <- match(
  cac_samples,
  ca_clinical$Sample_ID
)

ca_audit <- data.frame(
  SampleID = cac_samples,
  Clinical_match = !is.na(ca_match_index),
  stringsAsFactors = FALSE
)

matched_ca <- ca_clinical[
  ca_match_index,
  ,
  drop = FALSE
]

ca_audit$Stage_group <- matched_ca$Stage_group
ca_audit$DFS_status_raw <- matched_ca$DFS_status
ca_audit$DFS_time_months_raw <- suppressWarnings(
  as.numeric(
    matched_ca$DFS_time_months
  )
)

ca_audit$Surgery_date <- parse_date_safe(
  matched_ca$Surgery_date
)

ca_audit$Recurrence_date <- parse_date_safe(
  matched_ca$Recurrence_date
)


# ==============================================================================
# 7. Normalize AJCC stage
# ==============================================================================

stage_numeric <- suppressWarnings(
  as.integer(
    ca_audit$Stage_group
  )
)

ca_audit$Tumor_stage_plot <- NA_character_

ca_audit$Tumor_stage_plot[
  stage_numeric == 1
] <- "Stage I"

ca_audit$Tumor_stage_plot[
  stage_numeric == 2
] <- "Stage II"

ca_audit$Tumor_stage_plot[
  stage_numeric == 3
] <- "Stage III"


# ==============================================================================
# 8. Normalize recurrence
# ==============================================================================

ca_audit$Recurrence_plot <- normalize_recurrence_status(
  ca_audit$DFS_status_raw,
  ca_audit$Recurrence_date
)


# ==============================================================================
# 9. Recalculate recurrence-event DFS time where possible
#
# This is an AUDIT value only.
# The original DFS_time_months column is retained separately.
# ==============================================================================

ca_audit$DFS_time_from_dates_months <- NA_real_

valid_date_interval <- !is.na(
  ca_audit$Surgery_date
) &
  !is.na(
    ca_audit$Recurrence_date
  )

ca_audit$DFS_time_from_dates_months[
  valid_date_interval
] <- as.numeric(
  ca_audit$Recurrence_date[
    valid_date_interval
  ] -
    ca_audit$Surgery_date[
      valid_date_interval
    ]
) / 30.4375


ca_audit$DFS_time_difference_months <- NA_real_

valid_compare <- !is.na(
  ca_audit$DFS_time_months_raw
) &
  !is.na(
    ca_audit$DFS_time_from_dates_months
  )

ca_audit$DFS_time_difference_months[
  valid_compare
] <-
  ca_audit$DFS_time_months_raw[
    valid_compare
  ] -
  ca_audit$DFS_time_from_dates_months[
    valid_compare
  ]


# ==============================================================================
# 10. DFS validity
# ==============================================================================

ca_audit$DFS_valid <- !is.na(
  ca_audit$DFS_time_months_raw
) &
  ca_audit$DFS_time_months_raw >= 0 &
  !is.na(
    ca_audit$Recurrence_plot
  )


# ==============================================================================
# 11. CAC audit summaries
# ==============================================================================

stage_table <- table(
  factor(
    ca_audit$Tumor_stage_plot,
    levels = c(
      "Stage I",
      "Stage II",
      "Stage III"
    )
  ),
  useNA = "ifany"
)

recurrence_table <- table(
  factor(
    ca_audit$Recurrence_plot,
    levels = c(
      "No recurrence",
      "Recurrence"
    )
  ),
  useNA = "ifany"
)


cat(
  "Matched CAC clinical metadata: ",
  sum(ca_audit$Clinical_match),
  "/23\n",
  sep = ""
)

cat("\nTumor stage:\n")
print(stage_table)

cat("\nRecurrence:\n")
print(recurrence_table)

cat(
  "\nValid DFS records: ",
  sum(ca_audit$DFS_valid),
  "/23\n",
  sep = ""
)


# ==============================================================================
# 12. Read dysplasia clinical Excel
# ==============================================================================

section(
  "READ DYSPLASIA CLINICAL METADATA"
)

dys_clinical <- readxl::read_excel(
  dysplasia_clinical_file,
  sheet = 1
)

dys_clinical <- as.data.frame(
  dys_clinical,
  check.names = FALSE
)

required_dys_columns <- c(
  "Sample_ID",
  "Dysplasia_grade",
  "Dysplasia_grade_binary"
)

missing_dys_columns <- setdiff(
  required_dys_columns,
  colnames(dys_clinical)
)

if (length(missing_dys_columns) > 0) {
  stop(
    paste0(
      "Dysplasia clinical table missing required column(s): ",
      paste(missing_dys_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

dys_clinical$Sample_ID <- clean_sample_id(
  dys_clinical$Sample_ID
)

if (anyDuplicated(dys_clinical$Sample_ID) > 0) {
  stop(
    "Duplicated Sample_ID detected in dysplasia clinical table.",
    call. = FALSE
  )
}


# ==============================================================================
# 13. Match fixed Dysplasia17
# ==============================================================================

section(
  "MATCH DYSPLASIA17 TO CLINICAL METADATA"
)

dys_match_index <- match(
  dysplasia_samples,
  dys_clinical$Sample_ID
)

dys_audit <- data.frame(
  SampleID = dysplasia_samples,
  Clinical_match = !is.na(dys_match_index),
  stringsAsFactors = FALSE
)

matched_dys <- dys_clinical[
  dys_match_index,
  ,
  drop = FALSE
]

dys_audit$Dysplasia_grade_raw <- matched_dys$Dysplasia_grade
dys_audit$Dysplasia_grade_binary_raw <- matched_dys$Dysplasia_grade_binary


# ==============================================================================
# 14. Normalize dysplasia grade
# ==============================================================================

grade_text <- tolower(
  trimws(
    as.character(
      dys_audit$Dysplasia_grade_binary_raw
    )
  )
)

dys_audit$Dysplasia_grade_plot <- NA_character_

dys_audit$Dysplasia_grade_plot[
  grade_text %in%
    c(
      "high",
      "hgd",
      "high-grade",
      "high grade"
    )
] <- "High-grade dysplasia"

dys_audit$Dysplasia_grade_plot[
  grade_text %in%
    c(
      "low",
      "lgd",
      "low-grade",
      "low grade"
    )
] <- "Low-grade dysplasia"


# fallback to Dysplasia_grade column

raw_grade_text <- tolower(
  trimws(
    as.character(
      dys_audit$Dysplasia_grade_raw
    )
  )
)

dys_audit$Dysplasia_grade_plot[
  is.na(
    dys_audit$Dysplasia_grade_plot
  ) &
    raw_grade_text %in%
    c(
      "hgd",
      "high",
      "high-grade",
      "high grade"
    )
] <- "High-grade dysplasia"

dys_audit$Dysplasia_grade_plot[
  is.na(
    dys_audit$Dysplasia_grade_plot
  ) &
    raw_grade_text %in%
    c(
      "lgd",
      "low",
      "low-grade",
      "low grade"
    )
] <- "Low-grade dysplasia"


grade_table <- table(
  factor(
    dys_audit$Dysplasia_grade_plot,
    levels = c(
      "High-grade dysplasia",
      "Low-grade dysplasia"
    )
  ),
  useNA = "ifany"
)


cat(
  "Matched dysplasia clinical metadata: ",
  sum(dys_audit$Clinical_match),
  "/17\n",
  sep = ""
)

cat("\nDysplasia grade:\n")
print(grade_table)


# ==============================================================================
# 15. Save audit tables
# ==============================================================================

data.table::fwrite(
  ca_audit,
  file.path(
    output_dir,
    "CAC23_clinical_metadata_audit_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

data.table::fwrite(
  dys_audit,
  file.path(
    output_dir,
    "Dysplasia17_clinical_metadata_audit_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


summary_table <- data.frame(
  Metric = c(
    "CAC23_total",
    "CAC23_clinical_matched",
    "CAC23_stage_valid",
    "CAC23_recurrence_valid",
    "CAC23_DFS_valid",
    "Dysplasia17_total",
    "Dysplasia17_clinical_matched",
    "Dysplasia17_grade_valid"
  ),
  
  N = c(
    23,
    sum(ca_audit$Clinical_match),
    sum(!is.na(ca_audit$Tumor_stage_plot)),
    sum(!is.na(ca_audit$Recurrence_plot)),
    sum(ca_audit$DFS_valid),
    17,
    sum(dys_audit$Clinical_match),
    sum(!is.na(dys_audit$Dysplasia_grade_plot))
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  summary_table,
  file.path(
    output_dir,
    "Figure7_clinical_metadata_audit_summary_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 16. Final console summary
# ==============================================================================

section(
  "FINAL AUDIT SUMMARY"
)

print(
  summary_table,
  row.names = FALSE
)

cat("\nCAC stage distribution:\n")
print(stage_table)

cat("\nCAC recurrence distribution:\n")
print(recurrence_table)

cat("\nDysplasia grade distribution:\n")
print(grade_table)


cat("\nCAC samples with missing clinical match:\n")

print(
  ca_audit$SampleID[
    !ca_audit$Clinical_match
  ]
)


cat("\nCAC samples with missing stage:\n")

print(
  ca_audit$SampleID[
    is.na(
      ca_audit$Tumor_stage_plot
    )
  ]
)


cat("\nCAC samples with missing recurrence classification:\n")

print(
  ca_audit$SampleID[
    is.na(
      ca_audit$Recurrence_plot
    )
  ]
)


cat("\nCAC samples with invalid/missing DFS:\n")

print(
  ca_audit$SampleID[
    !ca_audit$DFS_valid
  ]
)


cat("\nDysplasia samples with missing clinical match:\n")

print(
  dys_audit$SampleID[
    !dys_audit$Clinical_match
  ]
)


cat("\nDysplasia samples with missing grade:\n")

print(
  dys_audit$SampleID[
    is.na(
      dys_audit$Dysplasia_grade_plot
    )
  ]
)


cat("\nOutput directory:\n")
cat(
  output_dir,
  "\n"
)