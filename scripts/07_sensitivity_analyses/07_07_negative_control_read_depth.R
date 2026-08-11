#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_07_negative_control_read_depth.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Read depth and sequencing source of each negative control.
############################################################



options(stringsAsFactors = FALSE)
options(width = 240)

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

BASE_DIR <- paste0(
  file.path(PROJECT_ROOT, "output/"),
  "decontamination_7controls_final"
)

DECONTAM_DIR <- file.path(
  BASE_DIR,
  "decontamination"
)

RETENTION_FILE <- file.path(
  DECONTAM_DIR,
  "sample_retention_fraction.tsv"
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

KB_EXCLUDED <- "kb2"


# ==============================================================================
# 2. Check authoritative without-kb2 retention table
# ==============================================================================

if (!file.exists(RETENTION_FILE)) {
  stop(
    paste0(
      "Retention file not found:\n",
      RETENTION_FILE
    )
  )
}

ret <- fread(
  RETENTION_FILE,
  data.table = FALSE,
  check.names = FALSE
)

required_cols <- c(
  "SampleID",
  "Raw_reads",
  "Reads_after_step2",
  "Retention_fraction"
)

missing_cols <- setdiff(
  required_cols,
  colnames(ret)
)

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Retention table missing column(s): ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}

ret$SampleID <- tolower(
  trimws(
    as.character(
      ret$SampleID
    )
  )
)


# ==============================================================================
# 3. Explicit kb2 audit
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("7KB NEGATIVE CONTROL READ-DEPTH AUDIT\n")
cat("============================================================\n")

cat(
  "Source:\n",
  RETENTION_FILE,
  "\n\n",
  sep = ""
)

if (KB_EXCLUDED %in% ret$SampleID) {
  stop(
    "ERROR: kb2 is present in the WITHOUT-kb2 retention table."
  )
}

cat("kb2 present in retention table: NO\n")


# ==============================================================================
# 4. Extract seven retained negative controls
# ==============================================================================

kb <- ret[
  ret$SampleID %in% KB_KEEP,
  ,
  drop = FALSE
]

missing_kb <- setdiff(
  KB_KEEP,
  kb$SampleID
)

if (length(missing_kb) > 0) {
  stop(
    paste0(
      "Missing retained control(s): ",
      paste(
        missing_kb,
        collapse = ", "
      )
    )
  )
}

if (nrow(kb) != 7) {
  stop(
    paste0(
      "Expected exactly 7 negative controls, found ",
      nrow(kb),
      "."
    )
  )
}

kb$SampleID <- factor(
  kb$SampleID,
  levels = KB_KEEP
)

kb <- kb[
  order(kb$SampleID),
  ,
  drop = FALSE
]

kb$SampleID <- as.character(
  kb$SampleID
)


# ==============================================================================
# 5. Read-depth summary
# ==============================================================================

kb_summary <- data.frame(
  Metric = c(
    "Number_of_negative_controls",
    "Raw_reads_min",
    "Raw_reads_median",
    "Raw_reads_max",
    "Post_step2_reads_min",
    "Post_step2_reads_median",
    "Post_step2_reads_max",
    "Retention_fraction_min",
    "Retention_fraction_median",
    "Retention_fraction_max"
  ),
  Value = c(
    nrow(kb),
    min(kb$Raw_reads),
    median(kb$Raw_reads),
    max(kb$Raw_reads),
    min(kb$Reads_after_step2),
    median(kb$Reads_after_step2),
    max(kb$Reads_after_step2),
    min(kb$Retention_fraction),
    median(kb$Retention_fraction),
    max(kb$Retention_fraction)
  ),
  stringsAsFactors = FALSE
)


# ==============================================================================
# 6. Console report
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("SEVEN NEGATIVE CONTROLS\n")
cat("============================================================\n")

display_cols <- intersect(
  c(
    "SampleID",
    "Raw_reads",
    "Reads_after_step2",
    "Retention_fraction",
    "IsBlank_final",
    "High_contamination_risk",
    "Remove_sample_final",
    "Removal_reason"
  ),
  colnames(kb)
)

print(
  kb[
    ,
    display_cols,
    drop = FALSE
  ],
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("READ-DEPTH SUMMARY\n")
cat("============================================================\n")

print(
  kb_summary,
  row.names = FALSE
)


# ==============================================================================
# 7. Locate count tables from the SAME without-kb2 branch
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("COUNT TABLE CANDIDATES FROM WITHOUT-kb2 BRANCH\n")
cat("============================================================\n")

all_files <- list.files(
  BASE_DIR,
  recursive = TRUE,
  full.names = TRUE
)

count_candidates <- all_files[
  grepl(
    "count.*step|step.*count|asv.*count|count.*asv",
    basename(all_files),
    ignore.case = TRUE
  ) &
    grepl(
      "\\.(tsv|csv)$",
      all_files,
      ignore.case = TRUE
    )
]

if (length(count_candidates) == 0) {
  
  cat("NONE\n")
  
} else {
  
  for (f in count_candidates) {
    cat(f, "\n")
  }
}


# ==============================================================================
# 8. Inspect candidate count tables strictly
# ==============================================================================

candidate_audit <- data.frame()

if (length(count_candidates) > 0) {
  
  for (f in count_candidates) {
    
    dat <- tryCatch(
      suppressWarnings(
        fread(
          f,
          data.table = FALSE,
          check.names = FALSE,
          nrows = 5,
          showProgress = FALSE
        )
      ),
      error = function(e) NULL
    )
    
    if (is.null(dat)) {
      next
    }
    
    cn <- tolower(
      trimws(
        colnames(dat)
      )
    )
    
    n_keep <- sum(
      KB_KEEP %in% cn
    )
    
    has_kb2 <- KB_EXCLUDED %in% cn
    
    tmp <- data.frame(
      File = f,
      N_retained_KB_columns = n_keep,
      Has_all_7_KB = n_keep == 7,
      Has_kb2 = has_kb2,
      N_columns = ncol(dat),
      First_column = colnames(dat)[1],
      stringsAsFactors = FALSE
    )
    
    candidate_audit <- rbind(
      candidate_audit,
      tmp
    )
  }
}

if (nrow(candidate_audit) > 0) {
  
  candidate_audit <- candidate_audit[
    order(
      -candidate_audit$Has_all_7_KB,
      candidate_audit$Has_kb2,
      candidate_audit$File
    ),
    ,
    drop = FALSE
  ]
  
  print(
    candidate_audit,
    row.names = FALSE
  )
}


# ==============================================================================
# 9. Identify strongest 7KB count-table candidates
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("STRICT 7KB COUNT TABLE CANDIDATES\n")
cat("============================================================\n")

strict_candidates <- candidate_audit[
  candidate_audit$Has_all_7_KB &
    !candidate_audit$Has_kb2,
  ,
  drop = FALSE
]

if (nrow(strict_candidates) == 0) {
  
  cat("NONE\n")
  
} else {
  
  print(
    strict_candidates,
    row.names = FALSE
  )
}


# ==============================================================================
# 10. Save outputs
# ==============================================================================

fwrite(
  kb,
  file.path(
    OUT_DIR,
    "negative_controls_7KB_read_depth_after_decontamination.tsv"
  ),
  sep = "\t"
)

fwrite(
  kb_summary,
  file.path(
    OUT_DIR,
    "negative_controls_7KB_read_depth_summary.tsv"
  ),
  sep = "\t"
)

if (nrow(candidate_audit) > 0) {
  
  fwrite(
    candidate_audit,
    file.path(
      OUT_DIR,
      "count_table_candidate_audit.tsv"
    ),
    sep = "\t"
  )
}


# ==============================================================================
# 11. Finish
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("OUTPUT DIRECTORY\n")
cat("============================================================\n")

cat(
  OUT_DIR,
  "\n"
)

cat("\nDone.\n")