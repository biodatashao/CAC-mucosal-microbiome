#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_09_negative_control_predecontam_composition.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Taxonomic composition of the negative controls before decontamination.
############################################################



options(stringsAsFactors = FALSE)
options(width = 260)

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# 1. Paths
# ==============================================================================

OLD_MASTER <- paste0(
  file.path(PROJECT_ROOT, "output/"),
  "decontamination_8controls"
)

OLD_DECONTAM <- file.path(
  OLD_MASTER,
  "decontamination"
)

WITHOUT_KB2 <- paste0(
  file.path(PROJECT_ROOT, "output/"),
  "decontamination_7controls_final"
)

WITHOUT_KB2_DECONTAM <- file.path(
  WITHOUT_KB2,
  "decontamination"
)

STEP0_FILE <- file.path(
  OLD_DECONTAM,
  "count_step0_non_target_removed_with_blanks.tsv"
)

RETENTION_FILE <- file.path(
  WITHOUT_KB2_DECONTAM,
  "sample_retention_fraction.tsv"
)

GENUS_BLANK_FILE <- file.path(
  WITHOUT_KB2_DECONTAM,
  "genus_blank_ratio_table.tsv"
)

OUT_DIR <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "05_manuscript_audit/negative_controls_7KB/"
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

KEY_GENERA <- c(
  "UCG-005",
  "Desulfovibrio",
  "Lactococcus",
  "Mediterraneibacter"
)


# ==============================================================================
# 2. Basic file checks
# ==============================================================================

required_files <- c(
  STEP0_FILE,
  RETENTION_FILE,
  GENUS_BLANK_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "Required file(s) missing:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# ==============================================================================
# 3. Read Step0 pre-decontamination matrix
# ==============================================================================

step0 <- fread(
  STEP0_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (!"ASV" %in% colnames(step0)) {
  stop("Step0 count table does not contain column 'ASV'.")
}

missing_kb_cols <- setdiff(
  KB_KEEP,
  colnames(step0)
)

if (length(missing_kb_cols) > 0) {
  
  stop(
    paste0(
      "Step0 table missing retained control columns: ",
      paste(
        missing_kb_cols,
        collapse = ", "
      )
    )
  )
}

if (!KB_EXCLUDED %in% colnames(step0)) {
  warning(
    "kb2 column was not found in historical Step0 matrix. Continuing."
  )
}

cat("\n")
cat("============================================================\n")
cat("7KB NEGATIVE CONTROL PRE-DECONTAMINATION AUDIT\n")
cat("============================================================\n")

cat(
  "Step0 source:\n",
  STEP0_FILE,
  "\n\n",
  sep = ""
)

cat(
  "Retained negative controls: ",
  paste(
    KB_KEEP,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Excluded negative control: ",
  KB_EXCLUDED,
  "\n",
  sep = ""
)

cat(
  "Step0 ASVs = ",
  nrow(step0),
  "\n",
  sep = ""
)


# ==============================================================================
# 4. Explicitly construct 7KB-only Step0 matrix
# ==============================================================================

step0_7kb <- step0[
  ,
  c(
    "ASV",
    KB_KEEP
  ),
  drop = FALSE
]

count_mat <- as.matrix(
  step0_7kb[
    ,
    KB_KEEP,
    drop = FALSE
  ]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- as.character(
  step0_7kb$ASV
)

if (any(!is.finite(count_mat))) {
  stop("Non-finite counts detected in Step0 control matrix.")
}

if (any(count_mat < 0)) {
  stop("Negative counts detected in Step0 control matrix.")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop("Step0 matrix contains non-integer counts.")
}

count_mat <- round(
  count_mat
)

step0_depth <- colSums(
  count_mat
)

detected_asv <- colSums(
  count_mat > 0
)


# ==============================================================================
# 5. Read current without-kb2 raw/post-step2 read summary
# ==============================================================================

ret <- fread(
  RETENTION_FILE,
  data.table = FALSE,
  check.names = FALSE
)

ret$SampleID <- tolower(
  trimws(
    as.character(
      ret$SampleID
    )
  )
)

if (KB_EXCLUDED %in% ret$SampleID) {
  stop("kb2 unexpectedly present in current without-kb2 retention table.")
}

ret_kb <- ret[
  ret$SampleID %in% KB_KEEP,
  ,
  drop = FALSE
]

if (nrow(ret_kb) != 7) {
  
  stop(
    paste0(
      "Expected 7 controls in retention table, found ",
      nrow(ret_kb),
      "."
    )
  )
}

ret_order <- match(
  KB_KEEP,
  ret_kb$SampleID
)

if (any(is.na(ret_order))) {
  stop("Could not match all seven KB controls in retention table.")
}

ret_kb <- ret_kb[
  ret_order,
  ,
  drop = FALSE
]


# ==============================================================================
# 6. Locate taxonomy candidates
#
# We need taxonomy for the Step0 ASV set, not merely the final clean ASV set.
# Candidate selection is therefore based on actual ASV overlap with Step0.
# ==============================================================================

all_master_files <- list.files(
  OLD_MASTER,
  recursive = TRUE,
  full.names = TRUE
)

taxonomy_candidates <- all_master_files[
  grepl(
    "taxonom",
    basename(all_master_files),
    ignore.case = TRUE
  ) &
    grepl(
      "\\.(tsv|csv)$",
      all_master_files,
      ignore.case = TRUE
    )
]

cat("\n")
cat("============================================================\n")
cat("TAXONOMY CANDIDATES\n")
cat("============================================================\n")

if (length(taxonomy_candidates) == 0) {
  stop("No taxonomy candidate files found under original master branch.")
}

taxonomy_audit <- data.frame(
  File = character(0),
  N_rows = integer(0),
  ASV_column = character(0),
  Taxonomy_column = character(0),
  Genus_column = character(0),
  N_step0_ASV_overlap = integer(0),
  Step0_ASV_coverage = numeric(0),
  stringsAsFactors = FALSE
)

taxonomy_data_list <- list()

step0_asvs <- as.character(
  step0_7kb$ASV
)

for (i in seq_along(taxonomy_candidates)) {
  
  file_i <- taxonomy_candidates[i]
  
  dat_i <- tryCatch(
    suppressWarnings(
      fread(
        file_i,
        data.table = FALSE,
        check.names = FALSE,
        showProgress = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat_i) || nrow(dat_i) == 0) {
    next
  }
  
  cn <- colnames(dat_i)
  
  asv_hits <- cn[
    tolower(cn) %in% c(
      "asv",
      "asv_id",
      "featureid",
      "feature_id"
    )
  ]
  
  if (length(asv_hits) == 0) {
    next
  }
  
  asv_col <- asv_hits[1]
  
  taxonomy_hits <- cn[
    grepl(
      "^taxonomy$|taxon",
      cn,
      ignore.case = TRUE
    )
  ]
  
  genus_hits <- cn[
    grepl(
      "^genus$|genus_name|genus_key",
      cn,
      ignore.case = TRUE
    )
  ]
  
  taxonomy_col <- if (
    length(taxonomy_hits) > 0
  ) {
    taxonomy_hits[1]
  } else {
    NA_character_
  }
  
  genus_col <- if (
    length(genus_hits) > 0
  ) {
    genus_hits[1]
  } else {
    NA_character_
  }
  
  overlap_n <- sum(
    step0_asvs %in% as.character(
      dat_i[, asv_col]
    )
  )
  
  coverage <- overlap_n /
    length(step0_asvs)
  
  taxonomy_audit <- rbind(
    taxonomy_audit,
    data.frame(
      File = file_i,
      N_rows = nrow(dat_i),
      ASV_column = asv_col,
      Taxonomy_column = taxonomy_col,
      Genus_column = genus_col,
      N_step0_ASV_overlap = overlap_n,
      Step0_ASV_coverage = coverage,
      stringsAsFactors = FALSE
    )
  )
  
  taxonomy_data_list[length(taxonomy_data_list) + 1] <- list(
    dat_i
  )
}

if (nrow(taxonomy_audit) == 0) {
  stop("No usable taxonomy table containing an ASV identifier column was found.")
}

taxonomy_audit <- taxonomy_audit[
  order(
    -taxonomy_audit$Step0_ASV_coverage,
    -taxonomy_audit$N_step0_ASV_overlap
  ),
  ,
  drop = FALSE
]

print(
  taxonomy_audit,
  row.names = FALSE
)

best_tax_file <- taxonomy_audit$File[1]

best_tax <- fread(
  best_tax_file,
  data.table = FALSE,
  check.names = FALSE
)

best_asv_col <- taxonomy_audit$ASV_column[1]
best_taxonomy_col <- taxonomy_audit$Taxonomy_column[1]
best_genus_col <- taxonomy_audit$Genus_column[1]
best_coverage <- taxonomy_audit$Step0_ASV_coverage[1]

cat("\n")
cat("============================================================\n")
cat("SELECTED TAXONOMY SOURCE\n")
cat("============================================================\n")

cat(
  "File:\n",
  best_tax_file,
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Step0 ASV coverage = %.2f%%\n",
    100 * best_coverage
  )
)

if (best_coverage < 0.95) {
  
  warning(
    paste0(
      "Selected taxonomy source covers only ",
      sprintf("%.2f%%", 100 * best_coverage),
      " of Step0 ASVs. Review before manuscript use."
    )
  )
}


# ==============================================================================
# 7. Genus parsing helpers
# ==============================================================================

clean_genus_label <- function(x) {
  
  x <- as.character(x)
  
  x <- trimws(x)
  
  x <- gsub(
    "^g__",
    "",
    x
  )
  
  x <- gsub(
    "^g_",
    "",
    x
  )
  
  x[
    is.na(x) |
      x == "" |
      tolower(x) %in% c(
        "na",
        "nan",
        "none",
        "unassigned"
      )
  ] <- "Unclassified"
  
  x
}


extract_genus_from_taxonomy <- function(x) {
  
  x <- as.character(x)
  
  out <- rep(
    "Unclassified",
    length(x)
  )
  
  for (i in seq_along(x)) {
    
    value_i <- x[i]
    
    if (is.na(value_i) || value_i == "") {
      next
    }
    
    pieces <- unlist(
      strsplit(
        value_i,
        ";",
        fixed = TRUE
      )
    )
    
    pieces <- trimws(
      pieces
    )
    
    genus_piece <- pieces[
      grepl(
        "^g__",
        pieces
      )
    ]
    
    if (length(genus_piece) > 0) {
      
      genus_value <- sub(
        "^g__",
        "",
        genus_piece[1]
      )
      
      if (
        !is.na(genus_value) &&
        genus_value != ""
      ) {
        out[i] <- genus_value
      }
    }
  }
  
  clean_genus_label(
    out
  )
}


normalize_genus_key <- function(x) {
  
  x <- tolower(
    trimws(
      as.character(x)
    )
  )
  
  x <- gsub(
    "\\[|\\]",
    "",
    x
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
  
  x <- gsub(
    "^g--",
    "",
    x
  )
  
  x <- gsub(
    "^g__",
    "",
    x
  )
  
  x
}


# ==============================================================================
# 8. Build Step0 ASV-to-genus map
# ==============================================================================

tax_map <- data.frame(
  ASV = as.character(
    best_tax[, best_asv_col]
  ),
  stringsAsFactors = FALSE
)

if (
  !is.na(best_genus_col) &&
  best_genus_col %in% colnames(best_tax)
) {
  
  tax_map$Genus <- clean_genus_label(
    best_tax[, best_genus_col]
  )
  
} else if (
  !is.na(best_taxonomy_col) &&
  best_taxonomy_col %in% colnames(best_tax)
) {
  
  tax_map$Genus <- extract_genus_from_taxonomy(
    best_tax[, best_taxonomy_col]
  )
  
} else {
  
  stop(
    paste0(
      "Selected taxonomy source contains neither a usable Genus column ",
      "nor a usable Taxonomy column."
    )
  )
}

tax_map <- tax_map[
  !duplicated(tax_map$ASV),
  ,
  drop = FALSE
]

tax_idx <- match(
  step0_asvs,
  tax_map$ASV
)

step0_genus <- tax_map$Genus[
  tax_idx
]

step0_genus[
  is.na(step0_genus)
] <- "Unclassified"

genus_key <- normalize_genus_key(
  step0_genus
)


# ==============================================================================
# 9. Genus-level composition per control
# ==============================================================================

genus_levels <- unique(
  step0_genus
)

genus_rows <- list()

for (sample_i in KB_KEEP) {
  
  counts_i <- count_mat[
    ,
    sample_i
  ]
  
  agg_i <- aggregate(
    counts_i,
    by = list(
      Genus = step0_genus
    ),
    FUN = sum
  )
  
  colnames(agg_i)[2] <- "Reads"
  
  agg_i <- agg_i[
    agg_i$Reads > 0,
    ,
    drop = FALSE
  ]
  
  total_i <- sum(
    agg_i$Reads
  )
  
  agg_i$SampleID <- sample_i
  agg_i$Total_step0_reads <- total_i
  
  if (total_i > 0) {
    
    agg_i$Relative_abundance <- agg_i$Reads /
      total_i
    
  } else {
    
    agg_i$Relative_abundance <- NA_real_
  }
  
  agg_i$Percent <- 100 *
    agg_i$Relative_abundance
  
  agg_i <- agg_i[
    order(
      -agg_i$Reads,
      agg_i$Genus
    ),
    ,
    drop = FALSE
  ]
  
  agg_i$Rank <- seq_len(
    nrow(agg_i)
  )
  
  genus_rows[length(genus_rows) + 1] <- list(
    agg_i
  )
}

genus_composition <- do.call(
  rbind,
  genus_rows
)

genus_composition <- genus_composition[
  ,
  c(
    "SampleID",
    "Rank",
    "Genus",
    "Reads",
    "Total_step0_reads",
    "Relative_abundance",
    "Percent"
  ),
  drop = FALSE
]


# ==============================================================================
# 10. Top 10 genera per control
# ==============================================================================

top10 <- genus_composition[
  genus_composition$Rank <= 10,
  ,
  drop = FALSE
]


# ==============================================================================
# 11. Key taxa audit
# ==============================================================================

key_genus_keys <- normalize_genus_key(
  KEY_GENERA
)

key_rows <- list()

for (sample_i in KB_KEEP) {
  
  sample_comp <- genus_composition[
    genus_composition$SampleID == sample_i,
    ,
    drop = FALSE
  ]
  
  sample_keys <- normalize_genus_key(
    sample_comp$Genus
  )
  
  for (j in seq_along(KEY_GENERA)) {
    
    target <- KEY_GENERA[j]
    target_key <- key_genus_keys[j]
    
    hit <- which(
      sample_keys == target_key
    )
    
    if (length(hit) == 0) {
      
      reads_j <- 0
      ra_j <- 0
      pct_j <- 0
      
    } else {
      
      reads_j <- sum(
        sample_comp$Reads[
          hit
        ]
      )
      
      ra_j <- sum(
        sample_comp$Relative_abundance[
          hit
        ]
      )
      
      pct_j <- 100 * ra_j
    }
    
    tmp <- data.frame(
      SampleID = sample_i,
      Genus = target,
      Reads = reads_j,
      Relative_abundance = ra_j,
      Percent = pct_j,
      stringsAsFactors = FALSE
    )
    
    key_rows[length(key_rows) + 1] <- list(
      tmp
    )
  }
}

key_taxa <- do.call(
  rbind,
  key_rows
)


# ==============================================================================
# 12. Per-control summary
# ==============================================================================

control_summary <- data.frame(
  SampleID = KB_KEEP,
  Raw_reads = ret_kb$Raw_reads,
  Step0_reads = as.numeric(
    step0_depth[
      KB_KEEP
    ]
  ),
  Detected_ASVs_step0 = as.numeric(
    detected_asv[
      KB_KEEP
    ]
  ),
  Reads_after_step2 = ret_kb$Reads_after_step2,
  Retention_fraction_step2 = ret_kb$Retention_fraction,
  stringsAsFactors = FALSE
)

for (target in KEY_GENERA) {
  
  target_key <- normalize_genus_key(
    target
  )
  
  key_sub <- key_taxa[
    normalize_genus_key(
      key_taxa$Genus
    ) == target_key,
    ,
    drop = FALSE
  ]
  
  idx <- match(
    control_summary$SampleID,
    key_sub$SampleID
  )
  
  safe_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    target
  )
  
  safe_name <- gsub(
    "^_|_$",
    "",
    safe_name
  )
  
  control_summary[
    ,
    paste0(
      safe_name,
      "_reads"
    )
  ] <- key_sub$Reads[
    idx
  ]
  
  control_summary[
    ,
    paste0(
      safe_name,
      "_percent"
    )
  ] <- key_sub$Percent[
    idx
  ]
}


# ==============================================================================
# 13. Seven-control summary statistics for key taxa
# ==============================================================================

key_summary_rows <- list()

for (target in KEY_GENERA) {
  
  target_sub <- key_taxa[
    key_taxa$Genus == target,
    ,
    drop = FALSE
  ]
  
  tmp <- data.frame(
    Genus = target,
    Controls_positive = sum(
      target_sub$Reads > 0
    ),
    Controls_total = 7,
    Total_reads_across_controls = sum(
      target_sub$Reads
    ),
    Median_reads = median(
      target_sub$Reads
    ),
    Max_reads = max(
      target_sub$Reads
    ),
    Median_percent = median(
      target_sub$Percent
    ),
    Max_percent = max(
      target_sub$Percent
    ),
    stringsAsFactors = FALSE
  )
  
  key_summary_rows[length(key_summary_rows) + 1] <- list(
    tmp
  )
}

key_summary <- do.call(
  rbind,
  key_summary_rows
)


# ==============================================================================
# 14. Read current without-kb2 genus blank-vs-biological evidence
# ==============================================================================

blank_genus <- fread(
  GENUS_BLANK_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (!"Genus" %in% colnames(blank_genus)) {
  stop("genus_blank_ratio_table.tsv lacks Genus column.")
}

blank_genus$Genus_normalized_key <- normalize_genus_key(
  blank_genus$Genus
)

evidence_rows <- list()

for (target in KEY_GENERA) {
  
  target_key <- normalize_genus_key(
    target
  )
  
  hit <- which(
    blank_genus$Genus_normalized_key ==
      target_key
  )
  
  if (length(hit) == 0) {
    
    tmp <- data.frame(
      Requested_genus = target,
      Match_found = FALSE,
      stringsAsFactors = FALSE
    )
    
  } else {
    
    tmp <- blank_genus[
      hit,
      ,
      drop = FALSE
    ]
    
    tmp$Requested_genus <- target
    tmp$Match_found <- TRUE
    
    tmp <- tmp[
      ,
      c(
        "Requested_genus",
        "Match_found",
        setdiff(
          colnames(tmp),
          c(
            "Requested_genus",
            "Match_found",
            "Genus_normalized_key"
          )
        )
      ),
      drop = FALSE
    ]
  }
  
  evidence_rows[length(evidence_rows) + 1] <- list(
    tmp
  )
}

blank_vs_bio_key_evidence <- rbindlist(
  evidence_rows,
  fill = TRUE
)

blank_vs_bio_key_evidence <- as.data.frame(
  blank_vs_bio_key_evidence
)


# ==============================================================================
# 15. Special UCG-005 console interpretation inputs
# ==============================================================================

ucg_control <- key_taxa[
  key_taxa$Genus == "UCG-005",
  ,
  drop = FALSE
]

ucg_positive_n <- sum(
  ucg_control$Reads > 0
)

ucg_total_reads <- sum(
  ucg_control$Reads
)

ucg_max_percent <- max(
  ucg_control$Percent
)


# ==============================================================================
# 16. Console output
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("PER-CONTROL PRE-DECONTAMINATION SUMMARY\n")
cat("============================================================\n")

print(
  control_summary,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("KEY TAXA IN 7 NEGATIVE CONTROLS - STEP0\n")
cat("============================================================\n")

print(
  key_taxa,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("KEY TAXA SUMMARY ACROSS 7 NEGATIVE CONTROLS\n")
cat("============================================================\n")

print(
  key_summary,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("UCG-005 NEGATIVE CONTROL SUMMARY\n")
cat("============================================================\n")

cat(
  "Controls positive for UCG-005 = ",
  ucg_positive_n,
  "/7\n",
  sep = ""
)

cat(
  "Total UCG-005 reads across 7 controls = ",
  ucg_total_reads,
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Maximum UCG-005 relative abundance in any control = %.8f%%\n",
    ucg_max_percent
  )
)


cat("\n")
cat("============================================================\n")
cat("WITHOUT-kb2 BLANK VS BIOLOGICAL GENUS EVIDENCE\n")
cat("============================================================\n")

print(
  blank_vs_bio_key_evidence,
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("TOP 10 GENERA PER NEGATIVE CONTROL\n")
cat("============================================================\n")

print(
  top10[
    ,
    c(
      "SampleID",
      "Rank",
      "Genus",
      "Reads",
      "Percent"
    ),
    drop = FALSE
  ],
  row.names = FALSE
)


# ==============================================================================
# 17. Save outputs
# ==============================================================================

fwrite(
  taxonomy_audit,
  file.path(
    OUT_DIR,
    "26b4_taxonomy_candidate_overlap_audit_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  control_summary,
  file.path(
    OUT_DIR,
    "26b4_negative_controls_predecontam_summary_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  genus_composition,
  file.path(
    OUT_DIR,
    "26b4_negative_controls_predecontam_genus_composition_full_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  top10,
  file.path(
    OUT_DIR,
    "26b4_negative_controls_predecontam_top10_genera_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  key_taxa,
  file.path(
    OUT_DIR,
    "26b4_negative_controls_key_taxa_by_control_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  key_summary,
  file.path(
    OUT_DIR,
    "26b4_negative_controls_key_taxa_summary_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  blank_vs_bio_key_evidence,
  file.path(
    OUT_DIR,
    "26b4_key_taxa_blank_vs_biological_evidence.tsv"
  ),
  sep = "\t"
)


# ==============================================================================
# 18. Save method/audit note
# ==============================================================================

method_note <- c(
  "Negative-control audit definition:",
  "",
  paste0(
    "Retained controls: ",
    paste(
      KB_KEEP,
      collapse = ", "
    )
  ),
  "kb2 was explicitly excluded.",
  "",
  paste0(
    "Pre-decontamination community composition was calculated from: ",
    STEP0_FILE
  ),
  "",
  paste0(
    "This Step0 matrix represents counts after non-target removal but before ",
    "decontam-based ASV removal and genus blacklist removal."
  ),
  "",
  paste0(
    "Raw read counts and post-Step2 retained read counts were obtained from: ",
    RETENTION_FILE
  ),
  "",
  paste0(
    "Blank-versus-biological genus evidence after excluding kb2 was obtained from: ",
    GENUS_BLANK_FILE
  ),
  "",
  paste0(
    "Selected taxonomy source: ",
    best_tax_file
  ),
  "",
  paste0(
    "Selected taxonomy source Step0 ASV coverage: ",
    sprintf(
      "%.4f%%",
      100 * best_coverage
    )
  )
)

writeLines(
  method_note,
  con = file.path(
    OUT_DIR,
    "26b4_negative_control_analysis_method_note_7KB.txt"
  )
)


# ==============================================================================
# 19. Session info
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_DIR,
    "26b4_sessionInfo_7KB.txt"
  )
)


# ==============================================================================
# 20. Finish
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