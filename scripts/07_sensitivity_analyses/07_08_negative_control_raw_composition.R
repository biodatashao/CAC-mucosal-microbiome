#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_08_negative_control_raw_composition.R
##
## Module 07 - Sensitivity analyses and negative controls
##
## Taxonomic composition of the negative controls, raw reads.
############################################################



options(stringsAsFactors = FALSE)
options(width = 260)

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

META_MAIN_FILE <- file.path(
  RAW_DIR,
  "metadata_main.tsv"
)

META_FEATURE_FILE <- file.path(
  RAW_DIR,
  "metadata_from_featuretable.tsv"
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

KEY_GENERA <- c(
  "UCG-005",
  "Desulfovibrio",
  "Lactococcus",
  "Mediterraneibacter"
)

# ==============================================================================
# 2. Check files
# ==============================================================================

required_files <- c(
  COUNT_FILE,
  TAX_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required raw file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

# ==============================================================================
# 3. Read raw ASV count table
# ==============================================================================

count_df <- fread(
  COUNT_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(count_df) < 2) {
  stop("Raw ASV count table has fewer than 2 columns.")
}

# Detect ASV column
asv_col_candidates <- colnames(count_df)[
  tolower(colnames(count_df)) %in% c(
    "asv",
    "asv_id",
    "featureid",
    "feature_id",
    "#otu id",
    "otu"
  )
]

if (length(asv_col_candidates) == 0) {
  asv_col <- colnames(count_df)[1]
} else {
  asv_col <- asv_col_candidates[1]
}

sample_cols <- setdiff(
  colnames(count_df),
  asv_col
)

missing_kb <- setdiff(
  KB_KEEP,
  sample_cols
)

if (length(missing_kb) > 0) {
  stop(
    paste0(
      "Raw count table is missing retained KB columns: ",
      paste(
        missing_kb,
        collapse = ", "
      )
    )
  )
}

if (!KB_EXCLUDED %in% sample_cols) {
  warning("kb2 not found in raw count table.")
}

cat("\n")
cat("============================================================\n")
cat("7KB RAW NEGATIVE CONTROL COMPOSITION AUDIT\n")
cat("============================================================\n")

cat(
  "Count source:\n",
  COUNT_FILE,
  "\n",
  sep = ""
)

cat(
  "Taxonomy source:\n",
  TAX_FILE,
  "\n",
  sep = ""
)

cat(
  "ASV ID column = ",
  asv_col,
  "\n",
  sep = ""
)

cat(
  "Raw ASVs = ",
  nrow(count_df),
  "\n",
  sep = ""
)

cat(
  "Retained controls = ",
  paste(
    KB_KEEP,
    collapse = ", "
  ),
  "\n",
  sep = ""
)

cat(
  "Excluded control = ",
  KB_EXCLUDED,
  "\n",
  sep = ""
)

# ==============================================================================
# 4. Construct 7KB raw count matrix
# ==============================================================================

raw_7kb <- count_df[
  ,
  c(
    asv_col,
    KB_KEEP
  ),
  drop = FALSE
]

colnames(raw_7kb)[1] <- "ASV"

count_mat <- as.matrix(
  raw_7kb[
    ,
    KB_KEEP,
    drop = FALSE
  ]
)

storage.mode(count_mat) <- "numeric"

rownames(count_mat) <- as.character(
  raw_7kb$ASV
)

if (any(!is.finite(count_mat))) {
  stop("Non-finite values detected in raw count matrix.")
}

if (any(count_mat < 0)) {
  stop("Negative counts detected in raw count matrix.")
}

if (any(abs(count_mat - round(count_mat)) > 1e-8)) {
  stop("Raw count matrix contains non-integer values.")
}

count_mat <- round(
  count_mat
)

raw_depth <- colSums(
  count_mat
)

detected_asv <- colSums(
  count_mat > 0
)

# ==============================================================================
# 5. Read taxonomy
# ==============================================================================

tax_df <- fread(
  TAX_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(tax_df) < 2) {
  stop("Taxonomy table has fewer than 2 columns.")
}

tax_asv_candidates <- colnames(tax_df)[
  tolower(colnames(tax_df)) %in% c(
    "asv",
    "asv_id",
    "featureid",
    "feature_id",
    "#otu id",
    "otu"
  )
]

if (length(tax_asv_candidates) == 0) {
  tax_asv_col <- colnames(tax_df)[1]
} else {
  tax_asv_col <- tax_asv_candidates[1]
}

taxonomy_candidates <- colnames(tax_df)[
  grepl(
    "taxonomy|taxon",
    colnames(tax_df),
    ignore.case = TRUE
  )
]

genus_candidates <- colnames(tax_df)[
  grepl(
    "^genus$|genus_name|genus_key",
    colnames(tax_df),
    ignore.case = TRUE
  )
]

cat("\n")
cat("============================================================\n")
cat("RAW TAXONOMY STRUCTURE\n")
cat("============================================================\n")

cat(
  "Taxonomy ASV column = ",
  tax_asv_col,
  "\n",
  sep = ""
)

cat(
  "Taxonomy-like columns = ",
  ifelse(
    length(taxonomy_candidates) == 0,
    "NONE",
    paste(
      taxonomy_candidates,
      collapse = " | "
    )
  ),
  "\n",
  sep = ""
)

cat(
  "Genus-like columns = ",
  ifelse(
    length(genus_candidates) == 0,
    "NONE",
    paste(
      genus_candidates,
      collapse = " | "
    )
  ),
  "\n",
  sep = ""
)

# ==============================================================================
# 6. Genus helpers
# ==============================================================================

clean_genus <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x <- sub(
    "^g__",
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
    
    pieces <- trimws(
      unlist(
        strsplit(
          value_i,
          ";",
          fixed = TRUE
        )
      )
    )
    
    genus_hit <- pieces[
      grepl(
        "^g__",
        pieces
      )
    ]
    
    if (length(genus_hit) > 0) {
      
      genus_i <- sub(
        "^g__",
        "",
        genus_hit[1]
      )
      
      if (!is.na(genus_i) && genus_i != "") {
        out[i] <- genus_i
      }
    }
  }
  
  clean_genus(
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
  
  x <- sub(
    "^g__",
    "",
    x
  )
  
  x
}

# ==============================================================================
# 7. Build taxonomy map
# ==============================================================================

tax_map <- data.frame(
  ASV = as.character(
    tax_df[, tax_asv_col]
  ),
  stringsAsFactors = FALSE
)

if (length(genus_candidates) > 0) {
  
  genus_col <- genus_candidates[1]
  
  tax_map$Genus <- clean_genus(
    tax_df[, genus_col]
  )
  
} else if (length(taxonomy_candidates) > 0) {
  
  taxonomy_col <- taxonomy_candidates[1]
  
  tax_map$Genus <- extract_genus_from_taxonomy(
    tax_df[, taxonomy_col]
  )
  
} else {
  
  stop(
    "Taxonomy file contains neither a usable genus column nor taxonomy string column."
  )
}

tax_map <- tax_map[
  !duplicated(tax_map$ASV),
  ,
  drop = FALSE
]

tax_idx <- match(
  rownames(count_mat),
  tax_map$ASV
)

genus_vec <- tax_map$Genus[
  tax_idx
]

genus_vec[
  is.na(genus_vec)
] <- "Unclassified"

taxonomy_coverage <- mean(
  !is.na(tax_idx)
)

cat("\n")
cat("============================================================\n")
cat("RAW TAXONOMY COVERAGE\n")
cat("============================================================\n")

cat(
  sprintf(
    "ASV taxonomy coverage = %.4f%%\n",
    100 * taxonomy_coverage
  )
)

cat(
  "ASVs without taxonomy match = ",
  sum(is.na(tax_idx)),
  "\n",
  sep = ""
)

# ==============================================================================
# 8. Per-control raw summary
# ==============================================================================

control_summary <- data.frame(
  SampleID = KB_KEEP,
  Raw_ASV_reads = as.numeric(
    raw_depth[
      KB_KEEP
    ]
  ),
  Detected_ASVs = as.numeric(
    detected_asv[
      KB_KEEP
    ]
  ),
  stringsAsFactors = FALSE
)

# ==============================================================================
# 9. Build full genus composition
# ==============================================================================

genus_rows <- list()

for (sample_i in KB_KEEP) {
  
  counts_i <- count_mat[
    ,
    sample_i
  ]
  
  tmp <- data.frame(
    Genus = genus_vec,
    Reads = counts_i,
    stringsAsFactors = FALSE
  )
  
  agg <- aggregate(
    Reads ~ Genus,
    data = tmp,
    FUN = sum
  )
  
  agg <- agg[
    agg$Reads > 0,
    ,
    drop = FALSE
  ]
  
  total_i <- sum(
    agg$Reads
  )
  
  agg$SampleID <- sample_i
  agg$Total_reads <- total_i
  agg$Relative_abundance <- agg$Reads / total_i
  agg$Percent <- 100 *
    agg$Relative_abundance
  
  agg <- agg[
    order(
      -agg$Reads,
      agg$Genus
    ),
    ,
    drop = FALSE
  ]
  
  agg$Rank <- seq_len(
    nrow(agg)
  )
  
  genus_rows[length(genus_rows) + 1] <- list(
    agg
  )
}

genus_comp <- do.call(
  rbind,
  genus_rows
)

genus_comp <- genus_comp[
  ,
  c(
    "SampleID",
    "Rank",
    "Genus",
    "Reads",
    "Total_reads",
    "Relative_abundance",
    "Percent"
  ),
  drop = FALSE
]

# ==============================================================================
# 10. Top 10 genera per control
# ==============================================================================

top10 <- genus_comp[
  genus_comp$Rank <= 10,
  ,
  drop = FALSE
]

# ==============================================================================
# 11. Key taxa audit
# ==============================================================================

key_rows <- list()

for (sample_i in KB_KEEP) {
  
  comp_i <- genus_comp[
    genus_comp$SampleID == sample_i,
    ,
    drop = FALSE
  ]
  
  comp_keys <- normalize_genus_key(
    comp_i$Genus
  )
  
  for (target in KEY_GENERA) {
    
    target_key <- normalize_genus_key(
      target
    )
    
    hit <- which(
      comp_keys == target_key
    )
    
    if (length(hit) == 0) {
      
      reads_i <- 0
      ra_i <- 0
      pct_i <- 0
      
    } else {
      
      reads_i <- sum(
        comp_i$Reads[
          hit
        ]
      )
      
      ra_i <- sum(
        comp_i$Relative_abundance[
          hit
        ]
      )
      
      pct_i <- 100 * ra_i
    }
    
    key_rows[length(key_rows) + 1] <- list(
      data.frame(
        SampleID = sample_i,
        Genus = target,
        Reads = reads_i,
        Relative_abundance = ra_i,
        Percent = pct_i,
        stringsAsFactors = FALSE
      )
    )
  }
}

key_taxa <- do.call(
  rbind,
  key_rows
)

# ==============================================================================
# 12. Add key taxa to per-control summary
# ==============================================================================

for (target in KEY_GENERA) {
  
  target_sub <- key_taxa[
    key_taxa$Genus == target,
    ,
    drop = FALSE
  ]
  
  idx <- match(
    control_summary$SampleID,
    target_sub$SampleID
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
  ] <- target_sub$Reads[
    idx
  ]
  
  control_summary[
    ,
    paste0(
      safe_name,
      "_percent"
    )
  ] <- target_sub$Percent[
    idx
  ]
}

# ==============================================================================
# 13. Key taxa summary across 7 controls
# ==============================================================================

key_summary_rows <- list()

for (target in KEY_GENERA) {
  
  x <- key_taxa[
    key_taxa$Genus == target,
    ,
    drop = FALSE
  ]
  
  key_summary_rows[length(key_summary_rows) + 1] <- list(
    data.frame(
      Genus = target,
      Controls_positive = sum(
        x$Reads > 0
      ),
      Controls_total = 7,
      Total_reads_across_controls = sum(
        x$Reads
      ),
      Median_reads = median(
        x$Reads
      ),
      Max_reads = max(
        x$Reads
      ),
      Median_percent = median(
        x$Percent
      ),
      Max_percent = max(
        x$Percent
      ),
      stringsAsFactors = FALSE
    )
  )
}

key_summary <- do.call(
  rbind,
  key_summary_rows
)

# ==============================================================================
# 14. UCG-005 ASV-level audit
# ==============================================================================

ucg_key <- normalize_genus_key(
  "UCG-005"
)

ucg_asv_mask <- normalize_genus_key(
  genus_vec
) == ucg_key

ucg_asvs <- rownames(count_mat)[
  ucg_asv_mask
]

ucg_asv_count <- count_mat[
  ucg_asv_mask,
  ,
  drop = FALSE
]

if (length(ucg_asvs) == 0) {
  
  ucg_asv_long <- data.frame(
    ASV = character(0),
    SampleID = character(0),
    Reads = numeric(0),
    stringsAsFactors = FALSE
  )
  
} else {
  
  ucg_asv_long <- data.frame()
  
  for (sample_i in KB_KEEP) {
    
    tmp <- data.frame(
      ASV = ucg_asvs,
      SampleID = sample_i,
      Reads = as.numeric(
        ucg_asv_count[
          ,
          sample_i
        ]
      ),
      stringsAsFactors = FALSE
    )
    
    ucg_asv_long <- rbind(
      ucg_asv_long,
      tmp
    )
  }
  
  ucg_asv_long <- ucg_asv_long[
    ucg_asv_long$Reads > 0,
    ,
    drop = FALSE
  ]
}

# ==============================================================================
# 15. Console output
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("PER-CONTROL RAW SUMMARY\n")
cat("============================================================\n")

print(
  control_summary,
  row.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("KEY TAXA BY NEGATIVE CONTROL\n")
cat("============================================================\n")

print(
  key_taxa,
  row.names = FALSE
)

cat("\n")
cat("============================================================\n")
cat("KEY TAXA SUMMARY ACROSS 7 CONTROLS\n")
cat("============================================================\n")

print(
  key_summary,
  row.names = FALSE
)

ucg_summary <- key_summary[
  key_summary$Genus == "UCG-005",
  ,
  drop = FALSE
]

cat("\n")
cat("============================================================\n")
cat("UCG-005 RAW NEGATIVE CONTROL RESULT\n")
cat("============================================================\n")

cat(
  "UCG-005 ASVs in raw taxonomy = ",
  length(ucg_asvs),
  "\n",
  sep = ""
)

cat(
  "Controls positive = ",
  ucg_summary$Controls_positive,
  "/7\n",
  sep = ""
)

cat(
  "Total UCG-005 reads across controls = ",
  ucg_summary$Total_reads_across_controls,
  "\n",
  sep = ""
)

cat(
  sprintf(
    "Median UCG-005 abundance = %.8f%%\n",
    ucg_summary$Median_percent
  )
)

cat(
  sprintf(
    "Maximum UCG-005 abundance = %.8f%%\n",
    ucg_summary$Max_percent
  )
)

cat("\n")
cat("============================================================\n")
cat("TOP 10 GENERA PER CONTROL\n")
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
# 16. Save outputs
# ==============================================================================

fwrite(
  control_summary,
  file.path(
    OUT_DIR,
    "26b4_RAW_negative_controls_summary_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  genus_comp,
  file.path(
    OUT_DIR,
    "26b4_RAW_negative_controls_genus_composition_full_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  top10,
  file.path(
    OUT_DIR,
    "26b4_RAW_negative_controls_top10_genera_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  key_taxa,
  file.path(
    OUT_DIR,
    "26b4_RAW_negative_controls_key_taxa_by_control_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  key_summary,
  file.path(
    OUT_DIR,
    "26b4_RAW_negative_controls_key_taxa_summary_7KB.tsv"
  ),
  sep = "\t"
)

fwrite(
  ucg_asv_long,
  file.path(
    OUT_DIR,
    "26b4_RAW_UCG005_ASV_reads_in_negative_controls_7KB.tsv"
  ),
  sep = "\t"
)

# ==============================================================================
# 17. Metadata audit, if available
# ==============================================================================

metadata_audit <- data.frame()

meta_files <- c(
  META_MAIN_FILE,
  META_FEATURE_FILE
)

meta_files <- meta_files[
  file.exists(meta_files)
]

for (f in meta_files) {
  
  dat <- tryCatch(
    fread(
      f,
      data.table = FALSE,
      check.names = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(dat) || nrow(dat) == 0) {
    next
  }
  
  sample_candidates <- colnames(dat)[
    grepl(
      "sample",
      colnames(dat),
      ignore.case = TRUE
    )
  ]
  
  if (length(sample_candidates) == 0) {
    next
  }
  
  sample_col <- sample_candidates[1]
  
  sample_vals <- tolower(
    trimws(
      as.character(
        dat[, sample_col]
      )
    )
  )
  
  n_keep_found <- sum(
    KB_KEEP %in% sample_vals
  )
  
  has_kb2 <- KB_EXCLUDED %in% sample_vals
  
  metadata_audit <- rbind(
    metadata_audit,
    data.frame(
      File = f,
      Sample_column = sample_col,
      N_retained_KB_found = n_keep_found,
      kb2_present = has_kb2,
      stringsAsFactors = FALSE
    )
  )
}

if (nrow(metadata_audit) > 0) {
  
  fwrite(
    metadata_audit,
    file.path(
      OUT_DIR,
      "26b4_RAW_metadata_KB_audit_7KB.tsv"
    ),
    sep = "\t"
  )
}

# ==============================================================================
# 18. Method note
# ==============================================================================

method_note <- c(
  "Negative-control raw-data audit",
  "",
  paste0(
    "Retained negative controls: ",
    paste(
      KB_KEEP,
      collapse = ", "
    )
  ),
  paste0(
    "Excluded control: ",
    KB_EXCLUDED
  ),
  "",
  paste0(
    "ASV count source: ",
    COUNT_FILE
  ),
  paste0(
    "Taxonomy source: ",
    TAX_FILE
  ),
  "",
  paste0(
    "Raw negative-control composition was calculated directly from the raw ASV count table ",
    "before decontam-based contaminant removal and before genus blacklist removal."
  ),
  "",
  paste0(
    "Taxonomy coverage of raw ASVs: ",
    sprintf(
      "%.4f%%",
      100 * taxonomy_coverage
    )
  )
)

writeLines(
  method_note,
  con = file.path(
    OUT_DIR,
    "26b4_RAW_negative_control_method_note_7KB.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT_DIR,
    "26b4_RAW_sessionInfo_7KB.txt"
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