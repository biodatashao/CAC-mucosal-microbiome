


############################################################
## 08_02_prepare_PICRUSt2_input_progression127.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
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
  library(stringr)
})

############################################################
## Paths
############################################################

PROJECT_DIR <- PROJECT_ROOT

INPUT_DIR <- file.path(
  OUTPUT_ROOT,
  "00_clean_data",
  "progression127"
)

RAW_DIR <- file.path(PROJECT_DIR, "raw")

OUT_DIR <- file.path(
  OUTPUT_ROOT,
  "08_PICRUSt2_pathways",
  "PICRUSt2_progression127_7KB_input"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COUNT_FILE <- file.path(INPUT_DIR, "asv_count_7KB_progression127.tsv")
META_FILE  <- file.path(INPUT_DIR, "metadata_7KB_progression127.tsv")
RAW_FASTA_FILE <- file.path(RAW_DIR, "sequences.fasta")

OUT_FEATURE_TABLE <- file.path(
  OUT_DIR,
  "feature_table_progression127_7KB.tsv"
)

OUT_METADATA <- file.path(
  OUT_DIR,
  "metadata_progression127_7KB.tsv"
)

OUT_FASTA <- file.path(
  OUT_DIR,
  "sequences_progression127_7KB.fasta"
)

OUT_ASV_AUDIT <- file.path(
  OUT_DIR,
  "ASV_sequence_matching_audit_progression127_7KB.tsv"
)

OUT_MANIFEST <- file.path(
  OUT_DIR,
  "PICRUSt2_input_manifest_progression127_7KB.tsv"
)

############################################################
## Helper functions
############################################################

stop_if_missing <- function(files) {
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("Missing required file(s):\n", paste(missing_files, collapse = "\n"))
  }
}

read_fasta_simple <- function(fasta_file) {
  lines <- readLines(fasta_file, warn = FALSE)
  
  header_idx <- which(str_starts(lines, ">"))
  
  if (length(header_idx) == 0) {
    stop("No FASTA headers found in:\n", fasta_file)
  }
  
  end_idx <- c(header_idx[-1] - 1, length(lines))
  
  out <- vector("list", length(header_idx))
  
  for (i in seq_along(header_idx)) {
    h <- lines[header_idx[i]]
    seq_lines <- lines[(header_idx[i] + 1):end_idx[i]]
    
    id <- str_remove(h, "^>")
    id <- str_split(id, "\\s+", simplify = TRUE)[1]
    
    out[[i]] <- tibble(
      ASV = id,
      Header = h,
      Sequence = paste(seq_lines, collapse = "")
    )
  }
  
  bind_rows(out)
}

write_fasta_simple <- function(fasta_df, output_file) {
  con <- file(output_file, open = "wt")
  on.exit(close(con), add = TRUE)
  
  for (i in seq_len(nrow(fasta_df))) {
    writeLines(paste0(">", fasta_df$ASV[i]), con)
    seq <- fasta_df$Sequence[i]
    
    chunks <- str_extract_all(seq, ".{1,80}")[[1]]
    writeLines(chunks, con)
  }
}

############################################################
## Read input files
############################################################

stop_if_missing(c(COUNT_FILE, META_FILE, RAW_FASTA_FILE))

count_df <- read_tsv(COUNT_FILE, show_col_types = FALSE)
meta_df  <- read_tsv(META_FILE, show_col_types = FALSE)

fasta_df <- read_fasta_simple(RAW_FASTA_FILE)

############################################################
## Check and prepare count table
############################################################

asv_col <- colnames(count_df)[1]

count_clean <- count_df %>%
  rename(ASV = all_of(asv_col))

clean_asvs <- count_clean$ASV

GROUP_ORDER <- c("Polyp", "UC_remission", "UC_active", "Dysplasia", "CA")

meta_use <- meta_df %>%
  filter(Progression5 %in% GROUP_ORDER) %>%
  mutate(
    Progression5 = factor(Progression5, levels = GROUP_ORDER)
  ) %>%
  arrange(Progression5, SampleID)

sample_keep <- meta_use$SampleID

missing_samples <- setdiff(sample_keep, colnames(count_clean))

if (length(missing_samples) > 0) {
  stop(
    "Metadata samples missing from count table:\n",
    paste(missing_samples, collapse = ", ")
  )
}

count_clean <- count_clean %>%
  select(ASV, all_of(sample_keep))

############################################################
## Match ASVs to FASTA sequences
############################################################

fasta_df <- fasta_df %>%
  distinct(ASV, .keep_all = TRUE)

asv_audit <- tibble(
  ASV = clean_asvs,
  Present_in_raw_fasta = clean_asvs %in% fasta_df$ASV
)

write_tsv(asv_audit, OUT_ASV_AUDIT)

missing_asv <- asv_audit %>%
  filter(!Present_in_raw_fasta) %>%
  pull(ASV)

if (length(missing_asv) > 0) {
  warning(
    "Some clean ASVs were not found in raw FASTA and will be excluded from PICRUSt2 input:\n",
    paste(head(missing_asv, 50), collapse = ", "),
    ifelse(length(missing_asv) > 50, "\n... truncated", "")
  )
}

matched_asvs <- asv_audit %>%
  filter(Present_in_raw_fasta) %>%
  pull(ASV)

count_picrust <- count_clean %>%
  filter(ASV %in% matched_asvs)

fasta_picrust <- fasta_df %>%
  filter(ASV %in% matched_asvs) %>%
  arrange(match(ASV, count_picrust$ASV))

count_picrust <- count_picrust %>%
  arrange(match(ASV, fasta_picrust$ASV))

stopifnot(identical(count_picrust$ASV, fasta_picrust$ASV))

############################################################
## Write outputs
############################################################

write_tsv(count_picrust, OUT_FEATURE_TABLE)
write_tsv(meta_use, OUT_METADATA)
write_fasta_simple(fasta_picrust, OUT_FASTA)

manifest_tbl <- tibble(
  Item = c(
    "Analysis",
    "Project directory",
    "Run ID",
    "Input clean ASV count",
    "Input clean metadata",
    "Input raw FASTA",
    "Output PICRUSt2 feature table",
    "Output PICRUSt2 FASTA",
    "Output PICRUSt2 metadata",
    "Output ASV audit",
    "Number of progression127 samples",
    "Group counts",
    "Clean ASVs in count table",
    "Raw FASTA sequences",
    "Matched ASVs used for PICRUSt2",
    "Missing clean ASVs from FASTA"
  ),
  Value = c(
    "PICRUSt2 input preparation for progression127 7KB",
    PROJECT_DIR,
    OUTPUT_ROOT,
    COUNT_FILE,
    META_FILE,
    RAW_FASTA_FILE,
    OUT_FEATURE_TABLE,
    OUT_FASTA,
    OUT_METADATA,
    OUT_ASV_AUDIT,
    as.character(nrow(meta_use)),
    paste(
      meta_use %>%
        count(Progression5, name = "n") %>%
        mutate(x = paste(Progression5, n, sep = "=")) %>%
        pull(x),
      collapse = "; "
    ),
    as.character(length(clean_asvs)),
    as.character(nrow(fasta_df)),
    as.character(length(matched_asvs)),
    as.character(length(missing_asv))
  )
)

write_tsv(manifest_tbl, OUT_MANIFEST)

############################################################
## Print summary
############################################################

cat("\n============================================================\n")
cat("PICRUSt2 input preparation finished\n")
cat("============================================================\n\n")

cat("Group counts:\n")
print(meta_use %>% count(Progression5, name = "n"), n = Inf, width = Inf)

cat("\nASV matching summary:\n")
cat("Clean ASVs in count table:        ", length(clean_asvs), "\n")
cat("Raw FASTA sequences:              ", nrow(fasta_df), "\n")
cat("Matched ASVs used for PICRUSt2:   ", length(matched_asvs), "\n")
cat("Missing clean ASVs from FASTA:    ", length(missing_asv), "\n")

cat("\nFiles written:\n")
cat("Feature table: ", OUT_FEATURE_TABLE, "\n")
cat("FASTA:         ", OUT_FASTA, "\n")
cat("Metadata:      ", OUT_METADATA, "\n")
cat("ASV audit:     ", OUT_ASV_AUDIT, "\n")
cat("Manifest:      ", OUT_MANIFEST, "\n")

cat("\nDone.\n")