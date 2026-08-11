
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_06_prepare_PICRUSt2_input_CAC_vs_nonCAC.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Build the PICRUSt2 input feature table, representative-sequence FASTA and
## metadata for the paired CAC23 vs nonCAC23 cohort.
## PICRUSt2 itself is run at the command line; the command is printed at the end.
############################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

PROJECT_ROOT <- PROJECT_ROOT
MAIN_OUT <- OUTPUT_ROOT  ## 7KB rerun output root

KEY_DIR <- file.path(
  OUTPUT_ROOT,
  "00_clean_data",
  "CA23_nonCA23_paired"
)

RAW_FASTA <- file.path(PROJECT_ROOT, "raw/sequences.fasta")

OUT_DIR <- file.path(
  MAIN_OUT,
  "PICRUSt2_CA23_vs_nonCA23_7KB_input"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_FEATURE_TABLE <- file.path(OUT_DIR, "feature_table_CA23_vs_nonCA23_7KB.tsv")
OUT_BIOM_INPUT <- file.path(OUT_DIR, "feature_table_CA23_vs_nonCA23_7KB_biom_input.tsv")
OUT_FASTA <- file.path(OUT_DIR, "sequences_CA23_vs_nonCA23_7KB.fasta")
OUT_METADATA <- file.path(OUT_DIR, "metadata_CA23_vs_nonCA23_7KB.tsv")
OUT_AUDIT <- file.path(OUT_DIR, "ASV_sequence_matching_audit_CA23_vs_nonCA23_7KB.tsv")
OUT_MANIFEST <- file.path(OUT_DIR, "PICRUSt2_input_manifest_CA23_vs_nonCA23_7KB.tsv")
OUT_RUN_SH <- file.path(OUT_DIR, "run_PICRUSt2_CA23_vs_nonCA23_7KB.sh")

stopifnot(dir.exists(KEY_DIR))
stopifnot(file.exists(RAW_FASTA))

count_file <- list.files(
  KEY_DIR,
  pattern = "asv_count.*\\.tsv$",
  full.names = TRUE,
  ignore.case = TRUE
)
count_file <- count_file[!str_detect(basename(count_file), "relative|relab|ra")]
if (length(count_file) == 0) stop("Cannot find ASV count table in: ", KEY_DIR)
count_file <- count_file[1]

metadata_file <- list.files(
  KEY_DIR,
  pattern = "metadata.*\\.tsv$",
  full.names = TRUE,
  ignore.case = TRUE
)
if (length(metadata_file) == 0) stop("Cannot find metadata table in: ", KEY_DIR)
metadata_file <- metadata_file[1]

read_fasta <- function(path) {
  lines <- readLines(path)
  header_idx <- which(str_starts(lines, ">"))
  ids <- str_remove(lines[header_idx], "^>") %>%
    str_split_fixed("\\s+", 2) %>%
    .[, 1]
  
  seqs <- map_chr(seq_along(header_idx), function(i) {
    start <- header_idx[i] + 1
    end <- if (i < length(header_idx)) header_idx[i + 1] - 1 else length(lines)
    paste(lines[start:end], collapse = "")
  })
  
  tibble(ASV_ID = ids, Sequence = seqs)
}

count_raw <- read_tsv(count_file, show_col_types = FALSE, progress = FALSE)
meta_raw <- read_tsv(metadata_file, show_col_types = FALSE, progress = FALSE)

asv_col <- names(count_raw)[1]
count_tbl <- count_raw %>%
  rename(ASV_ID = all_of(asv_col)) %>%
  mutate(ASV_ID = as.character(ASV_ID))

sample_id_col <- intersect(c("Sample_ID", "SampleID", "sample_id", "sample"), names(meta_raw))[1]
if (is.na(sample_id_col)) stop("Cannot find sample ID column in metadata.")

group_candidates <- names(meta_raw)[map_lgl(meta_raw, function(x) {
  vals <- unique(as.character(x))
  all(c("CA", "nonCA") %in% vals) || all(c("CA", "NonCA") %in% vals)
})]

if (length(group_candidates) == 0) {
  stop("Cannot find CA/nonCA group column in metadata.")
}

group_col <- group_candidates[1]

meta <- meta_raw %>%
  transmute(
    Sample_ID = as.character(.data[[sample_id_col]]),
    CA_status = as.character(.data[[group_col]])
  ) %>%
  mutate(
    CA_status = case_when(
      CA_status == "CA" ~ "CA",
      str_to_lower(CA_status) == "nonca" ~ "nonCA",
      TRUE ~ CA_status
    )
  ) %>%
  filter(CA_status %in% c("nonCA", "CA")) %>%
  mutate(CA_status = factor(CA_status, levels = c("nonCA", "CA"))) %>%
  arrange(CA_status, Sample_ID)

sample_cols <- intersect(meta$Sample_ID, names(count_tbl))
missing_samples <- setdiff(meta$Sample_ID, sample_cols)
if (length(missing_samples) > 0) {
  stop("Metadata samples missing from count table: ", paste(missing_samples, collapse = ", "))
}

meta <- meta %>%
  filter(Sample_ID %in% sample_cols) %>%
  arrange(CA_status, Sample_ID)

ordered_samples <- meta$Sample_ID

count_clean <- count_tbl %>%
  select(ASV_ID, all_of(ordered_samples))

count_clean[ordered_samples] <- lapply(count_clean[ordered_samples], function(x) {
  replace_na(as.numeric(x), 0)
})

count_clean <- count_clean %>%
  filter(rowSums(across(all_of(ordered_samples))) > 0)

fasta_tbl <- read_fasta(RAW_FASTA)

matched_asvs <- intersect(count_clean$ASV_ID, fasta_tbl$ASV_ID)
missing_asvs <- setdiff(count_clean$ASV_ID, fasta_tbl$ASV_ID)

count_picrust <- count_clean %>%
  filter(ASV_ID %in% matched_asvs) %>%
  arrange(match(ASV_ID, matched_asvs))

fasta_picrust <- fasta_tbl %>%
  filter(ASV_ID %in% count_picrust$ASV_ID) %>%
  arrange(match(ASV_ID, count_picrust$ASV_ID))

feature_table <- count_picrust
biom_input <- feature_table %>%
  rename(`#OTU ID` = ASV_ID)

write_tsv(feature_table, OUT_FEATURE_TABLE)
write_tsv(biom_input, OUT_BIOM_INPUT)
write_tsv(meta, OUT_METADATA)

fasta_lines <- as.vector(rbind(
  paste0(">", fasta_picrust$ASV_ID),
  fasta_picrust$Sequence
))
writeLines(fasta_lines, OUT_FASTA)

audit <- tibble(
  Metric = c(
    "CA_n",
    "nonCA_n",
    "Clean_ASVs_in_count_table",
    "Raw_FASTA_sequences",
    "Matched_ASVs_used_for_PICRUSt2",
    "Missing_clean_ASVs_from_FASTA"
  ),
  Value = c(
    sum(meta$CA_status == "CA"),
    sum(meta$CA_status == "nonCA"),
    nrow(count_clean),
    nrow(fasta_tbl),
    nrow(count_picrust),
    length(missing_asvs)
  )
)
write_tsv(audit, OUT_AUDIT)

run_lines <- c(
  "#!/bin/bash",
  paste0("cd '", OUT_DIR, "'"),
  "",
  "conda run -n picrust2 biom convert \\",
  "  -i feature_table_CA23_vs_nonCA23_7KB_biom_input.tsv \\",
  "  -o feature_table_CA23_vs_nonCA23_7KB.biom \\",
  "  --to-hdf5 \\",
  "  --table-type='OTU table'",
  "",
  "conda run -n picrust2 picrust2_pipeline.py \\",
  "  -s sequences_CA23_vs_nonCA23_7KB.fasta \\",
  "  -i feature_table_CA23_vs_nonCA23_7KB.biom \\",
  "  -o picrust2_out_CA23_vs_nonCA23_7KB \\",
  "  -p 8",
  "",
  "conda run -n picrust2 add_descriptions.py \\",
  "  -i picrust2_out_CA23_vs_nonCA23_7KB/pathways_out/path_abun_unstrat.tsv.gz \\",
  "  -o picrust2_out_CA23_vs_nonCA23_7KB/pathways_out/path_abun_unstrat_descrip.tsv.gz \\",
  "  -m METACYC"
)
writeLines(run_lines, OUT_RUN_SH)

manifest <- tibble(
  item = c(
    "strategy",
    "count_file",
    "metadata_file",
    "raw_fasta",
    "feature_table",
    "biom_input",
    "fasta",
    "metadata",
    "audit",
    "run_script"
  ),
  value = c(
    "CA23 vs explicit nonCA23; 7KB decontaminated ASV table; raw FASTA filtered to retained clean ASVs only",
    count_file,
    metadata_file,
    RAW_FASTA,
    OUT_FEATURE_TABLE,
    OUT_BIOM_INPUT,
    OUT_FASTA,
    OUT_METADATA,
    OUT_AUDIT,
    OUT_RUN_SH
  )
)
write_tsv(manifest, OUT_MANIFEST)

cat("\n============================================================\n")
cat("PICRUSt2 input preparation for CA23 vs nonCA23 finished\n")
cat("============================================================\n\n")

cat("Group counts:\n")
print(meta %>% count(CA_status), n = Inf)

cat("\nASV matching summary:\n")
print(audit, n = Inf)

cat("\nFiles written:\n")
cat("Feature table: ", OUT_FEATURE_TABLE, "\n")
cat("BIOM input:    ", OUT_BIOM_INPUT, "\n")
cat("FASTA:         ", OUT_FASTA, "\n")
cat("Metadata:      ", OUT_METADATA, "\n")
cat("Audit:         ", OUT_AUDIT, "\n")
cat("Manifest:      ", OUT_MANIFEST, "\n")
cat("Run script:    ", OUT_RUN_SH, "\n\n")

cat("Run in Terminal:\n\n")
cat("bash ", OUT_RUN_SH, "\n", sep = "")
cat("\nDone.\n")