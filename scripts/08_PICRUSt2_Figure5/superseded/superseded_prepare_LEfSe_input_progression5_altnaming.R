


############################################################
## 13_prepare_LEfSe_input_MetaCyc_pathways_progression5_previous127_v0.5_step2.R
##
## Prepare LEfSe input for pathway-level MetaCyc analysis
##
## This script ONLY:
## - reads PICRUSt2 MetaCyc pathway abundance with descriptions
## - matches current previous127 5-group metadata
## - normalizes pathway abundance to CPM
## - filters low-prevalence pathways
## - writes LEfSe input
##
## It does NOT run LEfSe.
## Run LEfSe manually in Terminal using the commands printed at the end.
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
## Paths
############################################################

PROJECT_DIR <- PROJECT_ROOT
RUN_ID <- "decontamination_8controls"

PICRUST2_DIR <- file.path(
  PROJECT_DIR,
  "output",
  RUN_ID,
  "PICRUSt2_previous127_v0.5_step2_input",
  "picrust2_out_previous127_v0.5_step2"
)

PATHWAY_FILE <- file.path(
  PICRUST2_DIR,
  "pathways_out",
  "path_abun_unstrat_descrip.tsv.gz"
)

META_FILE <- file.path(
  PROJECT_DIR,
  "output",
  RUN_ID,
  "clean_tables_previous127_v0.5_step2_retention_flagged",
  "metadata_clean_v0.5_step2_previous127.tsv"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "output",
  RUN_ID,
  "LEfSe_MetaCyc_pathways_progression5_previous127_v0.5_step2"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

############################################################
## Parameters
############################################################

GROUP_ORDER <- c("Polyp", "UC_remission", "UC_active", "Dysplasia", "CA")
MIN_PREVALENCE_N <- 3
LDA_CUTOFF <- 2.0

############################################################
## Output files
############################################################

OUT_PATHWAY_CPM <- file.path(
  OUT_DIR,
  "MetaCyc_pathway_CPM_previous127_v0.5_step2.tsv"
)

OUT_PATHWAY_DESCRIPTION <- file.path(
  OUT_DIR,
  "MetaCyc_pathway_description_previous127_v0.5_step2.tsv"
)

OUT_FILTER_AUDIT <- file.path(
  OUT_DIR,
  "MetaCyc_pathway_filter_audit_previous127_v0.5_step2.tsv"
)

OUT_LEFSE_INPUT <- file.path(
  OUT_DIR,
  "LEfSe_input_MetaCyc_pathways_progression5_previous127_v0.5_step2.tsv"
)

OUT_TERMINAL_COMMANDS <- file.path(
  OUT_DIR,
  "run_LEfSe_MetaCyc_pathways_progression5_previous127_v0.5_step2.sh"
)

OUT_SUMMARY <- file.path(
  OUT_DIR,
  "summary_prepare_LEfSe_input_MetaCyc_pathways_progression5_previous127_v0.5_step2.tsv"
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

clean_pathway_id <- function(x) {
  str_trim(as.character(x))
}

clean_pathway_name <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\\|", "_")
  x <- str_replace_all(x, "\t", " ")
  x
}

make_feature_name <- function(pathway_id, description) {
  description <- if_else(
    is.na(description) | description == "",
    pathway_id,
    description
  )
  
  out <- paste0(pathway_id, "__", description)
  out <- str_replace_all(out, "[^A-Za-z0-9_\\.\\-]+", "_")
  out <- str_replace_all(out, "_+", "_")
  out <- str_replace_all(out, "^_|_$", "")
  out
}

############################################################
## Read data
############################################################

stop_if_missing(c(PATHWAY_FILE, META_FILE))

pathway_raw <- read_tsv(PATHWAY_FILE, show_col_types = FALSE)
meta_df <- read_tsv(META_FILE, show_col_types = FALSE)

meta_use <- meta_df %>%
  filter(Progression5 %in% GROUP_ORDER) %>%
  mutate(Progression5 = factor(Progression5, levels = GROUP_ORDER)) %>%
  arrange(Progression5, SampleID)

sample_keep <- meta_use$SampleID

############################################################
## Detect columns
############################################################

pathway_id_col <- colnames(pathway_raw)[1]

description_candidates <- c(
  "description",
  "Description",
  "pathway_description",
  "Pathway",
  "pathway"
)

description_col <- description_candidates[
  description_candidates %in% colnames(pathway_raw)
][1]

if (is.na(description_col)) {
  description_col <- tail(colnames(pathway_raw), 1)
  message("No standard description column detected; using last column as description: ", description_col)
}

sample_cols <- intersect(sample_keep, colnames(pathway_raw))
missing_samples <- setdiff(sample_keep, sample_cols)

if (length(missing_samples) > 0) {
  stop(
    "These previous127 samples are missing from pathway table:\n",
    paste(missing_samples, collapse = ", ")
  )
}

############################################################
## Prepare pathway matrix
############################################################

pathway_df <- pathway_raw %>%
  rename(Pathway_ID = all_of(pathway_id_col)) %>%
  mutate(
    Pathway_ID = clean_pathway_id(Pathway_ID),
    Description = clean_pathway_name(.data[[description_col]]),
    Feature = make_feature_name(Pathway_ID, Description)
  ) %>%
  select(Pathway_ID, Description, Feature, all_of(sample_keep))

pathway_mat <- pathway_df %>%
  select(Feature, all_of(sample_keep)) %>%
  column_to_rownames("Feature") %>%
  as.matrix()

storage.mode(pathway_mat) <- "numeric"
pathway_mat[is.na(pathway_mat)] <- 0

############################################################
## Normalize to CPM
############################################################

sample_sums <- colSums(pathway_mat)

if (any(sample_sums <= 0)) {
  stop(
    "Some samples have zero total pathway abundance:\n",
    paste(names(sample_sums)[sample_sums <= 0], collapse = ", ")
  )
}

pathway_cpm <- sweep(pathway_mat, 2, sample_sums, "/") * 1e6
pathway_cpm[is.na(pathway_cpm)] <- 0

write_tsv(
  pathway_cpm %>%
    as.data.frame(check.names = FALSE) %>%
    rownames_to_column("Feature"),
  OUT_PATHWAY_CPM
)

description_tbl <- pathway_df %>%
  select(Pathway_ID, Description, Feature)

write_tsv(description_tbl, OUT_PATHWAY_DESCRIPTION)

############################################################
## Filter low-prevalence pathways
############################################################

filter_audit <- tibble(
  Feature = rownames(pathway_cpm),
  prevalence_n = rowSums(pathway_cpm > 0),
  mean_CPM = rowMeans(pathway_cpm),
  max_CPM = apply(pathway_cpm, 1, max),
  keep_for_LEfSe = prevalence_n >= MIN_PREVALENCE_N & max_CPM > 0
) %>%
  left_join(description_tbl, by = "Feature") %>%
  select(Pathway_ID, Description, Feature, everything())

write_tsv(filter_audit, OUT_FILTER_AUDIT)

features_keep <- filter_audit %>%
  filter(keep_for_LEfSe) %>%
  pull(Feature)

pathway_cpm_keep <- pathway_cpm[features_keep, sample_keep, drop = FALSE]

############################################################
## Build LEfSe input
############################################################

class_row <- c("Class", as.character(meta_use$Progression5))

lefse_body <- pathway_cpm_keep %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Feature") %>%
  mutate(across(everything(), as.character))

lefse_input <- bind_rows(
  as_tibble_row(setNames(class_row, colnames(lefse_body))),
  lefse_body
)

write_tsv(lefse_input, OUT_LEFSE_INPUT, col_names = FALSE)

############################################################
## Write Terminal commands
############################################################

terminal_commands <- c(
  "#!/bin/bash",
  "conda activate lefse",
  paste0("cd ", shQuote(OUT_DIR)),
  "",
  "lefse_format_input.py \\",
  "  LEfSe_input_MetaCyc_pathways_progression5_previous127_v0.5_step2.tsv \\",
  "  LEfSe_formatted_MetaCyc_pathways_progression5_previous127_v0.5_step2.in \\",
  "  -c 1 \\",
  "  -o 1000000",
  "",
  "lefse_run.py \\",
  "  LEfSe_formatted_MetaCyc_pathways_progression5_previous127_v0.5_step2.in \\",
  "  LEfSe_results_MetaCyc_pathways_progression5_previous127_v0.5_step2.res \\",
  paste0("  -l ", LDA_CUTOFF)
)

writeLines(terminal_commands, OUT_TERMINAL_COMMANDS)

############################################################
## Summary
############################################################

summary_tbl <- tibble(
  Item = c(
    "Analysis",
    "Input pathway file",
    "Input metadata",
    "Number of samples",
    "Group counts",
    "Total pathways in PICRUSt2 table",
    "Pathways retained for LEfSe",
    "Minimum prevalence N",
    "LDA cutoff for downstream LEfSe",
    "LEfSe input",
    "Terminal command script",
    "Output directory"
  ),
  Value = c(
    "Prepare PICRUSt2 MetaCyc pathway LEfSe input, previous127 v0.5 Step2, 5-group progression",
    PATHWAY_FILE,
    META_FILE,
    as.character(nrow(meta_use)),
    paste(
      meta_use %>%
        count(Progression5, name = "n") %>%
        mutate(x = paste(Progression5, n, sep = "=")) %>%
        pull(x),
      collapse = "; "
    ),
    as.character(nrow(pathway_df)),
    as.character(length(features_keep)),
    as.character(MIN_PREVALENCE_N),
    as.character(LDA_CUTOFF),
    OUT_LEFSE_INPUT,
    OUT_TERMINAL_COMMANDS,
    OUT_DIR
  )
)

write_tsv(summary_tbl, OUT_SUMMARY)

############################################################
## Print summary
############################################################

cat("\n============================================================\n")
cat("MetaCyc pathway LEfSe input preparation finished\n")
cat("============================================================\n\n")

cat("Group counts:\n")
print(meta_use %>% count(Progression5, name = "n"), n = Inf, width = Inf)

cat("\nPathway counts:\n")
cat("Total pathways in PICRUSt2 table: ", nrow(pathway_df), "\n")
cat("Retained for LEfSe:               ", length(features_keep), "\n")

cat("\nLEfSe input written to:\n")
cat(OUT_LEFSE_INPUT, "\n")

cat("\nTerminal commands written to:\n")
cat(OUT_TERMINAL_COMMANDS, "\n")

cat("\nRun in Terminal:\n\n")
cat(paste(terminal_commands, collapse = "\n"))
cat("\n\nDone.\n")