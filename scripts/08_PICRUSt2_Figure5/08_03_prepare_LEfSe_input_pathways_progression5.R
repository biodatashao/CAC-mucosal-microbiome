
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_03_prepare_LEfSe_input_pathways_progression5.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Convert PICRUSt2 MetaCyc pathway abundances to CPM, filter low-prevalence
## pathways and write the LEfSe input for the five progression groups.
## LEfSe itself is run at the command line; the command is printed at the end.
############################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

PROJECT_ROOT <- PROJECT_ROOT
MAIN_OUT <- OUTPUT_ROOT  ## 7KB rerun output root

PICRUST2_PATHWAY_FILE <- file.path(
  MAIN_OUT,
  "PICRUSt2_progression127_7KB_input/picrust2_out_progression127_7KB/pathways_out/path_abun_unstrat_descrip.tsv.gz"
)

META_FILE <- file.path(
  MAIN_OUT,
  "00_clean_data/progression127/metadata_7KB_progression127.tsv"
)

OUT_DIR <- file.path(
  MAIN_OUT,
  "LEfSe_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_LEVELS <- c("Polyp", "UC_remission", "UC_active", "Dysplasia", "CA")

OUT_INPUT <- file.path(
  OUT_DIR,
  "LEfSe_input_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.tsv"
)
OUT_FEATURE_MAP <- file.path(
  OUT_DIR,
  "MetaCyc_pathway_feature_name_map_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_CPM <- file.path(
  OUT_DIR,
  "MetaCyc_pathway_CPM_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_META <- file.path(
  OUT_DIR,
  "metadata_LEfSe_pathway_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_MANIFEST <- file.path(
  OUT_DIR,
  "manifest_LEfSe_pathway_oldstyle_currentdata_progression127_7KB.tsv"
)
OUT_RUN_SH <- file.path(
  OUT_DIR,
  "run_LEfSe_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.sh"
)

stopifnot(file.exists(PICRUST2_PATHWAY_FILE))
stopifnot(file.exists(META_FILE))

clean_feature_name <- function(pathway, description) {
  paste(pathway, description, sep = "__") %>%
    str_replace_all("&beta;", "beta") %>%
    str_replace_all("&alpha;", "alpha") %>%
    str_replace_all("&gamma;", "gamma") %>%
    str_replace_all("&delta;", "delta") %>%
    str_replace_all("[^A-Za-z0-9_]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "")
}

pathway_raw <- read_tsv(PICRUST2_PATHWAY_FILE, show_col_types = FALSE, progress = FALSE)
meta_raw <- read_tsv(META_FILE, show_col_types = FALSE, progress = FALSE)

pathway_col <- names(pathway_raw)[1]
desc_col <- names(pathway_raw)[str_detect(str_to_lower(names(pathway_raw)), "description|desc")]
desc_col <- if (length(desc_col) > 0) desc_col[1] else names(pathway_raw)[2]

sample_id_col <- intersect(c("Sample_ID", "SampleID", "sample_id", "sample"), names(meta_raw))[1]
if (is.na(sample_id_col)) stop("Cannot find sample ID column in metadata.")

group_col <- intersect(c("Progression5", "Group5", "group5"), names(meta_raw))[1]
if (is.na(group_col)) stop("Cannot find Progression5/Group5 column in metadata.")

meta <- meta_raw %>%
  transmute(
    Sample_ID = as.character(.data[[sample_id_col]]),
    Group5 = as.character(.data[[group_col]])
  ) %>%
  filter(Group5 %in% GROUP_LEVELS) %>%
  mutate(Group5 = factor(Group5, levels = GROUP_LEVELS)) %>%
  arrange(Group5, Sample_ID)

sample_cols <- intersect(meta$Sample_ID, names(pathway_raw))
missing_samples <- setdiff(meta$Sample_ID, sample_cols)

if (length(missing_samples) > 0) {
  stop("Samples in metadata but missing from pathway table: ", paste(missing_samples, collapse = ", "))
}

meta <- meta %>%
  filter(Sample_ID %in% sample_cols) %>%
  arrange(Group5, Sample_ID)

ordered_samples <- meta$Sample_ID

pathway_tbl <- pathway_raw %>%
  transmute(
    Pathway_ID = as.character(.data[[pathway_col]]),
    Description = as.character(.data[[desc_col]]),
    across(all_of(ordered_samples), ~ suppressWarnings(as.numeric(.x)))
  ) %>%
  mutate(
    Description = if_else(is.na(Description) | Description == "", Pathway_ID, Description),
    Feature = clean_feature_name(Pathway_ID, Description)
  )

if (anyDuplicated(pathway_tbl$Feature) > 0) {
  pathway_tbl <- pathway_tbl %>%
    group_by(Feature) %>%
    mutate(Feature = if_else(n() > 1, paste0(Feature, "__dup", row_number()), Feature)) %>%
    ungroup()
}

abun <- pathway_tbl %>%
  select(Feature, all_of(ordered_samples))

abun[ordered_samples] <- lapply(abun[ordered_samples], function(x) replace_na(as.numeric(x), 0))

sample_sums <- colSums(abun[ordered_samples], na.rm = TRUE)
if (any(sample_sums <= 0)) {
  stop("Samples with zero total pathway abundance: ",
       paste(names(sample_sums)[sample_sums <= 0], collapse = ", "))
}

abun_cpm <- abun
abun_cpm[ordered_samples] <- sweep(abun[ordered_samples], 2, sample_sums, "/") * 1e6

feature_map <- pathway_tbl %>%
  select(Feature, Pathway_ID, Description)

abun_cpm_chr <- abun_cpm %>%
  mutate(across(all_of(ordered_samples), as.character))

lefse_rows <- bind_rows(
  tibble(Feature = "class", !!!set_names(as.list(as.character(meta$Group5)), ordered_samples)),
  tibble(Feature = "subclass", !!!set_names(as.list(as.character(meta$Group5)), ordered_samples)),
  tibble(Feature = "subject", !!!set_names(as.list(as.character(meta$Sample_ID)), ordered_samples)),
  abun_cpm_chr
)

write_tsv(lefse_rows, OUT_INPUT, col_names = FALSE)
write_csv(feature_map, OUT_FEATURE_MAP)
write_csv(abun_cpm, OUT_CPM)
write_csv(meta, OUT_META)

run_lines <- c(
  "#!/bin/bash",
  "conda activate lefse",
  paste0("cd '", OUT_DIR, "'"),
  "",
  "lefse_format_input.py \\",
  "  LEfSe_input_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.tsv \\",
  "  LEfSe_formatted_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.in \\",
  "  -c 1 \\",
  "  -s 2 \\",
  "  -u 3 \\",
  "  -o 1000000",
  "",
  "lefse_run.py \\",
  "  LEfSe_formatted_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.in \\",
  "  LEfSe_results_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.res \\",
  "  -l 2.0"
)

writeLines(run_lines, OUT_RUN_SH)

manifest <- tibble(
  item = c(
    "strategy",
    "data_version",
    "picrust2_pathway_file",
    "metadata_file",
    "lefse_input",
    "feature_map",
    "cpm_table",
    "metadata_out",
    "run_shell_script"
  ),
  value = c(
    "Current 7KB progression127 data + old-style LEfSe input: class/subclass/subject + CPM + pathway-description feature names",
    "7KB_progression127",
    PICRUST2_PATHWAY_FILE,
    META_FILE,
    OUT_INPUT,
    OUT_FEATURE_MAP,
    OUT_CPM,
    OUT_META,
    OUT_RUN_SH
  )
)
write_tsv(manifest, OUT_MANIFEST)

cat("\n============================================================\n")
cat("FINAL old-style MetaCyc pathway LEfSe input preparation finished\n")
cat("============================================================\n\n")

cat("Group counts:\n")
print(meta %>% count(Group5), n = Inf)

cat("\nPathway counts:\n")
cat("Total pathways:", nrow(pathway_tbl), "\n")

cat("\nFiles written:\n")
cat("LEfSe input:   ", OUT_INPUT, "\n")
cat("Feature map:   ", OUT_FEATURE_MAP, "\n")
cat("CPM table:     ", OUT_CPM, "\n")
cat("Metadata:      ", OUT_META, "\n")
cat("Manifest:      ", OUT_MANIFEST, "\n")
cat("Run script:    ", OUT_RUN_SH, "\n\n")

cat("Run in Terminal:\n\n")
cat(paste(run_lines, collapse = "\n"))
cat("\n\nDone.\n")