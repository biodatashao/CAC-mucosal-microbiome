
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_07_prepare_LEfSe_input_pathways_CAC_vs_nonCAC.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Prepare the LEfSe input for pathway-level comparison of CAC vs nonCAC.
## LEfSe itself is run at the command line; the command is printed at the end.
############################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

PROJECT_ROOT <- PROJECT_ROOT
MAIN_OUT <- OUTPUT_ROOT  ## 7KB rerun output root

PICRUST2_INPUT_DIR <- file.path(
  MAIN_OUT,
  "PICRUSt2_CA23_vs_nonCA23_7KB_input"
)

PATHWAY_FILE <- file.path(
  PICRUST2_INPUT_DIR,
  "picrust2_out_CA23_vs_nonCA23_7KB/pathways_out/path_abun_unstrat_descrip.tsv.gz"
)

METADATA_FILE <- file.path(
  PICRUST2_INPUT_DIR,
  "metadata_CA23_vs_nonCA23_7KB.tsv"
)

OUT_DIR <- file.path(
  MAIN_OUT,
  "LEfSe_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_INPUT <- file.path(
  OUT_DIR,
  "LEfSe_input_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.tsv"
)
OUT_FEATURE_MAP <- file.path(
  OUT_DIR,
  "MetaCyc_pathway_feature_name_map_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_CPM <- file.path(
  OUT_DIR,
  "MetaCyc_pathway_CPM_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_META <- file.path(
  OUT_DIR,
  "metadata_LEfSe_pathway_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_MANIFEST <- file.path(
  OUT_DIR,
  "manifest_LEfSe_pathway_CA23_vs_nonCA23_oldstyle_7KB.tsv"
)
OUT_RUN_SH <- file.path(
  OUT_DIR,
  "run_LEfSe_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.sh"
)

stopifnot(file.exists(PATHWAY_FILE))
stopifnot(file.exists(METADATA_FILE))

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

pathway_raw <- read_tsv(PATHWAY_FILE, show_col_types = FALSE, progress = FALSE)
meta_raw <- read_tsv(METADATA_FILE, show_col_types = FALSE, progress = FALSE)

pathway_col <- names(pathway_raw)[1]
desc_col <- names(pathway_raw)[str_detect(str_to_lower(names(pathway_raw)), "description|desc")]
desc_col <- if (length(desc_col) > 0) desc_col[1] else names(pathway_raw)[2]

sample_id_col <- intersect(c("Sample_ID", "SampleID", "sample_id", "sample"), names(meta_raw))[1]
if (is.na(sample_id_col)) stop("Cannot find sample ID column in metadata.")

group_col <- intersect(c("CA_status", "Group", "Group2", "CA_group"), names(meta_raw))[1]
if (is.na(group_col)) {
  group_candidates <- names(meta_raw)[map_lgl(meta_raw, function(x) {
    vals <- unique(as.character(x))
    all(c("CA", "nonCA") %in% vals) || all(c("CA", "NonCA") %in% vals)
  })]
  if (length(group_candidates) == 0) stop("Cannot find CA/nonCA group column.")
  group_col <- group_candidates[1]
}

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
    ),
    CA_status = factor(CA_status, levels = c("nonCA", "CA"))
  ) %>%
  filter(!is.na(CA_status)) %>%
  arrange(CA_status, Sample_ID)

sample_cols <- intersect(meta$Sample_ID, names(pathway_raw))
missing_samples <- setdiff(meta$Sample_ID, sample_cols)

if (length(missing_samples) > 0) {
  stop("Metadata samples missing from pathway table: ", paste(missing_samples, collapse = ", "))
}

meta <- meta %>%
  filter(Sample_ID %in% sample_cols) %>%
  arrange(CA_status, Sample_ID)

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

abun[ordered_samples] <- lapply(abun[ordered_samples], function(x) {
  replace_na(as.numeric(x), 0)
})

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
  tibble(Feature = "class", !!!set_names(as.list(as.character(meta$CA_status)), ordered_samples)),
  tibble(Feature = "subclass", !!!set_names(as.list(as.character(meta$CA_status)), ordered_samples)),
  tibble(Feature = "subject", !!!set_names(as.list(as.character(meta$Sample_ID)), ordered_samples)),
  abun_cpm_chr
)

write_tsv(lefse_rows, OUT_INPUT, col_names = FALSE)
write_csv(feature_map, OUT_FEATURE_MAP)
write_csv(abun_cpm, OUT_CPM)
write_csv(meta, OUT_META)

run_lines <- c(
  "#!/bin/bash",
  paste0("cd '", OUT_DIR, "'"),
  "",
  "conda run -n lefse lefse_format_input.py \\",
  "  LEfSe_input_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.tsv \\",
  "  LEfSe_formatted_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.in \\",
  "  -c 1 \\",
  "  -s 2 \\",
  "  -u 3 \\",
  "  -o 1000000",
  "",
  "conda run -n lefse lefse_run.py \\",
  "  LEfSe_formatted_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.in \\",
  "  LEfSe_results_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.res \\",
  "  -l 2.0"
)

writeLines(run_lines, OUT_RUN_SH)

manifest <- tibble(
  item = c(
    "strategy",
    "comparison",
    "pathway_file",
    "metadata_file",
    "lefse_input",
    "feature_map",
    "cpm_table",
    "metadata_out",
    "run_script"
  ),
  value = c(
    "CA23 vs explicit nonCA23; old-style LEfSe input: class/subclass/subject + CPM + pathway-description feature names",
    "nonCA vs CA",
    PATHWAY_FILE,
    METADATA_FILE,
    OUT_INPUT,
    OUT_FEATURE_MAP,
    OUT_CPM,
    OUT_META,
    OUT_RUN_SH
  )
)

write_tsv(manifest, OUT_MANIFEST)

cat("\n============================================================\n")
cat("CA vs nonCA MetaCyc pathway LEfSe input preparation finished\n")
cat("============================================================\n\n")

cat("Group counts:\n")
print(meta %>% count(CA_status), n = Inf)

cat("\nPathway counts:\n")
cat("Total pathways:", nrow(pathway_tbl), "\n")

cat("\nFiles written:\n")
cat("LEfSe input: ", OUT_INPUT, "\n")
cat("Feature map: ", OUT_FEATURE_MAP, "\n")
cat("CPM table:   ", OUT_CPM, "\n")
cat("Metadata:    ", OUT_META, "\n")
cat("Manifest:    ", OUT_MANIFEST, "\n")
cat("Run script:  ", OUT_RUN_SH, "\n\n")

cat("Run in Terminal:\n\n")
cat("bash ", OUT_RUN_SH, "\n", sep = "")
cat("\nDone.\n")