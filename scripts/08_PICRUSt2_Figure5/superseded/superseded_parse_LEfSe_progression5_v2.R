
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

PROJECT_ROOT <- PROJECT_ROOT
MAIN_OUT <- file.path(PROJECT_ROOT, "output/decontamination_8controls")

LEFSE_DIR <- file.path(
  MAIN_OUT,
  "LEfSe_MetaCyc_pathways_progression5_previous127_v0.5_step2"
)

PICRUST2_DIR <- file.path(
  MAIN_OUT,
  "PICRUSt2_previous127_v0.5_step2_input/picrust2_out_previous127_v0.5_step2/pathways_out"
)

RES_FILE <- file.path(
  LEFSE_DIR,
  "LEfSe_results_MetaCyc_pathways_progression5_previous127_v0.5_step2.res"
)

PATH_DESC_FILE <- file.path(
  PICRUST2_DIR,
  "path_abun_unstrat_descrip.tsv.gz"
)

GROUP_LEVELS <- c("Polyp", "UC_remission", "UC_active", "Dysplasia", "CA")
LDA_CUTOFF <- 2

OUT_CLEAN <- file.path(
  LEFSE_DIR,
  "FINAL_v2_LEfSe_clean_markers_MetaCyc_pathways_progression5_previous127_v0.5_step2.csv"
)
OUT_COUNTS <- file.path(
  LEFSE_DIR,
  "FINAL_v2_LEfSe_marker_counts_by_group_MetaCyc_pathways_progression5_previous127_v0.5_step2.csv"
)
OUT_TOP20 <- file.path(
  LEFSE_DIR,
  "FINAL_v2_LEfSe_top20_by_group_MetaCyc_pathways_progression5_previous127_v0.5_step2.csv"
)
OUT_TOP8 <- file.path(
  LEFSE_DIR,
  "FINAL_v2_LEfSe_top8_representative_candidates_MetaCyc_pathways_progression5_previous127_v0.5_step2.csv"
)
OUT_AUDIT <- file.path(
  LEFSE_DIR,
  "FINAL_v2_MetaCyc_description_matching_audit_previous127_v0.5_step2.csv"
)
OUT_SUMMARY <- file.path(
  LEFSE_DIR,
  "FINAL_v2_summary_parse_LEfSe_MetaCyc_pathways_progression5_previous127_v0.5_step2.tsv"
)

stopifnot(file.exists(RES_FILE))
stopifnot(file.exists(PATH_DESC_FILE))

make_key <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("&beta;", "beta") %>%
    str_replace_all("&alpha;", "alpha") %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    str_to_lower()
}

display_text <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("&beta;", "beta") %>%
    str_replace_all("&alpha;", "alpha") %>%
    str_replace_all("_+", " ") %>%
    str_squish()
}

classify_module <- function(pathway_id, description) {
  x <- str_to_lower(paste(pathway_id, description, sep = " "))
  
  case_when(
    str_detect(x, "ethanolamine|choline") ~ "Ethanolamine / choline metabolism",
    str_detect(x, "peptidoglycan|muramoyl|cell wall|lipid iva|kdo|phospholipid|diacylglycerol|fatty acid|hopanoid|membrane") ~ "Lipid / membrane / cell envelope",
    str_detect(x, "glycolysis|gluconeogenesis|pentose phosphate|pyruvate|fermentation|acetate|lactate|butanoate|butyrate|starch|sucrose|glycogen|carbon|methylcitrate|ketogenesis") ~ "Central carbon / fermentation",
    str_detect(x, "purine|pyrimidine|nucleotide|nucleobase|ribonucleotide|guanosine|inosine|trna") ~ "Nucleotide / translation",
    str_detect(x, "methionine|cysteine|tryptophan|tyrosine|isoleucine|lysine|threonine|aspartate|amino acid") ~ "Amino acid metabolism",
    str_detect(x, "sulfate|sulfur") ~ "Sulfur metabolism",
    str_detect(x, "nitrate|nitrogen|denitrification") ~ "Nitrogen metabolism",
    str_detect(x, "cobalamin|vitamin|cofactor|menaquinol|menaquinone|ubiquinol|heme") ~ "Cofactor / quinone metabolism",
    str_detect(x, "glycan|hexitol|galactarate") ~ "Carbohydrate / glycan metabolism",
    TRUE ~ "Other"
  )
}

lefse_raw <- read_tsv(
  RES_FILE,
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)

if (ncol(lefse_raw) < 5) {
  stop("LEfSe result file should contain at least 5 columns.")
}

lefse <- lefse_raw %>%
  select(1:5) %>%
  set_names(c("Feature", "log_max_mean", "Class", "LDA", "p_value")) %>%
  mutate(
    Feature = as.character(Feature),
    log_max_mean = suppressWarnings(as.numeric(log_max_mean)),
    Class = as.character(Class),
    LDA = suppressWarnings(as.numeric(LDA)),
    p_value = suppressWarnings(as.numeric(p_value)),
    abs_LDA = abs(LDA),
    Feature_key = make_key(Feature),
    row_id = row_number()
  ) %>%
  filter(
    !is.na(Feature),
    !is.na(Class),
    Class %in% GROUP_LEVELS,
    !is.na(LDA),
    abs_LDA >= LDA_CUTOFF
  )

path_desc_raw <- read_tsv(
  PATH_DESC_FILE,
  show_col_types = FALSE,
  progress = FALSE
)

id_col <- names(path_desc_raw)[1]
desc_col_candidates <- names(path_desc_raw)[str_detect(str_to_lower(names(path_desc_raw)), "description|desc")]
desc_col <- if (length(desc_col_candidates) > 0) desc_col_candidates[1] else names(path_desc_raw)[2]

path_desc <- path_desc_raw %>%
  transmute(
    Pathway_ID = as.character(.data[[id_col]]),
    Description = as.character(.data[[desc_col]])
  ) %>%
  filter(!is.na(Pathway_ID), Pathway_ID != "") %>%
  mutate(
    Description = if_else(is.na(Description) | Description == "", Pathway_ID, Description),
    Feature_original = paste0(
      Pathway_ID,
      "_",
      Description %>%
        str_replace_all("&beta;", "beta") %>%
        str_replace_all("&alpha;", "alpha") %>%
        str_replace_all("[^A-Za-z0-9]+", "_") %>%
        str_replace_all("_+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )

desc_keys <- path_desc %>%
  mutate(
    key_id_raw = make_key(Pathway_ID),
    key_feature_original = make_key(Feature_original),
    key_id_description = make_key(paste(Pathway_ID, Description, sep = "_")),
    key_id_description_no_pwy_suffix = make_key(
      paste(str_replace(Pathway_ID, "-PWY$", "_PWY"), Description, sep = "_")
    )
  ) %>%
  pivot_longer(
    cols = starts_with("key_"),
    names_to = "Key_type",
    values_to = "Match_key"
  ) %>%
  mutate(
    Match_priority = case_when(
      Key_type == "key_feature_original" ~ 1L,
      Key_type == "key_id_description" ~ 2L,
      Key_type == "key_id_description_no_pwy_suffix" ~ 3L,
      Key_type == "key_id_raw" ~ 4L,
      TRUE ~ 9L
    )
  ) %>%
  filter(!is.na(Match_key), Match_key != "") %>%
  arrange(Match_key, Match_priority) %>%
  distinct(Match_key, .keep_all = TRUE)

clean_markers <- lefse %>%
  left_join(desc_keys, by = c("Feature_key" = "Match_key")) %>%
  arrange(row_id, Match_priority) %>%
  group_by(row_id) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    Description_matched = !is.na(Pathway_ID),
    Pathway_ID = if_else(
      Description_matched,
      Pathway_ID,
      str_replace(Feature, "_.*$", "") %>% str_replace_all("_", "-")
    ),
    Description = if_else(
      Description_matched,
      Description,
      Feature %>%
        str_replace("^[A-Za-z0-9_]+_PWY_", "") %>%
        str_replace("^[A-Za-z0-9]+_[0-9]+_", "") %>%
        display_text()
    ),
    Description_display = display_text(Description),
    Pathway_module = classify_module(Pathway_ID, Description),
    Class = factor(Class, levels = GROUP_LEVELS)
  ) %>%
  select(
    Class,
    Pathway_ID,
    Description,
    Description_display,
    Pathway_module,
    Feature,
    Feature_original,
    Description_matched,
    log_max_mean,
    LDA,
    abs_LDA,
    p_value
  ) %>%
  arrange(Class, desc(abs_LDA), p_value)

marker_counts <- clean_markers %>%
  count(Class, name = "n_markers") %>%
  complete(Class = factor(GROUP_LEVELS, levels = GROUP_LEVELS), fill = list(n_markers = 0)) %>%
  arrange(Class)

top20_by_group <- clean_markers %>%
  group_by(Class) %>%
  arrange(desc(abs_LDA), p_value, .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()

top8_candidates <- clean_markers %>%
  group_by(Class) %>%
  arrange(desc(abs_LDA), p_value, .by_group = TRUE) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  mutate(
    Representative_selection_note =
      "Initial candidate: top 8 by absolute LDA within enriched group; manual redundancy/module review required"
  )

description_audit <- clean_markers %>%
  count(Description_matched, name = "n") %>%
  mutate(percent = round(100 * n / sum(n), 1))

description_summary <- tibble(
  Total_markers = nrow(clean_markers),
  Description_matched_n = sum(clean_markers$Description_matched, na.rm = TRUE),
  Description_unmatched_n = sum(!clean_markers$Description_matched, na.rm = TRUE),
  Description_matched_percent = round(
    100 * Description_matched_n / Total_markers,
    1
  )
)

summary_tbl <- tibble(
  metric = c(
    "Total LEfSe markers retained",
    "Description matched",
    "Description unmatched",
    "Description matched percent",
    "LDA cutoff"
  ),
  value = c(
    nrow(clean_markers),
    description_summary$Description_matched_n,
    description_summary$Description_unmatched_n,
    description_summary$Description_matched_percent,
    LDA_CUTOFF
  )
)

write_csv(clean_markers, OUT_CLEAN)
write_csv(marker_counts, OUT_COUNTS)
write_csv(top20_by_group, OUT_TOP20)
write_csv(top8_candidates, OUT_TOP8)
write_csv(description_audit, OUT_AUDIT)
write_tsv(summary_tbl, OUT_SUMMARY)

cat("\n============================================================\n")
cat("FINAL v2 MetaCyc pathway LEfSe parsing finished\n")
cat("============================================================\n\n")

cat("Description matching:\n")
print(description_summary, n = Inf, width = Inf)

cat("\nMarker counts by group:\n")
print(marker_counts, n = Inf, width = Inf)

cat("\nTop 20 pathways by group:\n")
print(top20_by_group, n = Inf, width = Inf)

cat("\nTop 8 representative candidates by group:\n")
print(top8_candidates, n = Inf, width = Inf)

cat("\nFiles written:\n")
cat("Final clean markers:  ", OUT_CLEAN, "\n")
cat("Marker counts:        ", OUT_COUNTS, "\n")
cat("Top20 by group:       ", OUT_TOP20, "\n")
cat("Top8 candidates:      ", OUT_TOP8, "\n")
cat("Description audit:    ", OUT_AUDIT, "\n")
cat("Summary:              ", OUT_SUMMARY, "\n\n")
cat("Done.\n")