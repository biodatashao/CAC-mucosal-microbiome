
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_08_parse_LEfSe_pathways_CAC_vs_nonCAC.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Parse the CAC vs nonCAC LEfSe .res output into clean marker tables.
############################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

PROJECT_ROOT <- PROJECT_ROOT
MAIN_OUT <- OUTPUT_ROOT  ## 7KB rerun output root

LEFSE_DIR <- file.path(
  MAIN_OUT,
  "LEfSe_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB"
)

RES_FILE <- file.path(
  LEFSE_DIR,
  "LEfSe_results_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.res"
)

FEATURE_MAP_FILE <- file.path(
  LEFSE_DIR,
  "MetaCyc_pathway_feature_name_map_CA23_vs_nonCA23_oldstyle_7KB.csv"
)

GROUP_LEVELS <- c("nonCA", "CA")
LDA_CUTOFF <- 2

OUT_CLEAN <- file.path(
  LEFSE_DIR,
  "FINAL_LEfSe_clean_markers_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_COUNTS <- file.path(
  LEFSE_DIR,
  "FINAL_LEfSe_marker_counts_by_group_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_TOP40 <- file.path(
  LEFSE_DIR,
  "FINAL_LEfSe_top40_by_group_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_CA_ONLY <- file.path(
  LEFSE_DIR,
  "FINAL_CA_enriched_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_NONCA_ONLY <- file.path(
  LEFSE_DIR,
  "FINAL_nonCA_enriched_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.csv"
)
OUT_SUMMARY <- file.path(
  LEFSE_DIR,
  "FINAL_summary_LEfSe_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.tsv"
)

stopifnot(file.exists(RES_FILE))
stopifnot(file.exists(FEATURE_MAP_FILE))

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
    str_detect(x, "peptidoglycan|muramoyl|cell wall|lipid iva|kdo|phospholipid|diacylglycerol|fatty acid|hopanoid|membrane|colanic|o-antigen") ~ "Lipid / membrane / cell envelope",
    str_detect(x, "glycolysis|gluconeogenesis|pentose phosphate|pyruvate|fermentation|acetate|lactate|butanoate|butyrate|starch|sucrose|glycogen|carbon|methylcitrate|ketogenesis|tca|glyoxylate") ~ "Central carbon / fermentation",
    str_detect(x, "purine|pyrimidine|nucleotide|nucleobase|ribonucleotide|guanosine|inosine|trna|adenosine|adenine|ump") ~ "Nucleotide / translation",
    str_detect(x, "methionine|cysteine|tryptophan|tyrosine|isoleucine|lysine|threonine|aspartate|arginine|ornithine|phenylalanine|amino acid|polyamine") ~ "Amino acid metabolism",
    str_detect(x, "sulfate|sulfur") ~ "Sulfur metabolism",
    str_detect(x, "nitrate|nitrogen|denitrification") ~ "Nitrogen metabolism",
    str_detect(x, "cobalamin|vitamin|cofactor|menaquinol|menaquinone|ubiquinol|heme|nad|thiamine|pyridoxal") ~ "Cofactor / quinone metabolism",
    str_detect(x, "glycan|hexitol|galactarate|fucose|rhamnose|galactose") ~ "Carbohydrate / glycan metabolism",
    TRUE ~ "Other"
  )
}

res_raw <- read_tsv(
  RES_FILE,
  col_names = FALSE,
  show_col_types = FALSE,
  progress = FALSE
)

if (ncol(res_raw) < 5) {
  stop("LEfSe result file should contain at least 5 columns.")
}

res <- res_raw %>%
  select(1:5) %>%
  set_names(c("Feature", "log_max_mean", "Class", "LDA", "p_value")) %>%
  mutate(
    Feature = as.character(Feature),
    log_max_mean = suppressWarnings(as.numeric(log_max_mean)),
    Class = as.character(Class),
    LDA = suppressWarnings(as.numeric(LDA)),
    p_value = suppressWarnings(as.numeric(p_value)),
    abs_LDA = abs(LDA)
  ) %>%
  filter(
    !is.na(Feature),
    !is.na(Class),
    Class != "-",
    Class %in% GROUP_LEVELS,
    !is.na(LDA),
    abs_LDA >= LDA_CUTOFF
  ) %>%
  mutate(
    q_value_BH = p.adjust(p_value, method = "BH")
  )

feature_map <- read_csv(FEATURE_MAP_FILE, show_col_types = FALSE) %>%
  mutate(
    Feature = as.character(Feature),
    Pathway_ID = as.character(Pathway_ID),
    Description = as.character(Description)
  )

clean_markers <- res %>%
  left_join(feature_map, by = "Feature") %>%
  mutate(
    Description_matched = !is.na(Pathway_ID),
    Description = if_else(is.na(Description) | Description == "", Feature, Description),
    Pathway_ID = if_else(is.na(Pathway_ID) | Pathway_ID == "", Feature, Pathway_ID),
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
    Description_matched,
    log_max_mean,
    LDA,
    abs_LDA,
    p_value,
    q_value_BH
  ) %>%
  arrange(Class, desc(abs_LDA), p_value)

marker_counts <- clean_markers %>%
  count(Class, name = "n_markers") %>%
  complete(Class = factor(GROUP_LEVELS, levels = GROUP_LEVELS), fill = list(n_markers = 0)) %>%
  arrange(Class)

top40_by_group <- clean_markers %>%
  group_by(Class) %>%
  arrange(desc(abs_LDA), p_value, .by_group = TRUE) %>%
  slice_head(n = 40) %>%
  ungroup()

ca_only <- clean_markers %>%
  filter(Class == "CA") %>%
  arrange(desc(abs_LDA), p_value)

nonca_only <- clean_markers %>%
  filter(Class == "nonCA") %>%
  arrange(desc(abs_LDA), p_value)

summary_tbl <- tibble(
  metric = c(
    "Total LEfSe markers retained",
    "CA markers",
    "nonCA markers",
    "Description matched",
    "Description unmatched",
    "LDA cutoff"
  ),
  value = c(
    nrow(clean_markers),
    sum(clean_markers$Class == "CA"),
    sum(clean_markers$Class == "nonCA"),
    sum(clean_markers$Description_matched),
    sum(!clean_markers$Description_matched),
    LDA_CUTOFF
  )
)

write_csv(clean_markers, OUT_CLEAN)
write_csv(marker_counts, OUT_COUNTS)
write_csv(top40_by_group, OUT_TOP40)
write_csv(ca_only, OUT_CA_ONLY)
write_csv(nonca_only, OUT_NONCA_ONLY)
write_tsv(summary_tbl, OUT_SUMMARY)

cat("\n============================================================\n")
cat("CA vs nonCA MetaCyc pathway LEfSe parsing finished\n")
cat("============================================================\n\n")

cat("Marker counts by group:\n")
print(marker_counts, n = Inf, width = Inf)

cat("\nTop 40 pathways by group:\n")
print(top40_by_group, n = Inf, width = Inf)

cat("\nTop CA-enriched pathways:\n")
print(ca_only %>% slice_head(n = 40), n = Inf, width = Inf)

cat("\nFiles written:\n")
cat("Final clean markers: ", OUT_CLEAN, "\n")
cat("Marker counts:       ", OUT_COUNTS, "\n")
cat("Top40 by group:      ", OUT_TOP40, "\n")
cat("CA enriched:         ", OUT_CA_ONLY, "\n")
cat("nonCA enriched:      ", OUT_NONCA_ONLY, "\n")
cat("Summary:             ", OUT_SUMMARY, "\n\n")
cat("Done.\n")