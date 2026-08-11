
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_04_parse_and_select_pathways_progression5.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Parse the LEfSe .res output, map pathway identifiers to descriptions and
## functional modules, and select representative pathways per group.
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
  "LEfSe_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB"
)

RES_FILE <- file.path(
  LEFSE_DIR,
  "LEfSe_results_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.res"
)

FEATURE_MAP_FILE <- file.path(
  LEFSE_DIR,
  "MetaCyc_pathway_feature_name_map_oldstyle_currentdata_progression127_7KB.csv"
)

GROUP_LEVELS <- c("Polyp", "UC_remission", "UC_active", "Dysplasia", "CA")
LDA_CUTOFF <- 2

OUT_CLEAN <- file.path(
  LEFSE_DIR,
  "FINAL_LEfSe_clean_markers_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_COUNTS <- file.path(
  LEFSE_DIR,
  "FINAL_LEfSe_marker_counts_by_group_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_MODULE_SUMMARY <- file.path(
  LEFSE_DIR,
  "FINAL_LEfSe_marker_counts_by_group_module_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_REP_CANDIDATES <- file.path(
  LEFSE_DIR,
  "FINAL_representative_candidates_by_group_module_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_MAIN_PLOT <- file.path(
  LEFSE_DIR,
  "FINAL_main_plot_selected_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.csv"
)
OUT_SUMMARY <- file.path(
  LEFSE_DIR,
  "FINAL_summary_select_representative_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.tsv"
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
    str_detect(x, "peptidoglycan|muramoyl|cell wall|lipid iva|kdo|phospholipid|diacylglycerol|fatty acid|hopanoid|membrane|o-antigen|colanic") ~ "Cell envelope / lipid metabolism",
    str_detect(x, "glycolysis|gluconeogenesis|pentose phosphate|pyruvate|fermentation|acetate|lactate|butanoate|butyrate|starch|sucrose|glycogen|carbon|methylcitrate|ketogenesis|tca|glyoxylate") ~ "Central carbon / fermentation",
    str_detect(x, "purine|pyrimidine|nucleotide|nucleobase|ribonucleotide|guanosine|inosine|trna|adenosine|adenine|ump") ~ "Nucleotide / translation",
    str_detect(x, "methionine|cysteine|tryptophan|tyrosine|isoleucine|lysine|threonine|aspartate|arginine|ornithine|phenylalanine|amino acid|polyamine") ~ "Amino acid metabolism",
    str_detect(x, "sulfate|sulfur") ~ "Sulfur metabolism",
    str_detect(x, "nitrate|nitrogen|denitrification") ~ "Nitrogen metabolism",
    str_detect(x, "cobalamin|vitamin|cofactor|menaquinol|menaquinone|ubiquinol|heme|nad|thiamine|pyridoxal|queuosine") ~ "Cofactor / quinone metabolism",
    str_detect(x, "glycan|hexitol|galactarate|glucarate|fucose|rhamnose|galactose") ~ "Carbohydrate / glycan metabolism",
    TRUE ~ "Other"
  )
}

low_interpretability_flag <- function(pathway_id, description) {
  x <- str_to_lower(paste(pathway_id, description, sep = " "))
  
  str_detect(
    x,
    "photosynthetic|photorespiration|rubisco|calvin|plant|plastidic|yeast|saccharomyces|invertebrate|chlamydia|mycobacteria|mammalian"
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
    Feature_join = str_remove(Feature, "^f_"),
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
  mutate(q_value_BH = p.adjust(p_value, method = "BH"))

feature_map <- read_csv(FEATURE_MAP_FILE, show_col_types = FALSE) %>%
  mutate(
    Feature = as.character(Feature),
    Feature_join = str_remove(Feature, "^f_"),
    Pathway_ID = as.character(Pathway_ID),
    Description = as.character(Description)
  )

clean_markers <- res %>%
  left_join(
    feature_map %>% select(Feature_join, Pathway_ID, Description),
    by = "Feature_join"
  ) %>%
  mutate(
    Description_matched = !is.na(Pathway_ID),
    Pathway_ID = if_else(is.na(Pathway_ID) | Pathway_ID == "", Feature_join, Pathway_ID),
    Description = if_else(is.na(Description) | Description == "", Feature_join, Description),
    Description_display = display_text(Description),
    Pathway_module = classify_module(Pathway_ID, Description),
    Low_interpretability_flag = low_interpretability_flag(Pathway_ID, Description),
    Class = factor(Class, levels = GROUP_LEVELS),
    Main_plot_eligible = !Low_interpretability_flag & q_value_BH <= 0.10
  ) %>%
  select(
    Class,
    Pathway_ID,
    Description,
    Description_display,
    Pathway_module,
    Feature,
    Description_matched,
    Low_interpretability_flag,
    Main_plot_eligible,
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

module_summary <- clean_markers %>%
  count(Class, Pathway_module, name = "n_markers") %>%
  arrange(Class, desc(n_markers), Pathway_module)

representative_candidates <- clean_markers %>%
  filter(Main_plot_eligible) %>%
  group_by(Class, Pathway_module) %>%
  arrange(desc(abs_LDA), q_value_BH, p_value, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup() %>%
  arrange(Class, Pathway_module, desc(abs_LDA))

manual_priority <- tribble(
  ~Pathway_ID, ~Manual_priority_note,
  "PWY0-1477", "CA: ethanolamine utilization; biologically interpretable",
  "PWY-1042", "CA: glycolysis IV",
  "PWY-8178", "CA: pentose phosphate pathway II",
  "NONOXIPENT-PWY", "CA: pentose phosphate pathway I",
  "PWY-5100", "CA: pyruvate fermentation to acetate and lactate",
  "PWY0-1586", "CA: peptidoglycan maturation",
  "PWY-5667", "CA: CDP-diacylglycerol biosynthesis",
  "PWY0-1319", "CA: CDP-diacylglycerol biosynthesis",
  "PHOSLIPSYN-PWY", "CA: phospholipid biosynthesis",
  "TRNA-CHARGING-PWY", "CA: tRNA charging / translation",
  "THRESYN-PWY", "CA: threonine biosynthesis",
  "PWY-7228", "UC_active: guanosine nucleotide biosynthesis",
  "PWY-7208", "UC_active: pyrimidine nucleobase salvage",
  "DENOVOPURINE2-PWY", "UC_active: purine nucleotide biosynthesis",
  "NAGLIPASYN-PWY", "UC_active: lipid IVA biosynthesis",
  "PWY-8073", "UC_active: lipid IVA biosynthesis",
  "PWY-6467", "UC_active: Kdo transfer to lipid IVA",
  "PWY-6748", "UC_active: nitrate reduction / denitrification",
  "PWY-5861", "UC_active: demethylmenaquinol biosynthesis",
  "PWY66-399", "UC_active or Polyp: gluconeogenesis III",
  "PWY0-781", "Dysplasia: aspartate superpathway",
  "P163-PWY", "Dysplasia: lysine fermentation",
  "PWY66-367", "Dysplasia: ketogenesis",
  "PWY-7373", "Dysplasia: demethylmenaquinol biosynthesis",
  "PWY-7094", "Polyp: fatty acid salvage",
  "PWY66-409", "Polyp/nonCA-like: purine nucleotide salvage",
  "PWY0-162", "Polyp/nonCA-like: pyrimidine de novo biosynthesis",
  "SO4ASSIM-PWY", "Polyp/UC-active-like: assimilatory sulfate reduction",
  "PWY0-1241", "UC_remission: ADP-L-glycero-beta-D-manno-heptose biosynthesis",
  "PWY-7072", "UC_remission: hopanoid biosynthesis",
  "PWY-7399", "UC_remission: methylphosphonate degradation"
)

manual_selected <- clean_markers %>%
  inner_join(manual_priority, by = "Pathway_ID") %>%
  filter(Main_plot_eligible | Pathway_ID %in% c("PWY0-1477", "PWY-1042", "PWY-5100")) %>%
  arrange(Class, desc(abs_LDA), q_value_BH)

auto_fill <- representative_candidates %>%
  anti_join(manual_selected %>% select(Pathway_ID, Class), by = c("Pathway_ID", "Class")) %>%
  group_by(Class) %>%
  arrange(desc(abs_LDA), q_value_BH, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  mutate(Manual_priority_note = "Auto-added: top representative within group/module")

main_plot_selected <- bind_rows(manual_selected, auto_fill) %>%
  distinct(Class, Pathway_ID, .keep_all = TRUE) %>%
  mutate(
    Class = factor(Class, levels = GROUP_LEVELS),
    Plot_label = Description_display
  ) %>%
  arrange(Class, Pathway_module, desc(abs_LDA), q_value_BH)

write_csv(clean_markers, OUT_CLEAN)
write_csv(marker_counts, OUT_COUNTS)
write_csv(module_summary, OUT_MODULE_SUMMARY)
write_csv(representative_candidates, OUT_REP_CANDIDATES)
write_csv(main_plot_selected, OUT_MAIN_PLOT)

summary_tbl <- tibble(
  metric = c(
    "Total LEfSe markers retained",
    "Main plot selected pathways",
    "Description matched",
    "Description unmatched",
    "LDA cutoff",
    "Main plot eligibility q cutoff"
  ),
  value = c(
    nrow(clean_markers),
    nrow(main_plot_selected),
    sum(clean_markers$Description_matched),
    sum(!clean_markers$Description_matched),
    LDA_CUTOFF,
    0.10
  )
)
write_tsv(summary_tbl, OUT_SUMMARY)

cat("\n============================================================\n")
cat("5-group MetaCyc pathway representative selection finished\n")
cat("============================================================\n\n")

cat("Marker counts by group:\n")
print(marker_counts, n = Inf, width = Inf)

cat("\nMarker counts by group and module:\n")
print(module_summary, n = Inf, width = Inf)

cat("\nMain-plot selected pathways:\n")
print(
  main_plot_selected %>%
    select(Class, Pathway_ID, Description_display, Pathway_module, LDA, p_value, q_value_BH, Manual_priority_note),
  n = Inf,
  width = Inf
)

cat("\nFiles written:\n")
cat("Clean markers:              ", OUT_CLEAN, "\n")
cat("Marker counts:              ", OUT_COUNTS, "\n")
cat("Module summary:             ", OUT_MODULE_SUMMARY, "\n")
cat("Representative candidates:  ", OUT_REP_CANDIDATES, "\n")
cat("Main plot selected:         ", OUT_MAIN_PLOT, "\n")
cat("Summary:                    ", OUT_SUMMARY, "\n\n")
cat("Done.\n")