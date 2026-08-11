


############################################################
## 03_07_LEfSe_CAC_vs_nonCAC_genus.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
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
## 1. Paths
############################################################

PROJECT_DIR <- PROJECT_ROOT

RUN_DIR <- file.path(
  PROJECT_DIR,
  "output/analysis"
)

INPUT_DIR <- file.path(
  RUN_DIR,
  "00_clean_data",
  "CA23_nonCA23_paired"
)

OUT_DIR <- file.path(
  RUN_DIR,
  "02_Figure3_taxa_LEfSe",
  "LEfSe_CA23_vs_nonCA23_7KB_GENUS_ONLY"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

COUNT_FILE <- file.path(
  INPUT_DIR,
  "asv_count_7KB_CA23_nonCA23_paired.tsv"
)

TAX_FILE <- file.path(
  INPUT_DIR,
  "taxonomy_7KB_CA23_nonCA23_paired.tsv"
)

META_FILE <- file.path(
  INPUT_DIR,
  "metadata_7KB_CA23_nonCA23_paired.tsv"
)

CONDA_SH <- "/opt/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV <- "lefse"

############################################################
## 2. Parameters
## Exactly inherited from original 06 script
############################################################

GROUP_ORDER <- c(
  "nonCA",
  "CA"
)

MIN_PREVALENCE_N <- 5
MIN_MEAN_RA <- 1e-4
LDA_CUTOFF <- 2.0

############################################################
## 3. Output files
############################################################

OUT_SAMPLE_AUDIT <- file.path(
  OUT_DIR,
  "sample_presence_audit_CA23_vs_nonCA23_7KB.tsv"
)

OUT_METADATA <- file.path(
  OUT_DIR,
  "metadata_CA23_vs_nonCA23_unpaired_7KB.tsv"
)

OUT_GENUS_RA_ALL <- file.path(
  OUT_DIR,
  "genus_relative_abundance_CA23_vs_nonCA23_all_genus_only_7KB.tsv"
)

OUT_GENUS_RA_LEFSE <- file.path(
  OUT_DIR,
  "genus_relative_abundance_CA23_vs_nonCA23_LEfSe_features_GENUS_ONLY_7KB.tsv"
)

OUT_GENUS_FILTER <- file.path(
  OUT_DIR,
  "genus_filtering_summary_CA23_vs_nonCA23_GENUS_ONLY_7KB.tsv"
)

OUT_CANDIDATE_AUDIT <- file.path(
  OUT_DIR,
  "candidate_marker_genera_audit_CA23_vs_nonCA23_GENUS_ONLY_7KB.tsv"
)

LEFSE_INPUT_FILE <- file.path(
  OUT_DIR,
  "LEfSe_input_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.tsv"
)

LEFSE_FORMATTED_FILE <- file.path(
  OUT_DIR,
  "LEfSe_formatted_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.in"
)

LEFSE_RESULT_FILE <- file.path(
  OUT_DIR,
  "LEfSe_results_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.res"
)

OUT_LEFSE_ALL <- file.path(
  OUT_DIR,
  "LEfSe_all_results_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.csv"
)

OUT_LEFSE_MARKERS <- file.path(
  OUT_DIR,
  "LEfSe_clean_markers_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.csv"
)

OUT_LEFSE_CANDIDATES <- file.path(
  OUT_DIR,
  "LEfSe_candidate_results_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.csv"
)

OUT_SUMMARY <- file.path(
  OUT_DIR,
  "summary_LEfSe_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.tsv"
)

############################################################
## 4. Helpers
############################################################

stop_if_missing <- function(files) {
  
  missing_files <- files[
    !file.exists(files)
  ]
  
  if (length(missing_files) > 0) {
    stop(
      "Missing required file(s):\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  }
}

detect_taxonomy_col <- function(df) {
  
  candidates <- c(
    "Taxon",
    "taxon",
    "Taxonomy",
    "taxonomy",
    "Feature.Taxon",
    "Consensus.Lineage",
    "Consensus.Lineage."
  )
  
  out <- candidates[
    candidates %in% colnames(df)
  ][1]
  
  if (!is.na(out)) {
    return(out)
  }
  
  genus_candidates <- c(
    "Genus",
    "genus",
    "Genus_final",
    "Genus_clean",
    "Genus_clean_LEfSe"
  )
  
  out2 <- genus_candidates[
    genus_candidates %in% colnames(df)
  ][1]
  
  if (!is.na(out2)) {
    return(out2)
  }
  
  stop(
    "Cannot detect taxonomy/genus column. Columns are:\n",
    paste(
      colnames(df),
      collapse = ", "
    )
  )
}

extract_genus_strict <- function(x) {
  
  x <- as.character(x)
  x <- str_trim(x)
  
  genus <- rep(
    NA_character_,
    length(x)
  )
  
  has_g <- str_detect(
    x,
    "g__"
  )
  
  if (any(has_g, na.rm = TRUE)) {
    
    g_part <- str_extract(
      x[has_g],
      "g__[^;]+"
    )
    
    g_part <- str_replace(
      g_part,
      "^g__",
      ""
    )
    
    g_part <- str_trim(
      g_part
    )
    
    g_part[
      g_part == "" |
        is.na(g_part)
    ] <- NA_character_
    
    genus[has_g] <- g_part
  }
  
  no_g_idx <- which(
    !has_g |
      is.na(has_g)
  )
  
  if (length(no_g_idx) > 0) {
    
    candidate <- x[
      no_g_idx
    ]
    
    candidate <- str_replace(
      candidate,
      "^g__",
      ""
    )
    
    candidate <- str_trim(
      candidate
    )
    
    looks_higher_tax <- str_detect(
      candidate,
      regex(
        "(^|;)\\s*[dkpcofs]__|k__|p__|c__|o__|f__|s__",
        ignore_case = TRUE
      )
    )
    
    candidate[
      looks_higher_tax
    ] <- NA_character_
    
    candidate[
      candidate == ""
    ] <- NA_character_
    
    genus[
      no_g_idx
    ] <- candidate
  }
  
  genus
}

is_bad_genus <- function(x) {
  
  x0 <- tolower(
    as.character(x)
  )
  
  is.na(x0) |
    x0 == "" |
    x0 == "na" |
    x0 == "unassigned" |
    x0 == "unclassified" |
    x0 == "uncultured" |
    x0 == "unknown" |
    x0 == "metagenome" |
    x0 == "norank" |
    x0 == "no_rank" |
    x0 == "ambiguous" |
    str_detect(x0, "uncultured") |
    str_detect(x0, "unclassified") |
    str_detect(x0, "unknown") |
    str_detect(x0, "metagenome") |
    str_detect(x0, "norank") |
    str_detect(x0, "no_rank") |
    str_detect(x0, "ambiguous") |
    str_detect(x0, "^[dkpcofs]__") |
    str_detect(x0, "k__|p__|c__|o__|f__|s__")
}

make_LEfSe_safe_name <- function(x) {
  
  x <- as.character(x)
  
  x <- str_replace_all(
    x,
    "-",
    "_"
  )
  
  x <- str_replace_all(
    x,
    "\\s+",
    "_"
  )
  
  x <- str_replace_all(
    x,
    "[^A-Za-z0-9_\\.]",
    "_"
  )
  
  x
}

############################################################
## 5. Read data
############################################################

stop_if_missing(
  c(
    COUNT_FILE,
    TAX_FILE,
    META_FILE,
    CONDA_SH
  )
)

count_df <- read_tsv(
  COUNT_FILE,
  show_col_types = FALSE
)

tax_df <- read_tsv(
  TAX_FILE,
  show_col_types = FALSE
)

meta_df <- read_tsv(
  META_FILE,
  show_col_types = FALSE
)

if (!"SampleID" %in% names(meta_df)) {
  stop(
    "Metadata does not contain SampleID."
  )
}

if (!"PairID" %in% names(meta_df)) {
  stop(
    "Metadata does not contain PairID."
  )
}

if (!"CA_vs_nonCA_explicit" %in% names(meta_df)) {
  stop(
    "Metadata does not contain CA_vs_nonCA_explicit."
  )
}

############################################################
## 6. Fixed 23 vs 23 metadata
############################################################

meta_analysis <- meta_df %>%
  transmute(
    PairID = as.character(PairID),
    SampleID = as.character(SampleID),
    Group = as.character(
      CA_vs_nonCA_explicit
    )
  ) %>%
  filter(
    Group %in% GROUP_ORDER
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = GROUP_ORDER
    ),
    
    LEfSe_class = as.character(
      Group
    )
  ) %>%
  arrange(
    Group,
    SampleID
  )

group_counts <- meta_analysis %>%
  count(
    Group,
    name = "n"
  ) %>%
  arrange(
    Group
  )

cat(
  "\nCA23 vs nonCA23 group counts:\n"
)

print(
  group_counts,
  n = Inf,
  width = Inf
)

if (nrow(meta_analysis) != 46) {
  stop(
    "Expected 46 samples, found ",
    nrow(meta_analysis),
    "."
  )
}

observed_counts <- setNames(
  group_counts$n,
  as.character(
    group_counts$Group
  )
)

if (
  observed_counts["nonCA"] != 23 ||
  observed_counts["CA"] != 23
) {
  stop(
    "Expected nonCA=23 and CA=23."
  )
}

pair_audit <- meta_analysis %>%
  count(
    PairID,
    name = "n_samples"
  )

if (
  nrow(pair_audit) != 23 ||
  any(pair_audit$n_samples != 2)
) {
  stop(
    "PairID audit failed: expected 23 pairs with 2 samples per pair."
  )
}

############################################################
## 7. Count matrix
############################################################

asv_col_count <- names(
  count_df
)[1]

count_mat <- count_df %>%
  column_to_rownames(
    asv_col_count
  ) %>%
  as.matrix()

storage.mode(
  count_mat
) <- "numeric"

############################################################
## 8. Taxonomy -> strict genus
############################################################

asv_col_tax <- names(
  tax_df
)[1]

taxonomy_col <- detect_taxonomy_col(
  tax_df
)

tax_df <- tax_df %>%
  rename(
    ASV = all_of(
      asv_col_tax
    )
  ) %>%
  mutate(
    Taxonomy_source_for_genus =
      as.character(
        .data[[taxonomy_col]]
      ),
    
    Genus_strict =
      extract_genus_strict(
        Taxonomy_source_for_genus
      ),
    
    Genus_is_valid =
      !is_bad_genus(
        Genus_strict
      )
  )

############################################################
## 9. Sample audit
############################################################

sample_audit <- meta_analysis %>%
  mutate(
    Present_in_count =
      SampleID %in% colnames(
        count_mat
      )
  )

write_tsv(
  sample_audit,
  OUT_SAMPLE_AUDIT
)

if (
  any(
    !sample_audit$Present_in_count
  )
) {
  
  print(
    sample_audit %>%
      filter(
        !Present_in_count
      ),
    n = Inf,
    width = Inf
  )
  
  stop(
    "Some requested CA/nonCA samples are missing from 7KB count table."
  )
}

write_tsv(
  meta_analysis,
  OUT_METADATA
)

############################################################
## 10. Subset count table
############################################################

sample_keep <- meta_analysis$SampleID

count_sub <- count_mat[
  ,
  sample_keep,
  drop = FALSE
]

count_sub <- count_sub[
  rowSums(
    count_sub
  ) > 0,
  ,
  drop = FALSE
]

############################################################
## 11. Align ASVs
############################################################

common_asvs <- intersect(
  rownames(
    count_sub
  ),
  tax_df$ASV
)

count_sub <- count_sub[
  common_asvs,
  ,
  drop = FALSE
]

tax_sub <- tax_df %>%
  filter(
    ASV %in% common_asvs
  ) %>%
  arrange(
    match(
      ASV,
      rownames(
        count_sub
      )
    )
  )

stopifnot(
  identical(
    tax_sub$ASV,
    rownames(
      count_sub
    )
  )
)

stopifnot(
  identical(
    colnames(
      count_sub
    ),
    meta_analysis$SampleID
  )
)

cat(
  "\nASV alignment:\n"
)

cat(
  "ASVs retained in selected 46 samples: ",
  nrow(count_sub),
  "\n",
  sep = ""
)

cat(
  "Common ASVs with taxonomy: ",
  length(common_asvs),
  "\n",
  sep = ""
)

############################################################
## 12. Collapse ASVs to strict genus
############################################################

tax_genus <- tax_sub %>%
  filter(
    Genus_is_valid
  ) %>%
  mutate(
    Genus = Genus_strict
  )

if (nrow(tax_genus) == 0) {
  stop(
    "No valid genus-level ASVs after strict genus filtering."
  )
}

count_genus_input <- count_sub[
  tax_genus$ASV,
  ,
  drop = FALSE
]

genus_count <- rowsum(
  count_genus_input,
  group = tax_genus$Genus,
  reorder = FALSE
)

genus_count <- genus_count[
  rowSums(
    genus_count
  ) > 0,
  ,
  drop = FALSE
]

sample_depth_genus <- colSums(
  genus_count
)

if (any(sample_depth_genus <= 0)) {
  stop(
    "At least one sample has zero genus-level reads."
  )
}

genus_ra <- sweep(
  genus_count,
  2,
  sample_depth_genus,
  "/"
)

genus_ra[
  is.na(
    genus_ra
  )
] <- 0

cat(
  "\nGenus collapsing:\n"
)

cat(
  "Valid genus-level ASVs: ",
  nrow(tax_genus),
  "\n",
  sep = ""
)

cat(
  "Genus features: ",
  nrow(genus_ra),
  "\n",
  sep = ""
)

############################################################
## 13. Feature filtering for LEfSe
## Exact original thresholds
############################################################

prevalence_n <- rowSums(
  genus_ra > 0
)

mean_ra <- rowMeans(
  genus_ra
)

genus_filter_table <- tibble(
  Genus = rownames(
    genus_ra
  ),
  
  prevalence_n =
    as.integer(
      prevalence_n
    ),
  
  mean_RA =
    as.numeric(
      mean_ra
    ),
  
  mean_percent =
    mean_RA * 100,
  
  keep_for_LEfSe =
    prevalence_n >= MIN_PREVALENCE_N &
    mean_ra >= MIN_MEAN_RA
) %>%
  arrange(
    desc(
      keep_for_LEfSe
    ),
    desc(
      mean_RA
    )
  )

write_tsv(
  genus_filter_table,
  OUT_GENUS_FILTER
)

keep_genera <- genus_filter_table %>%
  filter(
    keep_for_LEfSe
  ) %>%
  pull(
    Genus
  )

if (length(keep_genera) == 0) {
  stop(
    "Strict genus-only filtering retained 0 genera."
  )
}

genus_ra_filtered <- genus_ra[
  keep_genera,
  ,
  drop = FALSE
]

############################################################
## 14. Save abundance matrices
############################################################

write_tsv(
  genus_ra %>%
    as.data.frame(
      check.names = FALSE
    ) %>%
    rownames_to_column(
      "Genus"
    ),
  OUT_GENUS_RA_ALL
)

write_tsv(
  genus_ra_filtered %>%
    as.data.frame(
      check.names = FALSE
    ) %>%
    rownames_to_column(
      "Genus"
    ),
  OUT_GENUS_RA_LEFSE
)

############################################################
## 15. Candidate marker audit
## Same list as original 06
############################################################

candidate_genera <- c(
  "UCG-005",
  "UCG_005",
  "Desulfovibrio",
  "Dorea",
  "Mediterraneibacter",
  "Monoglobus",
  "Lactococcus",
  "Atopostipes",
  "Bacteroides",
  "Fusobacterium",
  "Escherichia-Shigella",
  "Escherichia_Shigella",
  "Enterococcus",
  "Streptococcus",
  "Blautia",
  "Faecalibacterium",
  "Romboutsia"
)

candidate_audit <- genus_filter_table %>%
  mutate(
    Genus_alt1 = str_replace_all(
      Genus,
      "-",
      "_"
    ),
    
    Genus_alt2 = str_replace_all(
      Genus,
      "_",
      "-"
    )
  ) %>%
  filter(
    Genus %in% candidate_genera |
      Genus_alt1 %in% candidate_genera |
      Genus_alt2 %in% candidate_genera
  ) %>%
  mutate(
    kept_for_LEfSe =
      keep_for_LEfSe
  ) %>%
  select(
    Genus,
    prevalence_n,
    mean_RA,
    mean_percent,
    kept_for_LEfSe
  ) %>%
  arrange(
    desc(
      mean_RA
    )
  )

write_tsv(
  candidate_audit,
  OUT_CANDIDATE_AUDIT
)

############################################################
## 16. Build LEfSe input
############################################################

lefse_feature_mat <-
  genus_ra_filtered * 100

rownames(
  lefse_feature_mat
) <- make_LEfSe_safe_name(
  rownames(
    lefse_feature_mat
  )
)

class_row <- as.character(
  meta_analysis$LEfSe_class
)

names(
  class_row
) <- meta_analysis$SampleID

lefse_input <- rbind(
  Class = class_row,
  lefse_feature_mat
)

write.table(
  lefse_input,
  file = LEFSE_INPUT_FILE,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = FALSE
)

############################################################
## 17. LEfSe input audit
############################################################

cat(
  "\nChecking LEfSe input:\n"
)

cat(
  "Input file: ",
  LEFSE_INPUT_FILE,
  "\n",
  sep = ""
)

print(
  system2(
    "wc",
    args = c(
      "-l",
      LEFSE_INPUT_FILE
    ),
    stdout = TRUE
  )
)

nf_cmd <- paste0(
  "awk -F'\\t' '{print NF}' ",
  shQuote(
    LEFSE_INPUT_FILE
  ),
  " | sort | uniq -c"
)

print(
  system(
    nf_cmd,
    intern = TRUE
  )
)

############################################################
## 18. Run LEfSe
## Same conda strategy as original
############################################################

cat(
  "\nRunning LEfSe inside R...\n"
)

unlink(
  c(
    LEFSE_FORMATTED_FILE,
    LEFSE_RESULT_FILE
  )
)

run_cmd <- paste(
  "unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_PROFILE R_ENVIRON &&",
  "source", shQuote(CONDA_SH), "&&",
  "conda activate", CONDA_ENV, "&&",
  "export R_HOME=$CONDA_PREFIX/lib/R &&",
  "export R_DEFAULT_PACKAGES='utils,stats,graphics,grDevices,datasets,methods' &&",
  "echo CONDA_PREFIX=$CONDA_PREFIX &&",
  "echo R_HOME=$R_HOME &&",
  "which lefse_format_input.py &&",
  "which lefse_run.py &&",
  "lefse_format_input.py",
  shQuote(LEFSE_INPUT_FILE),
  shQuote(LEFSE_FORMATTED_FILE),
  "-c 1 -o 1000000 &&",
  "lefse_run.py",
  shQuote(LEFSE_FORMATTED_FILE),
  shQuote(LEFSE_RESULT_FILE),
  "-l",
  as.character(
    LDA_CUTOFF
  )
)

lefse_log <- system2(
  command = "bash",
  args = c(
    "-lc",
    shQuote(
      run_cmd
    )
  ),
  stdout = TRUE,
  stderr = TRUE
)

cat(
  "\nLEfSe log:\n"
)

cat(
  paste(
    lefse_log,
    collapse = "\n"
  ),
  "\n"
)

cat(
  "\nFormatted file exists: ",
  file.exists(
    LEFSE_FORMATTED_FILE
  ),
  "\n",
  sep = ""
)

cat(
  "Result file exists: ",
  file.exists(
    LEFSE_RESULT_FILE
  ),
  "\n",
  sep = ""
)

if (
  !file.exists(
    LEFSE_RESULT_FILE
  )
) {
  stop(
    "LEfSe result was not generated. Check LEfSe log above."
  )
}

############################################################
## 19. Read LEfSe results
############################################################

lefse_all <- read_tsv(
  LEFSE_RESULT_FILE,
  col_names = c(
    "Genus",
    "log_max_mean",
    "Class",
    "LDA",
    "p_value"
  ),
  show_col_types = FALSE
) %>%
  mutate(
    Class = na_if(
      Class,
      ""
    ),
    
    Class = if_else(
      Class == "-",
      NA_character_,
      Class
    ),
    
    LDA = suppressWarnings(
      as.numeric(
        LDA
      )
    ),
    
    p_value =
      suppressWarnings(
        as.numeric(
          p_value
        )
      ),
    
    marker_status = if_else(
      !is.na(Class) &
        !is.na(LDA),
      "LEfSe_marker",
      "Not_marker"
    )
  )

lefse_markers <- lefse_all %>%
  filter(
    marker_status ==
      "LEfSe_marker"
  ) %>%
  arrange(
    factor(
      Class,
      levels = GROUP_ORDER
    ),
    desc(
      LDA
    )
  )

############################################################
## 20. Candidate result subset
############################################################

candidate_pattern <- paste(
  c(
    "UCG_005",
    "UCG-005",
    "Desulfovibrio",
    "Dorea",
    "Mediterraneibacter",
    "Monoglobus",
    "Lactococcus",
    "Atopostipes",
    "Bacteroides",
    "Fusobacterium",
    "Escherichia_Shigella",
    "Escherichia-Shigella",
    "Enterococcus",
    "Streptococcus",
    "Blautia",
    "Faecalibacterium",
    "Romboutsia"
  ),
  collapse = "|"
)

lefse_candidates <- lefse_all %>%
  filter(
    str_detect(
      Genus,
      regex(
        candidate_pattern,
        ignore_case = TRUE
      )
    )
  ) %>%
  arrange(
    desc(
      marker_status ==
        "LEfSe_marker"
    ),
    factor(
      Class,
      levels = GROUP_ORDER
    ),
    desc(
      LDA
    )
  )

write_csv(
  lefse_all,
  OUT_LEFSE_ALL
)

write_csv(
  lefse_markers,
  OUT_LEFSE_MARKERS
)

write_csv(
  lefse_candidates,
  OUT_LEFSE_CANDIDATES
)

############################################################
## 21. Summary
############################################################

marker_counts <- lefse_markers %>%
  count(
    Class,
    name = "n_markers"
  ) %>%
  arrange(
    factor(
      Class,
      levels = GROUP_ORDER
    )
  )

summary_tbl <- tibble(
  Item = c(
    "Analysis",
    "Input count table",
    "Input taxonomy table",
    "Input metadata table",
    "Comparison",
    "Pairing",
    "Samples",
    "Group counts",
    "Pairs audited",
    "ASVs in selected 46-sample subset",
    "ASVs with valid strict genus",
    "Genus-only features before LEfSe filtering",
    "Genus-only features kept for LEfSe",
    "Minimum prevalence n",
    "Minimum mean relative abundance",
    "LDA cutoff",
    "Predefined composition contaminant removal",
    "LEfSe input",
    "LEfSe result",
    "LEfSe all result CSV",
    "LEfSe clean marker CSV",
    "LEfSe candidate CSV",
    "Output directory"
  ),
  
  Value = c(
    "7KB CA23 vs nonCA23 unpaired strict genus-only LEfSe",
    COUNT_FILE,
    TAX_FILE,
    META_FILE,
    "nonCA vs CA",
    "Unpaired; PairID retained for audit only",
    as.character(
      nrow(
        meta_analysis
      )
    ),
    paste(
      group_counts$Group,
      group_counts$n,
      sep = "=",
      collapse = "; "
    ),
    as.character(
      nrow(
        pair_audit
      )
    ),
    as.character(
      nrow(
        count_sub
      )
    ),
    as.character(
      nrow(
        tax_genus
      )
    ),
    as.character(
      nrow(
        genus_ra
      )
    ),
    as.character(
      nrow(
        genus_ra_filtered
      )
    ),
    as.character(
      MIN_PREVALENCE_N
    ),
    as.character(
      MIN_MEAN_RA
    ),
    as.character(
      LDA_CUTOFF
    ),
    "No; original 06 LEfSe script did not remove Mesobacillus/Fictibacillus/Phaselicystis",
    LEFSE_INPUT_FILE,
    LEFSE_RESULT_FILE,
    OUT_LEFSE_ALL,
    OUT_LEFSE_MARKERS,
    OUT_LEFSE_CANDIDATES,
    OUT_DIR
  )
)

write_tsv(
  summary_tbl,
  OUT_SUMMARY
)

############################################################
## 22. Print results
############################################################

cat(
  "\n============================================================\n"
)

cat(
  "7KB CA23 vs nonCA23 unpaired genus-only LEfSe finished\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Group counts:\n"
)

print(
  group_counts,
  n = Inf,
  width = Inf
)

cat(
  "\nPair audit:\n"
)

print(
  pair_audit,
  n = Inf,
  width = Inf
)

cat(
  "\nGenus filtering summary:\n"
)

print(
  genus_filter_table %>%
    summarise(
      all_genus_only_features = n(),
      
      kept_for_LEfSe =
        sum(
          keep_for_LEfSe
        ),
      
      not_kept =
        sum(
          !keep_for_LEfSe
        )
    ),
  n = Inf,
  width = Inf
)

cat(
  "\nCandidate genus filtering audit:\n"
)

print(
  candidate_audit,
  n = Inf,
  width = Inf
)

cat(
  "\nAll LEfSe markers:\n"
)

print(
  lefse_markers,
  n = Inf,
  width = Inf
)

cat(
  "\nCandidate LEfSe results:\n"
)

print(
  lefse_candidates,
  n = Inf,
  width = Inf
)

cat(
  "\nMarker counts by class:\n"
)

print(
  marker_counts,
  n = Inf,
  width = Inf
)

cat(
  "\nFiles written:\n"
)

cat(
  "Sample audit:      ",
  OUT_SAMPLE_AUDIT,
  "\n"
)

cat(
  "Metadata:          ",
  OUT_METADATA,
  "\n"
)

cat(
  "Genus RA all:      ",
  OUT_GENUS_RA_ALL,
  "\n"
)

cat(
  "Genus RA LEfSe:    ",
  OUT_GENUS_RA_LEFSE,
  "\n"
)

cat(
  "LEfSe input:       ",
  LEFSE_INPUT_FILE,
  "\n"
)

cat(
  "LEfSe result:      ",
  LEFSE_RESULT_FILE,
  "\n"
)

cat(
  "Clean markers CSV: ",
  OUT_LEFSE_MARKERS,
  "\n"
)

cat(
  "Candidate CSV:     ",
  OUT_LEFSE_CANDIDATES,
  "\n"
)

cat(
  "Summary:           ",
  OUT_SUMMARY,
  "\n"
)

cat(
  "\nDone.\n"
)