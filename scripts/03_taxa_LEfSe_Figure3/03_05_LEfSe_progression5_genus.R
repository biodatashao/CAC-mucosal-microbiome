


############################################################
## 03_05_LEfSe_progression5_genus.R
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
## Paths
############################################################

PROJECT_DIR <- PROJECT_ROOT

RUN_DIR <- file.path(
  PROJECT_DIR,
  "output/analysis"
)

INPUT_DIR <- file.path(
  RUN_DIR,
  "00_clean_data",
  "progression127"
)

OUT_DIR <- file.path(
  RUN_DIR,
  "02_Figure3_taxa_LEfSe",
  "LEfSe_progression5_7KB_GENUS_ONLY"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

COUNT_FILE <- file.path(
  INPUT_DIR,
  "asv_count_7KB_progression127.tsv"
)

TAX_FILE <- file.path(
  INPUT_DIR,
  "taxonomy_7KB_progression127.tsv"
)

META_FILE <- file.path(
  INPUT_DIR,
  "metadata_7KB_progression127.tsv"
)

CONDA_SH <- "/opt/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV <- "lefse"

############################################################
## Parameters
## EXACTLY inherited from original LEfSe workflow
############################################################

LDA_CUTOFF <- 2.0

MIN_PREVALENCE_N <- 5

MIN_MEAN_RA <- 1e-4

GROUP_ORDER <- c(
  "Polyp",
  "UC_remission",
  "UC_active",
  "Dysplasia",
  "CA"
)

############################################################
## Helpers
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

detect_sample_col <- function(df) {
  
  candidates <- c(
    "SampleID",
    "Sample_ID",
    "sample_id",
    "sample",
    "Sample",
    "#SampleID"
  )
  
  out <- candidates[
    candidates %in% colnames(df)
  ][1]
  
  if (is.na(out)) {
    stop(
      "Cannot detect sample column. Columns are:\n",
      paste(
        colnames(df),
        collapse = ", "
      )
    )
  }
  
  out
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
## Read input
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

asv_col_count <- colnames(
  count_df
)[1]

asv_col_tax <- colnames(
  tax_df
)[1]

sample_col <- detect_sample_col(
  meta_df
)

taxonomy_col <- detect_taxonomy_col(
  tax_df
)

############################################################
## Count matrix
############################################################

count_mat <- count_df %>%
  column_to_rownames(
    asv_col_count
  ) %>%
  as.matrix()

storage.mode(
  count_mat
) <- "numeric"

############################################################
## Taxonomy
############################################################

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
## Metadata
############################################################

if (sample_col != "SampleID") {
  
  meta_df <- meta_df %>%
    rename(
      SampleID = all_of(
        sample_col
      )
    )
}

if (!"Progression5" %in% colnames(meta_df)) {
  stop(
    "Metadata does not contain Progression5."
  )
}

meta_prog <- meta_df %>%
  filter(
    Progression5 %in% GROUP_ORDER
  ) %>%
  mutate(
    Progression5 = factor(
      Progression5,
      levels = GROUP_ORDER
    ),
    
    LEfSe_class = as.character(
      Progression5
    )
  )

############################################################
## Strict progression127 audit
############################################################

group_counts <- meta_prog %>%
  count(
    Progression5,
    name = "n"
  ) %>%
  arrange(
    Progression5
  )

cat(
  "\nProgression127 group counts:\n"
)

print(
  group_counts,
  n = Inf,
  width = Inf
)

if (nrow(meta_prog) != 127) {
  stop(
    "Expected progression127 n=127, found ",
    nrow(meta_prog),
    "."
  )
}

expected_counts <- c(
  Polyp = 26,
  UC_remission = 36,
  UC_active = 25,
  Dysplasia = 17,
  CA = 23
)

observed_counts <- setNames(
  group_counts$n,
  as.character(
    group_counts$Progression5
  )
)

if (
  !identical(
    as.integer(
      observed_counts[
        names(expected_counts)
      ]
    ),
    as.integer(
      expected_counts
    )
  )
) {
  stop(
    "Progression5 counts do not match expected 26/36/25/17/23."
  )
}

############################################################
## Arrange samples
##
## Original LEfSe script used:
## arrange(Progression5, SampleID)
############################################################

meta_prog <- meta_prog %>%
  arrange(
    Progression5,
    SampleID
  )

sample_audit <- meta_prog %>%
  transmute(
    SampleID,
    Progression5,
    Present_in_count =
      SampleID %in% colnames(
        count_mat
      )
  )

write_tsv(
  sample_audit,
  file.path(
    OUT_DIR,
    "progression127_sample_presence_audit_7KB.tsv"
  )
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
    "Some progression127 samples are missing from the 7KB count table."
  )
}

############################################################
## Subset count table
############################################################

count_prog <- count_mat[
  ,
  meta_prog$SampleID,
  drop = FALSE
]

count_prog <- count_prog[
  rowSums(
    count_prog
  ) > 0,
  ,
  drop = FALSE
]

############################################################
## Align ASVs
############################################################

common_asvs <- intersect(
  rownames(
    count_prog
  ),
  tax_df$ASV
)

count_prog <- count_prog[
  common_asvs,
  ,
  drop = FALSE
]

tax_prog <- tax_df %>%
  filter(
    ASV %in% common_asvs
  ) %>%
  arrange(
    match(
      ASV,
      rownames(
        count_prog
      )
    )
  )

stopifnot(
  identical(
    rownames(
      count_prog
    ),
    tax_prog$ASV
  )
)

stopifnot(
  identical(
    colnames(
      count_prog
    ),
    meta_prog$SampleID
  )
)

write_tsv(
  group_counts,
  file.path(
    OUT_DIR,
    "progression127_group_counts_7KB.tsv"
  )
)

write_tsv(
  meta_prog,
  file.path(
    OUT_DIR,
    "metadata_progression5_7KB.tsv"
  )
)

############################################################
## Collapse ASVs to strict genus
############################################################

taxonomy_audit <- tibble(
  Item = c(
    "ASVs in progression127 7KB subset",
    "ASVs with valid strict genus",
    "ASVs excluded due unresolved/high-level taxonomy",
    "Taxonomy column used"
  ),
  Value = c(
    as.character(
      nrow(
        tax_prog
      )
    ),
    
    as.character(
      sum(
        tax_prog$Genus_is_valid
      )
    ),
    
    as.character(
      sum(
        !tax_prog$Genus_is_valid
      )
    ),
    
    taxonomy_col
  )
)

write_tsv(
  taxonomy_audit,
  file.path(
    OUT_DIR,
    "taxonomy_genus_only_audit_progression127_7KB.tsv"
  )
)

tax_genus <- tax_prog %>%
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

count_genus_input <- count_prog[
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

############################################################
## Genus relative abundance
############################################################

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

############################################################
## Original LEfSe feature filtering
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
  file.path(
    OUT_DIR,
    "genus_filtering_summary_progression127_7KB_GENUS_ONLY.tsv"
  )
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
## Save abundance matrices
############################################################

write_tsv(
  genus_ra %>%
    as.data.frame(
      check.names = FALSE
    ) %>%
    rownames_to_column(
      "Genus"
    ),
  file.path(
    OUT_DIR,
    "genus_relative_abundance_progression127_all_genus_only_7KB.tsv"
  )
)

write_tsv(
  genus_ra_filtered %>%
    as.data.frame(
      check.names = FALSE
    ) %>%
    rownames_to_column(
      "Genus"
    ),
  file.path(
    OUT_DIR,
    "genus_relative_abundance_progression127_LEfSe_features_GENUS_ONLY_7KB.tsv"
  )
)

############################################################
## Candidate audit
##
## Same candidate list as original script
############################################################

candidate_genera <- c(
  "UCG-005",
  "UCG_005",
  "Desulfovibrio",
  "Candidatus_Soleaferrea",
  "Mediterraneibacter",
  "Dorea",
  "Monoglobus",
  "Lactococcus",
  "Atopostipes",
  "Peptoclostridium",
  "f__Eubacterium__eligens_group",
  "Olsenella",
  "Neisseria",
  "Aerococcus",
  "Howardella",
  "Pseudoleptotrichia",
  "Succiniclasticum",
  "Treponema",
  "Lachnospira",
  "Lachnospiraceae_NK3A20_group",
  "Lachnospiraceae-NK3A20-group"
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
  file.path(
    OUT_DIR,
    "candidate_marker_genera_audit_progression127_7KB_GENUS_ONLY.tsv"
  )
)

############################################################
## Build LEfSe input
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
  meta_prog$LEfSe_class
)

names(
  class_row
) <- meta_prog$SampleID

lefse_input <- rbind(
  Class = class_row,
  lefse_feature_mat
)

LEFSE_INPUT_FILE <- file.path(
  OUT_DIR,
  "LEfSe_input_progression5_7KB_GENUS_ONLY.tsv"
)

LEFSE_FORMATTED_FILE <- file.path(
  OUT_DIR,
  "LEfSe_formatted_progression5_7KB_GENUS_ONLY.in"
)

LEFSE_RESULT_FILE <- file.path(
  OUT_DIR,
  "LEfSe_results_progression5_7KB_GENUS_ONLY.res"
)

OUT_ALL_RESULTS <- file.path(
  OUT_DIR,
  "LEfSe_all_results_progression5_7KB_GENUS_ONLY.csv"
)

OUT_CLEAN_MARKERS <- file.path(
  OUT_DIR,
  "LEfSe_clean_markers_progression5_7KB_GENUS_ONLY.csv"
)

OUT_CANDIDATES <- file.path(
  OUT_DIR,
  "LEfSe_candidate_results_progression5_7KB_GENUS_ONLY.csv"
)

OUT_SUMMARY <- file.path(
  OUT_DIR,
  "summary_LEfSe_progression5_7KB_GENUS_ONLY.tsv"
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
## Check LEfSe input
############################################################

cat(
  "\nChecking LEfSe input:\n"
)

cat(
  "Input file:",
  LEFSE_INPUT_FILE,
  "\n"
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
## Run LEfSe
##
## EXACTLY same conda strategy as original script
############################################################

cat(
  "\nRunning LEfSe inside R with clean conda environment...\n"
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
  "which R &&",
  "R --version | head -n 1 &&",
  "which lefse_format_input.py &&",
  "which lefse_run.py &&",
  "lefse_format_input.py",
  shQuote(
    LEFSE_INPUT_FILE
  ),
  shQuote(
    LEFSE_FORMATTED_FILE
  ),
  "-c 1 -o 1000000 &&",
  "lefse_run.py",
  shQuote(
    LEFSE_FORMATTED_FILE
  ),
  shQuote(
    LEFSE_RESULT_FILE
  ),
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
  "\nFormatted file exists:",
  file.exists(
    LEFSE_FORMATTED_FILE
  ),
  "\n"
)

cat(
  "Result file exists:",
  file.exists(
    LEFSE_RESULT_FILE
  ),
  "\n"
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
## Read and clean LEfSe result
##
## Same logic as original
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
    
    marker_status =
      if_else(
        !is.na(
          Class
        ) &
          !is.na(
            LDA
          ),
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
## Candidate result subset
############################################################

candidate_pattern <- paste(
  c(
    "UCG_005",
    "UCG-005",
    "Desulfovibrio",
    "Candidatus_Soleaferrea",
    "Mediterraneibacter",
    "Dorea",
    "Monoglobus",
    "Lactococcus",
    "Atopostipes",
    "Peptoclostridium",
    "Eubacterium__eligens",
    "Olsenella",
    "Neisseria",
    "Aerococcus",
    "Howardella",
    "Pseudoleptotrichia",
    "Succiniclasticum",
    "Treponema",
    "Lachnospira",
    "Lachnospiraceae_NK3A20"
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

dorea_desulfovibrio_status <- lefse_all %>%
  filter(
    Genus %in% c(
      "Dorea",
      "Desulfovibrio"
    )
  ) %>%
  select(
    Genus,
    log_max_mean,
    Class,
    LDA,
    p_value,
    marker_status
  )

############################################################
## Write results
############################################################

write_csv(
  lefse_all,
  OUT_ALL_RESULTS
)

write_csv(
  lefse_markers,
  OUT_CLEAN_MARKERS
)

write_csv(
  lefse_candidates,
  OUT_CANDIDATES
)

############################################################
## Summary
############################################################

summary_tbl <- tibble(
  Item = c(
    "Analysis",
    "Input count table",
    "Input taxonomy table",
    "Input metadata table",
    "Sample set",
    "Comparison",
    "Samples requested",
    "Group counts",
    "ASVs in progression127 subset",
    "ASVs with valid strict genus",
    "Genus-only features before LEfSe filtering",
    "Genus-only features kept for LEfSe",
    "Minimum prevalence n",
    "Minimum mean relative abundance",
    "LDA cutoff",
    "LEfSe input file",
    "LEfSe result file",
    "LEfSe all results CSV",
    "LEfSe clean markers CSV",
    "LEfSe candidate CSV",
    "Output directory"
  ),
  
  Value = c(
    "7KB progression127 5-group strict genus-only LEfSe",
    COUNT_FILE,
    TAX_FILE,
    META_FILE,
    "fixed 7KB progression127",
    "Polyp / UC_remission / UC_active / Dysplasia / CA",
    as.character(
      nrow(
        meta_prog
      )
    ),
    paste(
      group_counts$Progression5,
      group_counts$n,
      sep = "=",
      collapse = "; "
    ),
    as.character(
      nrow(
        count_prog
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
    LEFSE_INPUT_FILE,
    LEFSE_RESULT_FILE,
    OUT_ALL_RESULTS,
    OUT_CLEAN_MARKERS,
    OUT_CANDIDATES,
    OUT_DIR
  )
)

write_tsv(
  summary_tbl,
  OUT_SUMMARY
)

############################################################
## Print key results
############################################################

cat(
  "\n============================================================\n"
)

cat(
  "7KB progression127 5-group genus-only LEfSe finished\n"
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
  "\nTaxonomy audit:\n"
)

print(
  taxonomy_audit,
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
  "\nCandidate marker genus filtering audit:\n"
)

print(
  candidate_audit,
  n = Inf,
  width = Inf
)

cat(
  "\nALL LEfSe markers:\n"
)

print(
  lefse_markers,
  n = Inf,
  width = Inf
)

cat(
  "\nCandidate genera results:\n"
)

print(
  lefse_candidates,
  n = Inf,
  width = Inf
)

cat(
  "\nDorea / Desulfovibrio status:\n"
)

print(
  dorea_desulfovibrio_status,
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
  "All results CSV:   ",
  OUT_ALL_RESULTS,
  "\n"
)

cat(
  "Clean markers CSV: ",
  OUT_CLEAN_MARKERS,
  "\n"
)

cat(
  "Candidate CSV:     ",
  OUT_CANDIDATES,
  "\n"
)

cat(
  "Summary:           ",
  OUT_SUMMARY,
  "\n"
)

cat(
  "Output directory:  ",
  OUT_DIR,
  "\n"
)

cat(
  "\nDone.\n"
)