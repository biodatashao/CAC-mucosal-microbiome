


############################################################
## 03_01_genus_composition_5groups_source.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
##
## Genus-level 5-group composition plot
## 7KB progression127 dataset
##
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
  library(ggplot2)
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
  "genus_composition_5groups_top20_remove_others_renormalized"
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

############################################################
## Parameters
############################################################

GROUP_ORDER <- c(
  "Polyp",
  "UC_remission",
  "UC_active",
  "Dysplasia",
  "CA"
)

TOP_N_GENERA <- 20

############################################################
## Output files
############################################################

OUT_GROUP_MEAN_ALL <- file.path(
  OUT_DIR,
  "genus_group_mean_all_genera_before_top20_7KB_progression127.csv"
)

OUT_GROUP_MEAN_TOP_WITH_OTHERS <- file.path(
  OUT_DIR,
  "genus_group_mean_top20_with_Others_7KB_progression127.csv"
)

OUT_GROUP_MEAN_TOP_NO_OTHERS_RENORM <- file.path(
  OUT_DIR,
  "genus_group_mean_top20_Others_removed_renormalized100_7KB_progression127.csv"
)

OUT_SAMPLE_RA <- file.path(
  OUT_DIR,
  "sample_level_genus_relative_abundance_7KB_progression127.csv"
)

OUT_TOP20 <- file.path(
  OUT_DIR,
  "top20_genera_overall_mean_abundance_7KB_progression127.csv"
)

OUT_PDF <- file.path(
  OUT_DIR,
  "Fig_genus_composition_5groups_top20_Others_removed_renormalized100_7KB_progression127.pdf"
)

OUT_PNG <- file.path(
  OUT_DIR,
  "Fig_genus_composition_5groups_top20_Others_removed_renormalized100_7KB_progression127.png"
)

OUT_TIFF <- file.path(
  OUT_DIR,
  "Fig_genus_composition_5groups_top20_Others_removed_renormalized100_7KB_progression127.tiff"
)

############################################################
## Helper functions
############################################################

stop_if_missing <- function(files) {
  missing_files <- files[!file.exists(files)]
  
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

############################################################
## IMPORTANT:
## The following taxonomy parsing logic is taken from the
## original progression LEfSe script.
############################################################

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
    str_detect(
      x0,
      "uncultured"
    ) |
    str_detect(
      x0,
      "unclassified"
    ) |
    str_detect(
      x0,
      "unknown"
    ) |
    str_detect(
      x0,
      "metagenome"
    ) |
    str_detect(
      x0,
      "norank"
    ) |
    str_detect(
      x0,
      "no_rank"
    ) |
    str_detect(
      x0,
      "ambiguous"
    ) |
    str_detect(
      x0,
      "^[dkpcofs]__"
    ) |
    str_detect(
      x0,
      "k__|p__|c__|o__|f__|s__"
    )
}

############################################################
## Read data
############################################################

stop_if_missing(
  c(
    COUNT_FILE,
    TAX_FILE,
    META_FILE
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

cat(
  "\nCount table columns:\n"
)

print(
  colnames(count_df)
)

cat(
  "\nTaxonomy table columns:\n"
)

print(
  colnames(tax_df)
)

cat(
  "\nMetadata table columns:\n"
)

print(
  colnames(meta_df)
)

############################################################
## Detect columns
############################################################

asv_col_count <- colnames(
  count_df
)[1]

asv_col_tax <- colnames(
  tax_df
)[1]

taxonomy_col <- detect_taxonomy_col(
  tax_df
)

cat(
  "\nTaxonomy column used: ",
  taxonomy_col,
  "\n",
  sep = ""
)

if (!"SampleID" %in% colnames(meta_df)) {
  stop(
    "metadata table does not contain SampleID column."
  )
}

if (!"Progression5" %in% colnames(meta_df)) {
  stop(
    "metadata table does not contain Progression5 column."
  )
}

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
## Taxonomy -> genus
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
    
    Genus_for_plot =
      extract_genus_strict(
        Taxonomy_source_for_genus
      ),
    
    Genus_valid_for_plot =
      !is_bad_genus(
        Genus_for_plot
      )
  )

cat(
  "\nTaxonomy audit:\n"
)

cat(
  "Total taxonomy ASVs: ",
  nrow(tax_df),
  "\n",
  sep = ""
)

cat(
  "ASVs with valid genus: ",
  sum(
    tax_df$Genus_valid_for_plot,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "ASVs without valid genus: ",
  sum(
    !tax_df$Genus_valid_for_plot |
      is.na(
        tax_df$Genus_valid_for_plot
      )
  ),
  "\n",
  sep = ""
)

############################################################
## Fixed progression127 metadata
############################################################

meta_127 <- meta_df %>%
  filter(
    Progression5 %in% GROUP_ORDER
  ) %>%
  mutate(
    Progression5 = factor(
      Progression5,
      levels = GROUP_ORDER
    )
  )

group_counts <- meta_127 %>%
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

if (nrow(meta_127) != 127) {
  stop(
    "progression127 metadata contains ",
    nrow(meta_127),
    " samples, expected 127."
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
    paste0(
      "Progression5 group counts do not match expected ",
      "26/36/25/17/23."
    )
  )
}

############################################################
## Align sample order
############################################################

sample_keep <- meta_127$SampleID

missing_samples <- setdiff(
  sample_keep,
  colnames(count_mat)
)

if (length(missing_samples) > 0) {
  
  stop(
    "These metadata samples are missing from count table:\n",
    paste(
      missing_samples,
      collapse = ", "
    )
  )
}

count_mat <- count_mat[
  ,
  sample_keep,
  drop = FALSE
]

############################################################
## Align ASVs
############################################################

common_asvs <- intersect(
  rownames(count_mat),
  tax_df$ASV
)

cat(
  "\nASV alignment:\n"
)

cat(
  "ASVs in count table: ",
  nrow(count_mat),
  "\n",
  sep = ""
)

cat(
  "ASVs in taxonomy table: ",
  nrow(tax_df),
  "\n",
  sep = ""
)

cat(
  "Common ASVs: ",
  length(common_asvs),
  "\n",
  sep = ""
)

count_mat <- count_mat[
  common_asvs,
  ,
  drop = FALSE
]

tax_df <- tax_df %>%
  filter(
    ASV %in% common_asvs
  ) %>%
  arrange(
    match(
      ASV,
      rownames(count_mat)
    )
  )

stopifnot(
  identical(
    tax_df$ASV,
    rownames(count_mat)
  )
)

stopifnot(
  identical(
    colnames(count_mat),
    meta_127$SampleID
  )
)

############################################################
## Collapse ASVs to genus
############################################################

tax_genus <- tax_df %>%
  filter(
    Genus_valid_for_plot
  ) %>%
  select(
    ASV,
    Genus_for_plot
  )

count_genus_input <- count_mat[
  tax_genus$ASV,
  ,
  drop = FALSE
]

genus_count <- rowsum(
  count_genus_input,
  group = tax_genus$Genus_for_plot,
  reorder = FALSE
)

genus_count <- genus_count[
  rowSums(genus_count) > 0,
  ,
  drop = FALSE
]

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
  "Genus features after collapsing: ",
  nrow(genus_count),
  "\n",
  sep = ""
)

############################################################
## Relative abundance
############################################################

sample_depth_genus <- colSums(
  genus_count
)

if (any(sample_depth_genus <= 0)) {
  
  stop(
    "At least one sample has zero genus-level reads after genus filtering."
  )
}

genus_ra <- sweep(
  genus_count,
  2,
  sample_depth_genus,
  "/"
)

genus_ra[
  is.na(genus_ra)
] <- 0

sample_genus_ra_long <- genus_ra %>%
  as.data.frame(
    check.names = FALSE
  ) %>%
  rownames_to_column(
    "Genus"
  ) %>%
  pivot_longer(
    cols = -Genus,
    names_to = "SampleID",
    values_to = "Relative_abundance"
  ) %>%
  left_join(
    meta_127 %>%
      select(
        SampleID,
        Progression5
      ),
    by = "SampleID"
  ) %>%
  mutate(
    Progression5 = factor(
      Progression5,
      levels = GROUP_ORDER
    ),
    
    Relative_abundance_percent =
      Relative_abundance * 100
  )

write_csv(
  sample_genus_ra_long,
  OUT_SAMPLE_RA
)

############################################################
## Group mean abundance across all genera
############################################################

group_mean_all <- sample_genus_ra_long %>%
  group_by(
    Progression5,
    Genus
  ) %>%
  summarise(
    mean_RA = mean(
      Relative_abundance
    ),
    
    mean_percent = mean(
      Relative_abundance_percent
    ),
    
    .groups = "drop"
  )

write_csv(
  group_mean_all,
  OUT_GROUP_MEAN_ALL
)

############################################################
## Top 20 genera by overall mean abundance
############################################################

top_genera_table <- sample_genus_ra_long %>%
  group_by(
    Genus
  ) %>%
  summarise(
    overall_mean_RA = mean(
      Relative_abundance
    ),
    
    overall_mean_percent = mean(
      Relative_abundance_percent
    ),
    
    prevalence_n = sum(
      Relative_abundance > 0
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(
      overall_mean_RA
    )
  ) %>%
  slice_head(
    n = TOP_N_GENERA
  )

top_genera <- top_genera_table$Genus

write_csv(
  top_genera_table,
  OUT_TOP20
)

cat(
  "\nTop 20 genera used for visualization:\n"
)

print(
  top_genera_table,
  n = Inf,
  width = Inf
)

############################################################
## Top20 + Others
############################################################

sample_top_with_others <- sample_genus_ra_long %>%
  mutate(
    Genus_plot = if_else(
      Genus %in% top_genera,
      Genus,
      "Others"
    )
  ) %>%
  group_by(
    SampleID,
    Progression5,
    Genus_plot
  ) %>%
  summarise(
    Relative_abundance = sum(
      Relative_abundance
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    Relative_abundance_percent =
      Relative_abundance * 100
  )

group_top_with_others <- sample_top_with_others %>%
  group_by(
    Progression5,
    Genus_plot
  ) %>%
  summarise(
    mean_RA = mean(
      Relative_abundance
    ),
    
    mean_percent = mean(
      Relative_abundance_percent
    ),
    
    .groups = "drop"
  )

write_csv(
  group_top_with_others,
  OUT_GROUP_MEAN_TOP_WITH_OTHERS
)

############################################################
## Remove Others and renormalize to 100%
############################################################

group_top_no_others_renorm <- group_top_with_others %>%
  filter(
    Genus_plot != "Others"
  ) %>%
  group_by(
    Progression5
  ) %>%
  mutate(
    sum_without_others_percent = sum(
      mean_percent
    ),
    
    renormalized_percent =
      100 *
      mean_percent /
      sum_without_others_percent
  ) %>%
  ungroup()

write_csv(
  group_top_no_others_renorm,
  OUT_GROUP_MEAN_TOP_NO_OTHERS_RENORM
)

############################################################
## Order genera exactly as original composition script:
## by overall displayed abundance
############################################################

genus_order <- group_top_no_others_renorm %>%
  group_by(
    Genus_plot
  ) %>%
  summarise(
    overall_plot_percent = mean(
      renormalized_percent
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    overall_plot_percent
  ) %>%
  pull(
    Genus_plot
  )

plot_df <- group_top_no_others_renorm %>%
  mutate(
    Progression5 = factor(
      Progression5,
      levels = GROUP_ORDER
    ),
    
    Genus_plot = factor(
      Genus_plot,
      levels = genus_order
    )
  )

############################################################
## Plot
############################################################

p <- ggplot(
  plot_df,
  aes(
    x = Progression5,
    y = renormalized_percent,
    fill = Genus_plot
  )
) +
  geom_col(
    width = 0.78,
    color = "white",
    linewidth = 0.12
  ) +
  scale_y_continuous(
    expand = c(
      0,
      0
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      100
    )
  ) +
  labs(
    x = NULL,
    
    y = "Relative abundance among displayed genera (%)",
    
    fill = "Genus",
    
    title = paste0(
      "Genus-level composition across progression groups\n",
      "7KB progression127 cohort; top 20 genera; ",
      "Others removed and re-normalized"
    )
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.title = element_text(
      size = 9
    ),
    
    plot.title = element_text(
      hjust = 0.5,
      size = 12
    )
  )

print(p)

############################################################
## Save figure
############################################################

ggsave(
  filename = OUT_PDF,
  plot = p,
  width = 9.5,
  height = 6.5,
  units = "in"
)

ggsave(
  filename = OUT_PNG,
  plot = p,
  width = 9.5,
  height = 6.5,
  units = "in",
  dpi = 300
)

ggsave(
  filename = OUT_TIFF,
  plot = p,
  width = 9.5,
  height = 6.5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

############################################################
## Audit output
############################################################

others_percent <- group_top_with_others %>%
  filter(
    Genus_plot == "Others"
  ) %>%
  select(
    Progression5,
    Others_mean_percent = mean_percent
  )

cat(
  "\n============================================================\n"
)

cat(
  "Genus composition finished\n"
)

cat(
  "============================================================\n"
)

cat(
  "\nGroup counts:\n"
)

print(
  group_counts,
  n = Inf,
  width = Inf
)

cat(
  "\nTop 20 genera:\n"
)

print(
  top_genera_table,
  n = Inf,
  width = Inf
)

cat(
  "\nOthers percentage before removal:\n"
)

print(
  others_percent,
  n = Inf,
  width = Inf
)

cat(
  "\nOutput directory:\n"
)

cat(
  OUT_DIR,
  "\n"
)

cat(
  "\nDone.\n"
)