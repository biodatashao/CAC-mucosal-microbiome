#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 03_06_plot_Figure3_LEfSe_progression5_markers.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
##
## 7KB progression127:
## A. LEfSe LDA score barplot
## B. Row-scaled group mean abundance heatmap
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})

############################################################
## 1. Paths
############################################################

base_dir <- file.path(PROJECT_ROOT, "output/analysis")

lefse_dir <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "LEfSe_progression5_7KB_GENUS_ONLY"
)

marker_file <- file.path(
  lefse_dir,
  "LEfSe_clean_markers_progression5_7KB_GENUS_ONLY.csv"
)

abundance_file <- file.path(
  lefse_dir,
  "genus_relative_abundance_progression127_LEfSe_features_GENUS_ONLY_7KB.tsv"
)

metadata_file <- file.path(
  lefse_dir,
  "metadata_progression5_7KB.tsv"
)

out_dir <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "Figure3_final_panels",
  "Figure3_LEfSe_progression5"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message(
  "Using marker file: ",
  marker_file
)

message(
  "Using abundance file: ",
  abundance_file
)

message(
  "Using metadata file: ",
  metadata_file
)

############################################################
## 2. Helpers
############################################################

pick_col <- function(dat, candidates, label) {
  
  hit <- intersect(
    candidates,
    names(dat)
  )
  
  if (length(hit) == 0) {
    stop(
      "Cannot identify ",
      label,
      " column.\n",
      "Available columns are: ",
      paste(
        names(dat),
        collapse = ", "
      )
    )
  }
  
  hit[1]
}

clean_genus <- function(x) {
  
  x <- as.character(x)
  
  x <- x %>%
    str_replace(
      "^k__.*\\|",
      ""
    ) %>%
    str_replace(
      "^.*g__",
      ""
    ) %>%
    str_replace(
      "^g__",
      ""
    ) %>%
    str_replace_all(
      "^f__Eubacterium__eligens_group$",
      "[Eubacterium] eligens group"
    ) %>%
    str_replace_all(
      "^f__Eubacterium_eligens_group$",
      "[Eubacterium] eligens group"
    ) %>%
    str_replace_all(
      "^f_Eubacterium_eligens_group$",
      "[Eubacterium] eligens group"
    ) %>%
    str_replace_all(
      "^f Eubacterium eligens group$",
      "[Eubacterium] eligens group"
    ) %>%
    str_replace_all(
      "^Eubacterium eligens group$",
      "[Eubacterium] eligens group"
    ) %>%
    str_replace_all(
      "_",
      " "
    ) %>%
    str_replace_all(
      "^f  Eubacterium  eligens group$",
      "[Eubacterium] eligens group"
    ) %>%
    str_replace_all(
      "^f Eubacterium eligens group$",
      "[Eubacterium] eligens group"
    ) %>%
    str_replace_all(
      "Oscillospiraceae UCG 005",
      "Oscillospiraceae UCG-005"
    ) %>%
    str_replace_all(
      "UCG 005",
      "UCG-005"
    ) %>%
    str_squish()
  
  x
}

format_group <- function(x) {
  
  recode(
    as.character(x),
    "Polyp" = "Polyp",
    "UC_remission" = "UC remission",
    "UC remission" = "UC remission",
    "UC_active" = "UC active",
    "UC active" = "UC active",
    "Dysplasia" = "Dysplasia",
    "CA" = "CAC",
    "CAC" = "CAC",
    .default = as.character(x)
  )
}

export_plot <- function(
    plot,
    filename,
    width,
    height,
    dpi = 600
) {
  
  pdf_file <- file.path(
    out_dir,
    paste0(
      filename,
      ".pdf"
    )
  )
  
  png_file <- file.path(
    out_dir,
    paste0(
      filename,
      ".png"
    )
  )
  
  tif_file <- file.path(
    out_dir,
    paste0(
      filename,
      ".tiff"
    )
  )
  
  ##########################################################
  ## PDF
  ##########################################################
  
  grDevices::pdf(
    pdf_file,
    width = width,
    height = height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  
  print(plot)
  
  grDevices::dev.off()
  
  ##########################################################
  ## PNG / TIFF
  ##########################################################
  
  if (
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ) {
    
    ragg::agg_png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      background = "white"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    ragg::agg_tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      background = "white"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
  } else {
    
    if (
      Sys.info()[["sysname"]] != "Darwin"
    ) {
      stop(
        "Package 'ragg' is required for PNG/TIFF export on non-macOS systems."
      )
    }
    
    grDevices::png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white",
      type = "quartz"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    grDevices::tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white",
      compression = "lzw",
      type = "quartz"
    )
    
    print(plot)
    
    grDevices::dev.off()
  }
}

############################################################
## 3. Check input files
############################################################

required_files <- c(
  marker_file,
  abundance_file,
  metadata_file
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
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

############################################################
## 4. Read LEfSe markers
############################################################

markers0 <- read_csv(
  marker_file,
  show_col_types = FALSE
)

message(
  "Marker columns: ",
  paste(
    names(markers0),
    collapse = ", "
  )
)

genus_col <- pick_col(
  markers0,
  c(
    "Genus",
    "genus",
    "Genus_plot",
    "Feature",
    "feature",
    "Taxon",
    "taxon",
    "lda_feature",
    "name"
  ),
  "marker/genus"
)

group_col <- pick_col(
  markers0,
  c(
    "Enriched_group",
    "enriched_group",
    "Group",
    "group",
    "Class",
    "class",
    "lda_class",
    "max_group"
  ),
  "enriched group"
)

lda_col <- pick_col(
  markers0,
  c(
    "LDA",
    "lda",
    "LDA_score",
    "lda_score",
    "LDAscore",
    "lda.log10",
    "log10_LDA"
  ),
  "LDA score"
)

markers <- markers0 %>%
  transmute(
    Genus_raw = .data[[genus_col]],
    Enriched_raw = .data[[group_col]],
    LDA_score = as.numeric(
      .data[[lda_col]]
    )
  ) %>%
  mutate(
    Genus = clean_genus(
      Genus_raw
    ),
    
    Enriched_group = format_group(
      Enriched_raw
    ),
    
    Enriched_group = factor(
      Enriched_group,
      levels = c(
        "Polyp",
        "UC remission",
        "UC active",
        "Dysplasia",
        "CAC"
      )
    )
  ) %>%
  filter(
    !is.na(Genus),
    !is.na(Enriched_group),
    !is.na(LDA_score)
  )

############################################################
## Keep all clean markers
##
## Same as original manuscript script
############################################################

n_per_group <- Inf

if (
  is.finite(
    n_per_group
  )
) {
  
  markers <- markers %>%
    group_by(
      Enriched_group
    ) %>%
    slice_max(
      order_by = LDA_score,
      n = n_per_group,
      with_ties = FALSE
    ) %>%
    ungroup()
}

marker_order <- markers %>%
  arrange(
    Enriched_group,
    LDA_score
  ) %>%
  pull(
    Genus
  )

markers <- markers %>%
  mutate(
    Genus = factor(
      Genus,
      levels = unique(
        marker_order
      )
    )
  )

message(
  "Number of LEfSe markers plotted: ",
  nrow(markers)
)

message(
  "Marker counts by enriched group:"
)

print(
  markers %>%
    count(
      Enriched_group,
      name = "n_markers"
    ),
  n = Inf
)

############################################################
## 5. Read abundance and metadata
############################################################

abund0 <- read_tsv(
  abundance_file,
  show_col_types = FALSE
)

meta0 <- read_tsv(
  metadata_file,
  show_col_types = FALSE
)

message(
  "Abundance columns: ",
  paste(
    names(abund0),
    collapse = ", "
  )
)

message(
  "Metadata columns: ",
  paste(
    names(meta0),
    collapse = ", "
  )
)

sample_col_meta <- pick_col(
  meta0,
  c(
    "SampleID",
    "Sample_ID",
    "sample_id",
    "sample",
    "Sample"
  ),
  "metadata sample ID"
)

group_col_meta <- pick_col(
  meta0,
  c(
    "Progression5",
    "Group",
    "group",
    "Class",
    "class"
  ),
  "metadata group"
)

meta <- meta0 %>%
  transmute(
    SampleID = as.character(
      .data[[sample_col_meta]]
    ),
    
    Group = format_group(
      .data[[group_col_meta]]
    )
  ) %>%
  mutate(
    Group = factor(
      Group,
      levels = c(
        "Polyp",
        "UC remission",
        "UC active",
        "Dysplasia",
        "CAC"
      )
    )
  )

############################################################
## 6. Parse abundance table
##
## 04a output format is fixed:
## first column = Genus
## remaining columns = 127 SampleIDs
############################################################

if (!"Genus" %in% names(abund0)) {
  stop(
    "Expected a 'Genus' column in abundance table.\n",
    "Available columns are: ",
    paste(names(abund0), collapse = ", ")
  )
}

sample_cols_abund <- setdiff(
  names(abund0),
  "Genus"
)

if (length(sample_cols_abund) != 127) {
  stop(
    "Expected 127 sample columns in abundance table, found ",
    length(sample_cols_abund),
    "."
  )
}

############################################################
## Check abundance sample IDs against metadata
############################################################

missing_in_metadata <- setdiff(
  sample_cols_abund,
  meta$SampleID
)

missing_in_abundance <- setdiff(
  meta$SampleID,
  sample_cols_abund
)

if (length(missing_in_metadata) > 0) {
  stop(
    "Abundance table contains sample(s) absent from metadata:\n",
    paste(missing_in_metadata, collapse = ", ")
  )
}

if (length(missing_in_abundance) > 0) {
  stop(
    "Metadata contains sample(s) absent from abundance table:\n",
    paste(missing_in_abundance, collapse = ", ")
  )
}

############################################################
## Matrix -> long format
############################################################

abund_long <- abund0 %>%
  mutate(
    Genus = clean_genus(Genus)
  ) %>%
  pivot_longer(
    cols = all_of(sample_cols_abund),
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  mutate(
    SampleID = as.character(SampleID),
    Abundance = as.numeric(Abundance)
  ) %>%
  left_join(
    meta,
    by = "SampleID"
  )

############################################################
## Audit
############################################################

if (any(is.na(abund_long$Group))) {
  bad_samples <- abund_long %>%
    filter(is.na(Group)) %>%
    distinct(SampleID) %>%
    pull(SampleID)
  
  stop(
    "Group could not be assigned for sample(s):\n",
    paste(bad_samples, collapse = ", ")
  )
}

cat(
  "\nAbundance parsing audit:\n"
)

cat(
  "Genus features: ",
  nrow(abund0),
  "\n",
  sep = ""
)

cat(
  "Sample columns: ",
  length(sample_cols_abund),
  "\n",
  sep = ""
)

cat(
  "Long-format rows: ",
  nrow(abund_long),
  "\n",
  sep = ""
)

cat(
  "Unique samples after pivot: ",
  dplyr::n_distinct(abund_long$SampleID),
  "\n",
  sep = ""
)

############################################################
## 7. Convert abundance if input is percent
############################################################

if (
  max(
    abund_long$Abundance,
    na.rm = TRUE
  ) > 1.5
) {
  
  abund_long <- abund_long %>%
    mutate(
      Abundance =
        Abundance / 100
    )
}

############################################################
## 8. Heatmap data
##
## Original method:
## group mean RA
## -> log10(mean_RA + 1e-6)
## -> row-wise z scaling
############################################################

heat_dat <- abund_long %>%
  filter(
    !is.na(Group),
    
    Group %in% c(
      "Polyp",
      "UC remission",
      "UC active",
      "Dysplasia",
      "CAC"
    ),
    
    Genus %in% as.character(
      markers$Genus
    )
  ) %>%
  group_by(
    Genus,
    Group
  ) %>%
  summarise(
    mean_RA = mean(
      Abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  complete(
    Genus = as.character(
      markers$Genus
    ),
    
    Group = factor(
      c(
        "Polyp",
        "UC remission",
        "UC active",
        "Dysplasia",
        "CAC"
      ),
      levels = c(
        "Polyp",
        "UC remission",
        "UC active",
        "Dysplasia",
        "CAC"
      )
    ),
    
    fill = list(
      mean_RA = 0
    )
  ) %>%
  group_by(
    Genus
  ) %>%
  mutate(
    scaled_abundance = as.numeric(
      scale(
        log10(
          mean_RA + 1e-6
        )
      )
    )
  ) %>%
  ungroup() %>%
  mutate(
    scaled_abundance = if_else(
      is.na(
        scaled_abundance
      ),
      0,
      scaled_abundance
    ),
    
    Genus = factor(
      Genus,
      levels = rev(
        levels(
          markers$Genus
        )
      )
    ),
    
    Group = factor(
      Group,
      levels = c(
        "Polyp",
        "UC remission",
        "UC active",
        "Dysplasia",
        "CAC"
      )
    )
  )

############################################################
## 9. Plot style
##
## Exact original colors.
############################################################

group_colors <- c(
  "Polyp" = "#7DA9B7",
  "UC remission" = "#8DBA91",
  "UC active" = "#E3A04F",
  "Dysplasia" = "#C77DA8",
  "CAC" = "#A32635"
)

base_theme <- theme_classic(
  base_size = 8,
  base_family = "Helvetica"
) +
  theme(
    text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    plot.title = element_text(
      size = 9,
      face = "bold",
      hjust = 0,
      margin = margin(
        b = 5
      )
    ),
    
    legend.title = element_text(
      size = 8.5,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )

############################################################
## 10. Panel A: LDA barplot
############################################################

pA_with_legend <- ggplot(
  markers,
  aes(
    x = LDA_score,
    y = Genus,
    fill = Enriched_group
  )
) +
  geom_col(
    width = 0.72
  ) +
  scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  labs(
    title = "LEfSe-identified genus markers across progression groups",
    x = "LDA score",
    y = NULL,
    fill = "Enriched group"
  ) +
  base_theme +
  theme(
    legend.position = "bottom",
    
    legend.justification = "left",
    
    legend.box.just = "left",
    
    legend.title = element_text(
      size = 8.5,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.size = unit(
      3.8,
      "mm"
    ),
    
    legend.spacing.x = unit(
      2.5,
      "mm"
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      5,
      5,
      4,
      5
    )
  ) +
  guides(
    fill = guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "left"
    )
  )

pA_legend <- cowplot::get_legend(
  pA_with_legend
)

pA_core <- pA_with_legend +
  theme(
    legend.position = "none",
    panel.grid = element_blank()
  )

pA <- cowplot::plot_grid(
  pA_core,
  pA_legend,
  ncol = 1,
  rel_heights = c(
    1,
    0.13
  ),
  align = "v",
  axis = "lr"
)

############################################################
## 11. Panel B: Heatmap
############################################################

pB <- ggplot(
  heat_dat,
  aes(
    x = Group,
    y = Genus,
    fill = scaled_abundance
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.45
  ) +
  scale_fill_gradient2(
    low = "#4E6FAE",
    mid = "white",
    high = "#B9443E",
    midpoint = 0,
    na.value = "white",
    name = "Row-scaled\nmean abundance"
  ) +
  labs(
    title = "Mean abundance patterns of LEfSe-identified genus markers",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(
    base_size = 8,
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      color = "black"
    ),
    
    plot.title = element_text(
      size = 9,
      face = "bold",
      hjust = 0,
      margin = margin(
        b = 5
      )
    ),
    
    axis.text.x = element_text(
      size = 8,
      color = "black",
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    panel.grid = element_blank(),
    
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    legend.title = element_text(
      size = 8.5,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.margin = margin(
      5,
      5,
      4,
      5
    )
  )

############################################################
## 12. Combine
############################################################

fig3 <- cowplot::plot_grid(
  pA,
  pB,
  ncol = 2,
  labels = c(
    "A",
    "B"
  ),
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = "Helvetica",
  label_x = c(
    0.00,
    0.00
  ),
  label_y = c(
    1.00,
    1.00
  ),
  hjust = c(
    -0.3,
    -0.3
  ),
  vjust = c(
    1.2,
    1.2
  ),
  rel_widths = c(
    1.18,
    1.0
  ),
  align = "h",
  axis = "tb"
)

############################################################
## 13. Export
############################################################

export_plot(
  fig3,
  "Figure3_LEfSe_progression5_genus_markers_7KB",
  width = 10.5,
  height = 5.4
)

export_plot(
  pA,
  "Figure3A_LEfSe_progression5_LDA_barplot_7KB",
  width = 5.4,
  height = 5.2
)

export_plot(
  pB,
  "Figure3B_LEfSe_progression5_marker_heatmap_7KB",
  width = 5.0,
  height = 5.2
)

############################################################
## 14. Save plotting data
############################################################

write_csv(
  markers,
  file.path(
    out_dir,
    "Figure3A_LEfSe_progression5_markers_plot_data_7KB.csv"
  )
)

write_csv(
  heat_dat,
  file.path(
    out_dir,
    "Figure3B_LEfSe_progression5_heatmap_plot_data_7KB.csv"
  )
)

############################################################
## 15. Save plot objects
############################################################

saveRDS(
  pA,
  file.path(
    out_dir,
    "Figure3A_LEfSe_progression5_LDA_barplot_7KB.rds"
  )
)

saveRDS(
  pB,
  file.path(
    out_dir,
    "Figure3B_LEfSe_progression5_marker_heatmap_7KB.rds"
  )
)

############################################################
## 16. Finish
############################################################

message(
  "Done. Outputs written to: ",
  out_dir
)

message(
  "Saved plot object: ",
  file.path(
    out_dir,
    "Figure3A_LEfSe_progression5_LDA_barplot_7KB.rds"
  )
)

message(
  "Saved plot object: ",
  file.path(
    out_dir,
    "Figure3B_LEfSe_progression5_marker_heatmap_7KB.rds"
  )
)