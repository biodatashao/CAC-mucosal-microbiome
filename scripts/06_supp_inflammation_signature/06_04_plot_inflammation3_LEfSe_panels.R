#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 06_04_plot_inflammation3_LEfSe_panels.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## Supplementary Figure:
##
## A. LEfSe LDA score barplot
## B. Row-scaled group mean abundance heatmap
##
## Inflammation spectrum:
## Polyp
## UC remission
## UC active
##
## Clean rerun version
############################################################


suppressPackageStartupMessages({
  
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(cowplot)
  
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
  "progression127"
)


LEFSE_DIR <- file.path(
  RUN_DIR,
  "06_LEfSe_UC_inflammation3_7KB_GENUS_ONLY"
)


MARKER_FILE <- file.path(
  LEFSE_DIR,
  "LEfSe_clean_markers_inflammation3_7KB_GENUS_ONLY.csv"
)


ASV_FILE <- file.path(
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


OUT_DIR <- file.path(
  LEFSE_DIR,
  "SupplementaryFigure_inflammation3_LEfSe_clean"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################
## 2. Parameters
############################################################


GROUP_ORDER <- c(
  "Polyp",
  "UC remission",
  "UC active"
)



REMOVE_GENUS <- c(
  "Aquitalea",
  "Oryzomicrobium",
  "Ornithinimicrobium",
  "Limnohabitans",
  "Fuscovulum"
)



############################################################
## 3. Helper functions
############################################################


format_group <- function(x){
  
  x <- as.character(x)
  
  x <- recode(
    x,
    "UC_remission" = "UC remission",
    "UC_active" = "UC active",
    .default = x
  )
  
  x
  
}



clean_taxon_name <- function(x){
  
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
      "_",
      " "
    ) %>%
    
    str_replace(
      "^f ",
      ""
    ) %>%
    
    str_squish()
  
  
  
  ## exact special taxonomy harmonization
  
  x <- case_when(
    
    x %in% c(
      "[Clostridium] innocuum group",
      "Clostridium innocuum",
      "Clostridium innocuum group",
      "Clostridium innocuum group group",
      "Clostridium innocuum group group group"
    ) ~
      "Clostridium innocuum group",
    
    
    x %in% c(
      "[Eubacterium] eligens group",
      "Eubacterium eligens group",
      "Eubacterium eligens group group",
      "Eubacterium eligens group group group"
    ) ~
      "Eubacterium eligens group",
    
    
    x %in% c(
      "Prevotellaceae UCG 003",
      "Prevotellaceae UCG-003"
    ) ~
      "Prevotellaceae UCG-003",
    
    
    x %in% c(
      "UCG 009",
      "UCG-009"
    ) ~
      "UCG-009",
    
    
    TRUE ~ x
    
  )
  
  
  x
  
}



export_plot <- function(
    plot,
    filename,
    width,
    height,
    dpi = 600
){
  
  pdf_file <- file.path(
    OUT_DIR,
    paste0(
      filename,
      ".pdf"
    )
  )
  
  
  png_file <- file.path(
    OUT_DIR,
    paste0(
      filename,
      ".png"
    )
  )
  
  
  tif_file <- file.path(
    OUT_DIR,
    paste0(
      filename,
      ".tiff"
    )
  )
  
  
  grDevices::pdf(
    pdf_file,
    width = width,
    height = height,
    family = "Helvetica",
    useDingbats = FALSE
  )
  
  print(plot)
  
  grDevices::dev.off()
  
  
  
  if(
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ){
    
    ragg::agg_png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    
    ragg::agg_tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
  }
  
}
############################################################
## 4. Read LEfSe markers
############################################################


markers_raw <- read_csv(
  MARKER_FILE,
  show_col_types = FALSE
)


cat(
  "Marker columns:\n"
)

print(
  names(markers_raw)
)



############################################################
## LEfSe csv structure:
##
## Class   = genus name
## LDA     = enriched group
## p_value = LDA score
############################################################


markers <- markers_raw %>%
  
  transmute(
    
    Genus = clean_taxon_name(
      Class
    ),
    
    Enriched_group = format_group(
      LDA
    ),
    
    LDA_score = as.numeric(
      p_value
    )
    
  ) %>%
  
  filter(
    
    !is.na(Genus),
    
    !is.na(Enriched_group),
    
    !is.na(LDA_score)
    
  )



############################################################
## Remove only five environmental genera
############################################################


markers <- markers %>%
  
  filter(
    !Genus %in% REMOVE_GENUS
  )



markers$Enriched_group <- factor(
  markers$Enriched_group,
  levels = GROUP_ORDER
)



cat(
  "Final marker number:\n"
)

print(
  markers %>%
    count(
      Enriched_group
    )
)



cat(
  "Total markers:",
  nrow(markers),
  "\n"
)



############################################################
## 5. Read ASV / taxonomy / metadata
############################################################


asv <- read.table(
  ASV_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)



taxonomy <- read.table(
  TAX_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)



metadata <- read.table(
  META_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)



############################################################
## 6. Select inflammation spectrum samples
############################################################


metadata <- metadata[
  metadata$Progression5 %in%
    c(
      "Polyp",
      "UC_remission",
      "UC_active"
    ),
]



metadata$Group <- factor(
  format_group(
    metadata$Progression5
  ),
  levels = GROUP_ORDER
)



sample_ids <- rownames(metadata)



asv <- asv[
  ,
  sample_ids
]



############################################################
## 7. Extract strict genus
############################################################


tax_string <- taxonomy$Taxonomy



genus <- sapply(
  strsplit(
    tax_string,
    ";"
  ),
  function(x){
    
    g <- x[
      grepl(
        "^g__",
        x
      )
    ]
    
    if(
      length(g) == 0
    ){
      return(NA)
    }
    
    g[1]
    
  }
)



genus_clean <- clean_taxon_name(
  genus
)



############################################################
## Remove invalid taxonomy
############################################################


keep_genus <- !
  
  grepl(
    "uncultured|unclassified|norank|unknown|ambiguous",
    genus_clean,
    ignore.case = TRUE
  )



keep_genus <- keep_genus &
  !is.na(genus_clean) &
  genus_clean != ""



asv_genus <- asv[
  keep_genus,
]



genus_clean <- genus_clean[
  keep_genus
]



############################################################
## 8. Aggregate ASV to genus
############################################################


genus_table <- rowsum(
  as.matrix(asv_genus),
  group = genus_clean
)



############################################################
## 9. Relative abundance
############################################################


genus_RA <- sweep(
  genus_table,
  2,
  colSums(genus_table),
  FUN = "/"
)



############################################################
## 10. Final genus harmonization
############################################################


rownames(genus_RA) <- clean_taxon_name(
  rownames(genus_RA)
)



############################################################
## Check marker matching
############################################################


marker_genus <- unique(
  markers$Genus
)



missing_marker <- setdiff(
  marker_genus,
  rownames(genus_RA)
)



cat(
  "\nMarkers not found in abundance table:\n"
)

print(
  missing_marker
)



if(
  length(missing_marker) > 0
){
  
  cat(
    "\nMissing markers:\n"
  )
  
  print(
    missing_marker
  )
  
  
  cat(
    "\nAvailable similar genus names:\n"
  )
  
  print(
    grep(
      paste(
        missing_marker,
        collapse = "|"
      ),
      rownames(genus_RA),
      value = TRUE
    )
  )
  
  
  stop(
    "Some LEfSe markers cannot be matched to abundance table."
  )
  
}



############################################################
## 11. Prepare heatmap long data
############################################################


heat_long <- as.data.frame(
  genus_RA
) %>%
  
  mutate(
    Genus = rownames(.)
  ) %>%
  
  filter(
    Genus %in% marker_genus
  ) %>%
  
  pivot_longer(
    cols = -Genus,
    names_to = "SampleID",
    values_to = "Abundance"
  ) %>%
  
  left_join(
    data.frame(
      SampleID = rownames(metadata),
      Group = metadata$Group
    ),
    by = "SampleID"
  )



############################################################
## Continue Part 3
############################################################
############################################################
## 12. Heatmap data
############################################################


heat_dat <- heat_long %>%
  
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
    Genus = marker_genus,
    Group = factor(
      GROUP_ORDER,
      levels = GROUP_ORDER
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
    
    scaled_abundance = ifelse(
      is.na(
        scaled_abundance
      ),
      0,
      scaled_abundance
    )
    
  )



############################################################
## Fix group order
############################################################


heat_dat$Group <- factor(
  heat_dat$Group,
  levels = GROUP_ORDER
)



############################################################
## Genus order
############################################################


marker_order <- markers %>%
  
  arrange(
    Enriched_group,
    LDA_score
  ) %>%
  
  pull(
    Genus
  )


markers$Genus <- factor(
  markers$Genus,
  levels = marker_order
)


heat_dat$Genus <- factor(
  heat_dat$Genus,
  levels = rev(
    unique(
      marker_order
    )
  )
)



############################################################
## 13. LDA barplot
############################################################


group_colors <- c(
  
  "Polyp" = "#7DA9B7",
  
  "UC remission" = "#8DBA91",
  
  "UC active" = "#E3A04F"
  
)



theme_base <- theme_classic(
  base_size = 8,
  base_family = "Helvetica"
) +
  
  theme(
    
    text = element_text(
      color = "black"
    ),
    
    axis.text = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.title = element_text(
      size = 9,
      face = "bold"
    )
    
  )



pA <- ggplot(
  
  markers,
  
  aes(
    
    x = LDA_score,
    
    y = Genus,
    
    fill = Enriched_group
    
  )
) +
  
  geom_col(
    width = 0.7
  ) +
  
  scale_fill_manual(
    values = group_colors
  ) +
  
  labs(
    
    title =
      "LEfSe-identified genus markers across inflammation spectrum",
    
    x =
      "LDA score",
    
    y =
      NULL,
    
    fill =
      "Enriched group"
    
  ) +
  
  theme_base +
  
  theme(
    
    legend.position = "bottom",
    
    panel.grid = element_blank()
    
  )



############################################################
## 14. Heatmap
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
    
    name =
      "Row-scaled\nmean abundance"
    
  ) +
  
  labs(
    
    title =
      "Mean abundance patterns of inflammation-associated LEfSe markers",
    
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
    
    axis.text.x = element_text(
      
      angle = 35,
      
      hjust = 1,
      
      color = "black",
      
      size = 8
      
    ),
    
    axis.text.y = element_text(
      
      color = "black",
      
      size = 8
      
    ),
    
    panel.grid = element_blank()
    
  )



############################################################
## 15. Export
############################################################


export_plot(
  
  pA,
  
  "SupplementaryFig_inflammation3_LEfSe_LDA_barplot_7KB",
  
  width = 5.5,
  
  height = 5.5
  
)



export_plot(
  
  pB,
  
  "SupplementaryFig_inflammation3_LEfSe_heatmap_7KB",
  
  width = 5.0,
  
  height = 5.5
  
)



############################################################
## 16. Save data
############################################################


write_csv(
  
  markers,
  
  file.path(
    
    OUT_DIR,
    
    "SupplementaryFig_inflammation3_LEfSe_marker_plot_data_7KB.csv"
    
  )
  
)



write_csv(
  
  heat_dat,
  
  file.path(
    
    OUT_DIR,
    
    "SupplementaryFig_inflammation3_LEfSe_heatmap_data_7KB.csv"
    
  )
  
)



saveRDS(
  
  pA,
  
  file.path(
    
    OUT_DIR,
    
    "SupplementaryFig_inflammation3_LEfSe_LDA_barplot_7KB.rds"
    
  )
  
)



saveRDS(
  
  pB,
  
  file.path(
    
    OUT_DIR,
    
    "SupplementaryFig_inflammation3_LEfSe_heatmap_7KB.rds"
    
  )
  
)



############################################################
## Finish
############################################################


message(
  "Done. Output written to: ",
  OUT_DIR
)