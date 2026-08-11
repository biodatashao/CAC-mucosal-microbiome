#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 03_08_plot_Figure3C_LEfSe_CAC_vs_nonCAC_markers.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
})

############################################################
## 1. Paths
############################################################

base_dir <- file.path(PROJECT_ROOT, "output/analysis")

lefse_dir <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "LEfSe_CA23_vs_nonCA23_7KB_GENUS_ONLY"
)

marker_file <- file.path(
  lefse_dir,
  "LEfSe_clean_markers_CA23_vs_nonCA23_unpaired_GENUS_ONLY_7KB.csv"
)

out_dir <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "Figure3_final_panels",
  "Figure3C_LEfSe_CA_vs_nonCA"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(marker_file)) {
  stop(
    "Marker file not found: ",
    marker_file
  )
}

message(
  "Using marker file: ",
  marker_file
)

############################################################
## 2. Helper: clean genus display names
##
## Same logic as original Figure3C.
############################################################

clean_genus <- function(x) {
  
  x %>%
    as.character() %>%
    str_replace(
      "^g__",
      ""
    ) %>%
    str_replace_all(
      "^f__Eubacterium__brachy_group$",
      "[Eubacterium] brachy group"
    ) %>%
    str_replace_all(
      "^f__Eubacterium_brachy_group$",
      "[Eubacterium] brachy group"
    ) %>%
    str_replace_all(
      "^f Eubacterium brachy group$",
      "[Eubacterium] brachy group"
    ) %>%
    str_replace_all(
      "_",
      " "
    ) %>%
    str_replace_all(
      "UCG 005",
      "UCG-005"
    ) %>%
    str_squish()
}

############################################################
## 3. Helper: group display names
##
## Same logic as original Figure3C.
############################################################

format_group <- function(x) {
  
  recode(
    as.character(x),
    
    "nonCA" = "nonCAC",
    "NonCA" = "nonCAC",
    "nonCAC" = "nonCAC",
    
    "CA" = "CAC",
    "CAC" = "CAC",
    
    .default = as.character(x)
  )
}

############################################################
## 4. Export
##
## Same dimensions and export logic as original.
## On this macOS workflow, do not use Cairo fallback.
############################################################

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
## 5. Read markers
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

required_cols <- c(
  "Genus",
  "Class",
  "LDA",
  "p_value"
)

missing_cols <- setdiff(
  required_cols,
  names(markers0)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required marker column(s): ",
    paste(
      missing_cols,
      collapse = ", "
    )
  )
}

############################################################
## 6. Predefined exclusions
##
## EXACT original Figure3C exclusion list.
############################################################

contaminant_genera <- c(
  "Noviherbaspirillum",
  "Marmoricola",
  "Pseudonocardia"
)

############################################################
## 7. Prepare plotting data
############################################################

plot_dat <- markers0 %>%
  transmute(
    Genus_raw = as.character(
      Genus
    ),
    
    Genus = clean_genus(
      Genus
    ),
    
    Enriched_group = format_group(
      Class
    ),
    
    LDA_score = as.numeric(
      LDA
    ),
    
    p_value = as.numeric(
      p_value
    )
  ) %>%
  filter(
    !Genus_raw %in% contaminant_genera
  ) %>%
  mutate(
    Enriched_group = factor(
      Enriched_group,
      levels = c(
        "nonCAC",
        "CAC"
      )
    )
  ) %>%
  arrange(
    Enriched_group,
    LDA_score
  ) %>%
  mutate(
    Genus = factor(
      Genus,
      levels = unique(
        Genus
      )
    )
  )

############################################################
## 8. Audit
############################################################

excluded_found <- markers0 %>%
  filter(
    Genus %in% contaminant_genera
  ) %>%
  select(
    Genus,
    Class,
    LDA,
    p_value
  )

message(
  "Raw LEfSe marker count: ",
  nrow(markers0)
)

message(
  "Excluded predefined likely contaminant markers: ",
  nrow(excluded_found)
)

message(
  "Displayed marker count: ",
  nrow(plot_dat)
)

message(
  "Displayed markers:"
)

print(
  plot_dat %>%
    select(
      Genus,
      Enriched_group,
      LDA_score,
      p_value
    ),
  n = Inf,
  width = Inf
)

message(
  "Excluded markers:"
)

print(
  excluded_found,
  n = Inf,
  width = Inf
)

############################################################
## Expected 7KB result audit
##
## Current 05a result:
## 7 raw markers
## - 6 nonCA
## - 1 CA
##
## After original predefined exclusion:
## - Noviherbaspirillum removed
## - Marmoricola removed
## - Pseudonocardia removed
##
## Expected displayed markers = 4
############################################################

if (nrow(markers0) != 7) {
  warning(
    "Raw LEfSe marker count is ",
    nrow(markers0),
    ", whereas the current audited 7KB result had 7 markers."
  )
}

if (nrow(plot_dat) != 4) {
  warning(
    "Displayed marker count is ",
    nrow(plot_dat),
    ", whereas applying the original Figure3C exclusion rule ",
    "to the current 7KB result should leave 4 markers."
  )
}

############################################################
## 9. Original colors
############################################################

group_colors <- c(
  "nonCAC" = "#A99E79",
  "CAC" = "#A32635"
)

############################################################
## 10. Plot
##
## EXACT styling from original Figure3C.
############################################################

p <- ggplot(
  plot_dat,
  aes(
    x = LDA_score,
    y = Genus,
    fill = Enriched_group
  )
) +
  geom_col(
    width = 0.62
  ) +
  scale_fill_manual(
    values = group_colors,
    drop = FALSE
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0,
        0.04
      )
    )
  ) +
  labs(
    title = "LEfSe-identified genus markers in CAC versus nonCAC",
    x = "LDA score",
    y = NULL,
    fill = "Enriched group"
  ) +
  theme_classic(
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
    
    axis.title.x = element_text(
      size = 9,
      face = "bold",
      margin = margin(
        t = 5
      )
    ),
    
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.text.y = element_text(
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
    
    panel.grid = element_blank(),
    
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
    
    plot.margin = margin(
      5,
      5,
      5,
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

print(p)

############################################################
## 11. Export
##
## Original Figure3C dimensions:
## 4.8 x 2.7 inch
############################################################

export_plot(
  p,
  "Figure3C_LEfSe_CA_vs_nonCA_genus_markers_filtered_7KB",
  width = 4.8,
  height = 2.7
)

############################################################
## 12. Save plot data
############################################################

write_csv(
  plot_dat,
  file.path(
    out_dir,
    "Figure3C_LEfSe_CA_vs_nonCA_genus_markers_filtered_plot_data_7KB.csv"
  )
)

write_csv(
  excluded_found,
  file.path(
    out_dir,
    "Figure3C_LEfSe_CA_vs_nonCA_excluded_contaminant_markers_7KB.csv"
  )
)

############################################################
## 13. Save plot object
############################################################

rds_file <- file.path(
  out_dir,
  "Figure3C_LEfSe_CA_vs_nonCA_genus_markers_filtered_7KB.rds"
)

saveRDS(
  p,
  rds_file
)

############################################################
## 14. Finish
############################################################

message(
  "Done. Outputs written to: ",
  out_dir
)

message(
  "Saved plot object: ",
  rds_file
)

message(
  "Displayed markers: ",
  paste(
    as.character(
      plot_dat$Genus
    ),
    collapse = ", "
  )
)

message(
  "Excluded likely contaminant genera: ",
  paste(
    contaminant_genera,
    collapse = ", "
  )
)