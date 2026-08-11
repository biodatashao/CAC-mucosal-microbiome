#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 03_04_plot_Figure3B_genus_composition_CAC_vs_nonCAC.R
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

color_mapping_file <- file.path(
  rprojroot::find_root(rprojroot::has_file("config.R")),
  "scripts", "03_taxa_LEfSe_Figure3", "genus_color_mapping.R"
)

if (!file.exists(color_mapping_file)) {
  stop(
    "Shared genus color mapping file not found: ",
    color_mapping_file
  )
}

source(color_mapping_file)

if (!exists("GENUS_COLORS")) {
  stop(
    "GENUS_COLORS was not created after sourcing: ",
    color_mapping_file
  )
}

input_file <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "genus_composition_CA23_vs_nonCA23_top20_remove_contaminants_remove_others_renormalized",
  "genus_group_mean_top20_Others_removed_renormalized100_7KB_CA23_nonCA23.csv"
)

out_dir <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "Figure3_final_panels",
  "Figure3B_genus_composition_CA23_vs_nonCA23"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file
  )
}

message(
  "Using input file: ",
  input_file
)

############################################################
## 2. Read data
############################################################

dat0 <- read_csv(
  input_file,
  show_col_types = FALSE
)

message(
  "Input columns: ",
  paste(
    names(dat0),
    collapse = ", "
  )
)

############################################################
## 3. Standardize columns
##
## Same logic as original Figure2B script.
############################################################

group_col <- intersect(
  c(
    "Group_pair",
    "Group",
    "group",
    "ComparisonGroup",
    "comparison_group"
  ),
  names(dat0)
)[1]

genus_col <- intersect(
  c(
    "Genus_plot",
    "Genus",
    "genus",
    "Taxon",
    "taxon"
  ),
  names(dat0)
)[1]

abund_col <- intersect(
  c(
    "renormalized_percent",
    "Renormalized_percent",
    "mean_percent",
    "Mean_percent",
    "mean_RA",
    "Mean_RA",
    "mean_ra",
    "Mean_RA_percent"
  ),
  names(dat0)
)[1]

if (
  is.na(group_col) ||
  is.na(genus_col) ||
  is.na(abund_col)
) {
  stop(
    "Cannot identify required columns.\n",
    "Detected group_col: ", group_col, "\n",
    "Detected genus_col: ", genus_col, "\n",
    "Detected abund_col: ", abund_col, "\n",
    "Available columns are: ",
    paste(
      names(dat0),
      collapse = ", "
    )
  )
}

############################################################
## 4. Prepare plotting data
############################################################

plot_dat <- dat0 %>%
  transmute(
    Group_raw = as.character(
      .data[[group_col]]
    ),
    Genus = as.character(
      .data[[genus_col]]
    ),
    Display_abundance = as.numeric(
      .data[[abund_col]]
    )
  ) %>%
  mutate(
    Group = recode(
      Group_raw,
      "nonCA" = "nonCAC",
      "NonCA" = "nonCAC",
      "nonCAC" = "nonCAC",
      "NonCAC" = "nonCAC",
      "CA" = "CAC",
      "CAC" = "CAC",
      .default = Group_raw
    ),
    
    Group = factor(
      Group,
      levels = c(
        "nonCAC",
        "CAC"
      )
    ),
    
    ########################################################
    ## Harmonize current 7KB SILVA naming with the original
    ## manuscript color mapping nomenclature.
    ##
    ## Naming aliases only; biological taxa are unchanged.
    ########################################################
    
    Genus = recode(
      Genus,
      "UCG-005" = "UCG_005",
      "UCG-009" = "UCG_009",
      "Escherichia-Shigella" = "Escherichia_Shigella",
      "Christensenellaceae_R-7_group" =
        "Christensenellaceae_R_7_group",
      .default = Genus
    ),
    
    ########################################################
    ## Original display convention
    ########################################################
    
    Genus = str_replace_all(
      Genus,
      "_",
      " "
    ),
    
    Genus = str_squish(
      Genus
    )
  ) %>%
  filter(
    !is.na(Group),
    !is.na(Genus),
    !is.na(Display_abundance)
  ) %>%
  filter(
    !str_detect(
      str_to_lower(Genus),
      "^others?$|^other$"
    )
  )

############################################################
## 5. Convert 0-1 abundance to percent if needed
##
## Original Figure2B logic.
############################################################

if (
  max(
    plot_dat$Display_abundance,
    na.rm = TRUE
  ) <= 1.5
) {
  plot_dat <- plot_dat %>%
    mutate(
      Display_abundance =
        Display_abundance * 100
    )
}

############################################################
## 6. Re-normalize each group to exactly 100%
##
## Original Figure2B logic.
############################################################

plot_dat <- plot_dat %>%
  group_by(
    Group
  ) %>%
  mutate(
    Display_abundance =
      Display_abundance /
      sum(
        Display_abundance,
        na.rm = TRUE
      ) *
      100
  ) %>%
  ungroup()

############################################################
## 7. Order genera by overall mean abundance
##
## Original Figure2B logic:
## legend follows stack from top to bottom.
############################################################

genus_order <- plot_dat %>%
  group_by(
    Genus
  ) %>%
  summarise(
    mean_abundance = mean(
      Display_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    mean_abundance
  ) %>%
  pull(
    Genus
  )

plot_dat <- plot_dat %>%
  mutate(
    Genus = factor(
      Genus,
      levels = genus_order
    )
  )

############################################################
## 8. Shared genus palette
##
## Reuse original manuscript palette unchanged.
############################################################

genus_palette_all <- GENUS_COLORS

names(genus_palette_all) <- str_replace_all(
  names(genus_palette_all),
  "_",
  " "
)

names(genus_palette_all) <- str_squish(
  names(genus_palette_all)
)

displayed_genera <- levels(
  plot_dat$Genus
)

missing_colors <- setdiff(
  displayed_genera,
  names(genus_palette_all)
)

if (length(missing_colors) > 0) {
  stop(
    "Missing shared colors for the following genera:\n",
    paste(
      missing_colors,
      collapse = ", "
    )
  )
}

genus_palette <- genus_palette_all[
  displayed_genera
]

font_family <- "Helvetica"

############################################################
## 9. Audit
############################################################

message(
  "Number of displayed genera: ",
  length(displayed_genera)
)

message(
  "Displayed genera:"
)

print(
  tibble::tibble(
    Order = seq_along(
      displayed_genera
    ),
    Genus = displayed_genera
  ),
  n = Inf
)

group_sum_check <- plot_dat %>%
  group_by(
    Group
  ) %>%
  summarise(
    Total_percent = sum(
      Display_abundance
    ),
    .groups = "drop"
  )

message(
  "Group abundance sums after final re-normalization:"
)

print(
  group_sum_check,
  n = Inf
)

############################################################
## 10. Plot
##
## Exact settings from original Figure2B.
############################################################

p <- ggplot(
  plot_dat,
  aes(
    x = Group,
    y = Display_abundance,
    fill = Genus
  )
) +
  geom_col(
    width = 0.58,
    color = "white",
    linewidth = 0.22
  ) +
  scale_fill_manual(
    values = genus_palette,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      100,
      25
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      100
    ),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Relative abundance of top 20 genera (%)",
    fill = "Genus"
  ) +
  theme_classic(
    base_size = 8,
    base_family = font_family
  ) +
  theme(
    text = element_text(
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(
        r = 6
      )
    ),
    
    axis.text.x = element_text(
      size = 9,
      color = "black",
      angle = 35,
      hjust = 1,
      vjust = 1
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
    
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 7.4
    ),
    
    legend.key.size = unit(
      3.3,
      "mm"
    ),
    
    legend.spacing.y = unit(
      0.8,
      "mm"
    ),
    
    plot.margin = margin(
      6,
      8,
      6,
      6
    )
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      reverse = FALSE,
      keyheight = unit(
        3.5,
        "mm"
      ),
      keywidth = unit(
        4.5,
        "mm"
      )
    )
  )

print(p)

############################################################
## 11. Export
##
## Same dimensions as original Figure2B:
## 6.8 x 4.4 inches, 600 dpi
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
    family = font_family,
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
        "Package 'ragg' is required on non-macOS systems."
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
## 12. Export final panel
############################################################

output_stem <- paste0(
  "Figure3B_genus_composition_CA23_vs_nonCA23_",
  "top20_3contaminants_removed_others_removed_",
  "renormalized_7KB"
)

export_plot(
  p,
  output_stem,
  width = 6.8,
  height = 4.4
)

############################################################
## 13. Save plotting data
############################################################

write_csv(
  plot_dat,
  file.path(
    out_dir,
    "Figure3B_genus_composition_CA23_vs_nonCA23_plot_data_7KB.csv"
  )
)

############################################################
## 14. Save final ggplot object
############################################################

figure3B_rds_file <- file.path(
  out_dir,
  "Figure3B_genus_composition_CA23_vs_nonCA23_shared_colors_7KB.rds"
)

saveRDS(
  p,
  figure3B_rds_file
)

############################################################
## 15. Finish
############################################################

message(
  "Done. Outputs written to: ",
  out_dir
)

message(
  "Saved Figure 3B plot object: ",
  figure3B_rds_file
)