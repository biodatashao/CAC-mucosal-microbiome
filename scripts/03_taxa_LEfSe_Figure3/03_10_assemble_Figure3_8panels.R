#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 03_10_assemble_Figure3_8panels.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
##
## Purpose:
## Assemble the complete 7KB rerun Figure 3 using the same final layout
## and formatting as the original manuscript assembly script:
##
## Final panel order:
## A. Genus composition across five progression groups
## B. Genus composition in CAC versus nonCAC
## C. LEfSe markers across progression groups
## D. Mean abundance heatmap of LEfSe markers
## E. LEfSe markers in CAC versus nonCAC
## F. Mediterraneibacter abundance
## G. UCG-005 abundance
## H. Lactococcus abundance
############################################################


# ==============================================================================
# 0. Required packages
# ==============================================================================

required_packages <- c(
  "ggplot2",
  "patchwork",
  "cowplot"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following R packages are missing:\n",
      paste(
        missing_packages,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(cowplot)
})


# ==============================================================================
# 1. Base directories
# ==============================================================================

base_dir <- file.path(
  PROJECT_ROOT,
  "output",
  "analysis"
)

final_panel_dir <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "Figure3_final_panels"
)

figure3A_dir <- file.path(
  final_panel_dir,
  "Figure3A_genus_composition_5groups"
)

figure3B_dir <- file.path(
  final_panel_dir,
  "Figure3B_genus_composition_CA23_vs_nonCA23"
)

figure3_progression_lefse_dir <- file.path(
  final_panel_dir,
  "Figure3_LEfSe_progression5"
)

figure3C_dir <- file.path(
  final_panel_dir,
  "Figure3C_LEfSe_CA_vs_nonCA"
)

figure3D_dir <- file.path(
  final_panel_dir,
  "Figure3D_representative_genus_trends"
)

output_dir <- file.path(
  final_panel_dir,
  "Combined_Figure3_8panels_7KB"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. RDS file paths
# ==============================================================================

figure3A_file <- file.path(
  figure3A_dir,
  "Figure3A_genus_composition_5groups_shared_colors_7KB.rds"
)

figure3B_file <- file.path(
  figure3B_dir,
  "Figure3B_genus_composition_CA23_vs_nonCA23_shared_colors_7KB.rds"
)

figure3C_file <- file.path(
  figure3_progression_lefse_dir,
  "Figure3A_LEfSe_progression5_LDA_barplot_7KB.rds"
)

figure3D_file <- file.path(
  figure3_progression_lefse_dir,
  "Figure3B_LEfSe_progression5_marker_heatmap_7KB.rds"
)

figure3E_file <- file.path(
  figure3C_dir,
  "Figure3C_LEfSe_CA_vs_nonCA_genus_markers_filtered_7KB.rds"
)

figure3F_file <- file.path(
  figure3D_dir,
  "Figure3D_panel_Mediterraneibacter_7KB.rds"
)

figure3G_file <- file.path(
  figure3D_dir,
  "Figure3D_panel_UCG_005_7KB.rds"
)

figure3H_file <- file.path(
  figure3D_dir,
  "Figure3D_panel_Lactococcus_7KB.rds"
)

rds_files <- c(
  Figure3A = figure3A_file,
  Figure3B = figure3B_file,
  Figure3C = figure3C_file,
  Figure3D = figure3D_file,
  Figure3E = figure3E_file,
  Mediterraneibacter = figure3F_file,
  UCG_005 = figure3G_file,
  Lactococcus = figure3H_file
)


# ==============================================================================
# 3. Confirm all required files exist
# ==============================================================================

missing_files <- rds_files[
  !file.exists(rds_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "The following required RDS files were not found:\n\n",
      paste(
        paste0(
          names(missing_files),
          ":\n",
          unname(missing_files)
        ),
        collapse = "\n\n"
      )
    ),
    call. = FALSE
  )
}

message("All required 7KB RDS files were found.")

message("")
message("RDS files used:")
for (nm in names(rds_files)) {
  message(
    nm,
    ": ",
    rds_files[[nm]]
  )
}


# ==============================================================================
# 4. Read plot objects
# ==============================================================================

p_A <- readRDS(
  figure3A_file
)

p_B <- readRDS(
  figure3B_file
)

p_C <- readRDS(
  figure3C_file
)

p_D <- readRDS(
  figure3D_file
)

p_E <- readRDS(
  figure3E_file
)

p_F <- readRDS(
  figure3F_file
)

p_G <- readRDS(
  figure3G_file
)

p_H <- readRDS(
  figure3H_file
)


# ==============================================================================
# 5. Panel-specific formatting for final assembly
#
# Copied from the original final manuscript assembly script.
# ==============================================================================


# ------------------------------------------------------------------------------
# Panel A: five-group composition
# Reduce legend size so that the stacked bars retain enough plotting area.
# ------------------------------------------------------------------------------

p_A <- p_A +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0,
      keyheight = grid::unit(
        3.2,
        "mm"
      ),
      keywidth = grid::unit(
        3.2,
        "mm"
      )
    )
  ) +
  theme(
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 7.5,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 6.5
    ),
    
    legend.spacing.y = grid::unit(
      0.5,
      "mm"
    ),
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 2,
      unit = "pt"
    ),
    
    plot.margin = margin(
      t = 6,
      r = 4,
      b = 4,
      l = 4,
      unit = "pt"
    )
  )


# ------------------------------------------------------------------------------
# Panel B: CAC versus nonCAC composition
# Use the same legend formatting as Panel A.
# ------------------------------------------------------------------------------

p_B <- p_B +
  guides(
    fill = guide_legend(
      title.position = "top",
      title.hjust = 0,
      keyheight = grid::unit(
        3.2,
        "mm"
      ),
      keywidth = grid::unit(
        3.2,
        "mm"
      )
    )
  ) +
  theme(
    plot.title = element_text(
      size = 10,
      face = "bold",
      hjust = 0.5
    ),
    
    axis.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      size = 7.5,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 6.5
    ),
    
    legend.spacing.y = grid::unit(
      0.5,
      "mm"
    ),
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 2,
      unit = "pt"
    ),
    
    plot.margin = margin(
      t = 6,
      r = 4,
      b = 4,
      l = 4,
      unit = "pt"
    )
  )


# ------------------------------------------------------------------------------
# Panel C: progression-group LEfSe bar plot
# ------------------------------------------------------------------------------

p_C <- p_C +
  theme(
    plot.title = element_text(
      size = 9.5,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        b = 4,
        unit = "pt"
      )
    ),
    
    axis.title = element_text(
      size = 8.5,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 7,
      color = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 7.5,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 6.5
    ),
    
    legend.key.width = grid::unit(
      3.5,
      "mm"
    ),
    
    legend.spacing.x = grid::unit(
      1.5,
      "mm"
    ),
    
    plot.margin = margin(
      t = 6,
      r = 5,
      b = 3,
      l = 5,
      unit = "pt"
    )
  )


# ------------------------------------------------------------------------------
# Panel D: LEfSe marker heatmap
# ------------------------------------------------------------------------------

p_D <- p_D +
  theme(
    plot.title = element_text(
      size = 9.5,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        b = 4,
        unit = "pt"
      )
    ),
    
    axis.title = element_blank(),
    
    axis.text.x = element_text(
      size = 7,
      color = "black",
      angle = 35,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 6.8,
      color = "black"
    ),
    
    legend.title = element_text(
      size = 7.5,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 6.5
    ),
    
    plot.margin = margin(
      t = 6,
      r = 5,
      b = 3,
      l = 5,
      unit = "pt"
    )
  )


# ------------------------------------------------------------------------------
# Panel E: CAC versus nonCAC LEfSe bar plot
# ------------------------------------------------------------------------------

p_E <- p_E +
  theme(
    plot.title = element_text(
      size = 9,
      face = "bold",
      hjust = 0.5,
      margin = margin(
        b = 4,
        unit = "pt"
      )
    ),
    
    axis.title = element_text(
      size = 8,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 6.8,
      color = "black"
    ),
    
    legend.position = "bottom",
    
    legend.title = element_text(
      size = 7,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 6.3
    ),
    
    legend.key.width = grid::unit(
      3.5,
      "mm"
    ),
    
    plot.margin = margin(
      t = 6,
      r = 5,
      b = 3,
      l = 5,
      unit = "pt"
    )
  )


# ------------------------------------------------------------------------------
# Panels F, G, H: representative-genus abundance plots
#
# Same original formatting.
# Remove repeated y-axis titles from G and H.
# ------------------------------------------------------------------------------

representative_genus_theme <- theme(
  plot.title = element_text(
    size = 9.5,
    face = "bold",
    hjust = 0.5,
    margin = margin(
      b = 4,
      unit = "pt"
    )
  ),
  
  axis.title.x = element_blank(),
  
  axis.title.y = element_text(
    size = 8,
    face = "bold"
  ),
  
  axis.text.x = element_text(
    size = 6.8,
    color = "black",
    angle = 40,
    hjust = 1
  ),
  
  axis.text.y = element_text(
    size = 7,
    color = "black"
  ),
  
  plot.margin = margin(
    t = 6,
    r = 4,
    b = 4,
    l = 4,
    unit = "pt"
  )
)

p_F <- p_F +
  representative_genus_theme +
  labs(
    y = "Relative abundance (%)"
  )

p_G <- p_G +
  representative_genus_theme +
  labs(
    y = NULL
  ) +
  theme(
    axis.title.y = element_blank()
  )

p_H <- p_H +
  representative_genus_theme +
  labs(
    y = NULL
  ) +
  theme(
    axis.title.y = element_blank()
  )


# ==============================================================================
# 5B. Force equal panel widths for F, G, and H
#
# Exact original strategy.
# ==============================================================================

aligned_representative_panels <- cowplot::align_plots(
  p_F,
  p_G,
  p_H,
  align = "hv",
  axis = "tblr"
)

representative_row <- cowplot::plot_grid(
  plotlist = aligned_representative_panels,
  nrow = 1,
  rel_widths = c(
    1,
    1,
    1
  ),
  labels = c(
    "F",
    "G",
    "H"
  ),
  label_fontfamily = "Helvetica",
  label_fontface = "bold",
  label_size = 14,
  label_x = 0,
  label_y = 1,
  hjust = 0,
  vjust = 1
)


# ==============================================================================
# 6. Define final 12-column layout
#
# Exact original layout:
#
# Row 1:
# A occupies eight columns; B occupies four columns.
#
# Row 2:
# C and D occupy six columns each.
#
# Row 3:
# E occupies four columns;
# F/G/H are contained in the remaining eight columns
# as one equal-width representative_row.
# ==============================================================================

figure_design <- "
AAAAAAAABBBB
CCCCCCDDDDDD
EEEERRRRRRRR
"


# ==============================================================================
# 7. Assemble complete figure
#
# Exact original patchwork strategy.
# ==============================================================================

combined_figure <- (
  p_A +
    p_B +
    p_C +
    p_D +
    p_E +
    wrap_elements(
      full = representative_row
    )
) +
  plot_layout(
    design = figure_design,
    heights = c(
      1.00,
      1.04,
      0.90
    )
  ) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(
        family = "Helvetica",
        face = "bold",
        size = 14,
        color = "black"
      ),
      
      plot.tag.position = c(
        0.004,
        0.996
      ),
      
      plot.margin = margin(
        t = 7,
        r = 7,
        b = 7,
        l = 7,
        unit = "pt"
      )
    )
  )


# ==============================================================================
# 8. Output filenames and dimensions
#
# Exact original dimensions.
# ==============================================================================

output_stem <- file.path(
  output_dir,
  "Figure3_8panels_final_7KB"
)

pdf_file <- paste0(
  output_stem,
  ".pdf"
)

png_file <- paste0(
  output_stem,
  ".png"
)

tiff_file <- paste0(
  output_stem,
  ".tif"
)

combined_rds_file <- paste0(
  output_stem,
  ".rds"
)

figure_width <- 16
figure_height <- 13
figure_dpi <- 600


# ==============================================================================
# 9. Export vector PDF
#
# Base PDF only; no Cairo.
# ==============================================================================

grDevices::pdf(
  file = pdf_file,
  width = figure_width,
  height = figure_height,
  family = "Helvetica",
  useDingbats = FALSE
)

print(
  combined_figure
)

grDevices::dev.off()


# ==============================================================================
# 10. Export high-resolution PNG and TIFF
#
# Same original strategy:
# ragg first; quartz fallback on macOS.
# ==============================================================================

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_png(
    filename = png_file,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    background = "white"
  )
  
  print(
    combined_figure
  )
  
  grDevices::dev.off()
  
  ragg::agg_tiff(
    filename = tiff_file,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    compression = "lzw",
    background = "white"
  )
  
  print(
    combined_figure
  )
  
  grDevices::dev.off()
  
} else {
  
  if (
    Sys.info()[["sysname"]] == "Darwin"
  ) {
    
    grDevices::png(
      filename = png_file,
      width = figure_width,
      height = figure_height,
      units = "in",
      res = figure_dpi,
      bg = "white",
      type = "quartz"
    )
    
    print(
      combined_figure
    )
    
    grDevices::dev.off()
    
    grDevices::tiff(
      filename = tiff_file,
      width = figure_width,
      height = figure_height,
      units = "in",
      res = figure_dpi,
      bg = "white",
      compression = "lzw",
      type = "quartz"
    )
    
    print(
      combined_figure
    )
    
    grDevices::dev.off()
    
  } else {
    
    stop(
      paste0(
        "Package 'ragg' is not installed and this workflow ",
        "uses quartz as the fallback only on macOS."
      ),
      call. = FALSE
    )
  }
}


# ==============================================================================
# 11. Save complete patchwork object
# ==============================================================================

saveRDS(
  combined_figure,
  combined_rds_file
)


# ==============================================================================
# 12. Save assembly manifest
# ==============================================================================

assembly_manifest <- data.frame(
  Panel = c(
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H"
  ),
  
  Content = c(
    "Genus composition across five progression groups",
    "Genus composition in CAC versus nonCAC",
    "LEfSe markers across progression groups",
    "Mean abundance heatmap of progression LEfSe markers",
    "LEfSe markers in CAC versus nonCAC",
    "Mediterraneibacter abundance across five groups",
    "UCG-005 abundance across five groups",
    "Lactococcus abundance across five groups"
  ),
  
  RDS_file = unname(
    rds_files[
      c(
        "Figure3A",
        "Figure3B",
        "Figure3C",
        "Figure3D",
        "Figure3E",
        "Mediterraneibacter",
        "UCG_005",
        "Lactococcus"
      )
    ]
  ),
  
  stringsAsFactors = FALSE
)

utils::write.table(
  assembly_manifest,
  file = file.path(
    output_dir,
    "Figure3_8panels_assembly_manifest_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# ==============================================================================
# 13. Completion messages
# ==============================================================================

message("")
message("==============================================================")
message("7KB final Figure 3 assembly completed.")
message("==============================================================")
message("")

message("Panel order:")
message("A. Five-group genus composition")
message("B. CAC versus nonCAC genus composition")
message("C. Progression LEfSe LDA")
message("D. Progression LEfSe heatmap")
message("E. CAC versus nonCAC LEfSe")
message("F. Mediterraneibacter")
message("G. UCG-005")
message("H. Lactococcus")
message("")

message(
  "Desulfovibrio was intentionally not included, ",
  "following the original final manuscript assembly."
)

message("")
message("PDF:")
message(pdf_file)

message("")
message("PNG:")
message(png_file)

message("")
message("TIFF:")
message(tiff_file)

message("")
message("Combined RDS:")
message(combined_rds_file)

message("")
message("Assembly manifest:")
message(
  file.path(
    output_dir,
    "Figure3_8panels_assembly_manifest_7KB.tsv"
  )
)