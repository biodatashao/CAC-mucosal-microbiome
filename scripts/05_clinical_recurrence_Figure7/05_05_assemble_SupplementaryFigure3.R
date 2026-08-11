#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 05_05_assemble_SupplementaryFigure3.R
##
## Module 05 - Clinical association and recurrence (Figure 7, Supplementary Figure 3)
############################################################
# Plot style is retained from the original manuscript script.
# Only:
#   1. input source is replaced by the 7KB source table;
#   2. current column names are adapted;
#   3. marker order is locked to the ORIGINAL final Supplementary Figure 3.
#
# Main Figure 7 contains:
#   Lactococcus
#   UCG-005
#
# Supplementary Figure 3 contains:
#   Mediterraneibacter
#   Candidatus Soleaferrea
#   Peptoclostridium
#   Monoglobus
#   Dorea
#   Atopostipes
#   [Eubacterium] eligens group
#   Desulfovibrio
# ==============================================================================


options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "ggplot2"
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
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

rerun_root <- file.path(
  project_root,
  "output",
  "analysis"
)

input_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "01_core_analysis",
  "Supplementary_Figure3"
)

marker_file <- file.path(
  input_dir,
  "Supplementary_Figure3_marker_abundance_source_7KB.tsv"
)

out_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "04_final_Supplementary_Figure3"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Input check
# ==============================================================================

if (!file.exists(marker_file)) {
  stop(
    paste0(
      "Marker source file does not exist:\n",
      marker_file
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 3. ORIGINAL style constants
# ==============================================================================

font_family <- "Helvetica"


# ==============================================================================
# 4. Marker label formatting
# Retains original manuscript labels.
# ==============================================================================

format_marker <- function(x) {
  
  dplyr::recode(
    as.character(x),
    
    "UCG_005" =
      "UCG-005",
    
    "Candidatus_Soleaferrea" =
      "Candidatus Soleaferrea",
    
    "Eubacterium_eligens_group" =
      "[Eubacterium] eligens group",
    
    "f__Eubacterium__eligens_group" =
      "[Eubacterium] eligens group",
    
    .default = as.character(x)
  )
}


# ==============================================================================
# 5. P-value formatting
#
# Visual final manuscript figure uses capital P.
# Numeric precision otherwise follows the values reported in the manuscript:
#   signif(p, 2)
# ==============================================================================

format_p <- function(p) {
  
  if (is.na(p)) {
    return("P = NA")
  }
  
  if (p < 0.001) {
    return(
      paste0(
        "P = ",
        formatC(
          p,
          format = "e",
          digits = 2
        )
      )
    )
  }
  
  paste0(
    "P = ",
    signif(
      p,
      2
    )
  )
}


# ==============================================================================
# 6. Read 7KB marker source
# ==============================================================================

dat0 <- readr::read_tsv(
  marker_file,
  show_col_types = FALSE
)


required_columns <- c(
  "Genus",
  "Recurrence_plot",
  "Relative_abundance_percent"
)

missing_columns <- setdiff(
  required_columns,
  colnames(dat0)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "7KB marker source is missing required column(s): ",
      paste(
        missing_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# Adapt current 7KB source column to original script naming.
dat0 <- dat0 %>%
  mutate(
    Marker_label = as.character(Genus)
  )


# ==============================================================================
# 7. Lock ORIGINAL Supplementary Figure 3 marker order
#
# Do NOT reorder by the new 7KB P values.
# This preserves the original 4 x 2 layout.
# ==============================================================================

marker_order_raw <- c(
  "Mediterraneibacter",
  "Candidatus_Soleaferrea",
  "Peptoclostridium",
  "Monoglobus",
  "Dorea",
  "Atopostipes",
  "Eubacterium_eligens_group",
  "Desulfovibrio"
)

marker_order_plot <- format_marker(
  marker_order_raw
)


missing_markers <- setdiff(
  marker_order_raw,
  unique(
    as.character(
      dat0$Marker_label
    )
  )
)

if (length(missing_markers) > 0) {
  stop(
    paste0(
      "Required Supplementary Figure 3 marker(s) absent from source:\n",
      paste(
        missing_markers,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 8. Build plot data
# ==============================================================================

plot_dat <- dat0 %>%
  filter(
    Marker_label %in%
      marker_order_raw
  ) %>%
  mutate(
    Recurrence_plot = factor(
      Recurrence_plot,
      levels = c(
        "No recurrence",
        "Recurrence"
      )
    ),
    
    Marker_label_plot = format_marker(
      Marker_label
    ),
    
    Marker_label_plot = factor(
      Marker_label_plot,
      levels = marker_order_plot
    ),
    
    Relative_abundance_percent = as.numeric(
      Relative_abundance_percent
    )
  ) %>%
  filter(
    !is.na(Recurrence_plot),
    !is.na(Relative_abundance_percent)
  )


# ==============================================================================
# 9. Recalculate Wilcoxon statistics exactly as original script
# ==============================================================================

stats_tbl <- plot_dat %>%
  group_by(
    Marker_label,
    Marker_label_plot
  ) %>%
  summarise(
    p_value = if (
      n_distinct(
        Recurrence_plot
      ) == 2
    ) {
      
      suppressWarnings(
        wilcox.test(
          Relative_abundance_percent ~ Recurrence_plot,
          exact = FALSE
        )$p.value
      )
      
    } else {
      
      NA_real_
    },
    
    .groups = "drop"
  ) %>%
  mutate(
    P_label = vapply(
      p_value,
      format_p,
      character(1)
    )
  )


# Preserve original display order.
stats_tbl$Marker_label_plot <- factor(
  stats_tbl$Marker_label_plot,
  levels = marker_order_plot
)

stats_tbl <- stats_tbl %>%
  arrange(
    Marker_label_plot
  )


# ==============================================================================
# 10. Attach P labels
# ==============================================================================

plot_dat <- plot_dat %>%
  left_join(
    stats_tbl %>%
      select(
        Marker_label,
        P_label
      ),
    by = "Marker_label"
  )


p_label_df <- plot_dat %>%
  distinct(
    Marker_label_plot,
    P_label
  )


# ==============================================================================
# 11. ORIGINAL recurrence colors
# ==============================================================================

recurrence_colors <- c(
  "No recurrence" = "#A99E79",
  "Recurrence" = "#A32635"
)


# ==============================================================================
# 12. Plot
#
# Geometry/theme parameters retained from original manuscript script.
# ==============================================================================

p <- ggplot(
  plot_dat,
  aes(
    x = Recurrence_plot,
    y = Relative_abundance_percent
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    color = "black",
    fill = "white",
    linewidth = 0.55
  ) +
  geom_jitter(
    aes(
      color = Recurrence_plot
    ),
    width = 0.11,
    size = 1.7,
    alpha = 0.82,
    stroke = 0
  ) +
  geom_text(
    data = p_label_df,
    aes(
      x = Inf,
      y = Inf,
      label = P_label
    ),
    inherit.aes = FALSE,
    hjust = 1.05,
    vjust = 1.35,
    size = 3.0,
    fontface = "bold",
    family = font_family
  ) +
  scale_color_manual(
    values = recurrence_colors,
    drop = FALSE
  ) +
  facet_wrap(
    ~ Marker_label_plot,
    scales = "free_y",
    ncol = 4
  ) +
  labs(
    x = NULL,
    y = "Relative abundance (%)"
  ) +
  theme_classic(
    base_size = 8,
    base_family = font_family
  ) +
  theme(
    text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      size = 9,
      face = "bold"
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
    
    axis.line = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    strip.text = element_text(
      size = 8.5,
      face = "bold"
    ),
    
    strip.background = element_blank(),
    
    panel.grid = element_blank(),
    
    legend.position = "none",
    
    plot.margin = margin(
      6,
      6,
      6,
      6
    )
  )


# ==============================================================================
# 13. Export
#
# Original dimensions:
#   8.6 x 5.2 inches
#   600 dpi
#
# PDF uses base grDevices::pdf() — no Cairo.
# ==============================================================================

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
      ".tif"
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # PDF
  # ---------------------------------------------------------------------------
  
  grDevices::pdf(
    pdf_file,
    width = width,
    height = height,
    family = font_family,
    useDingbats = FALSE
  )
  
  print(
    plot
  )
  
  grDevices::dev.off()
  
  
  # ---------------------------------------------------------------------------
  # PNG
  # ---------------------------------------------------------------------------
  
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
    
  } else {
    
    grDevices::png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white",
      type = "quartz"
    )
  }
  
  print(
    plot
  )
  
  grDevices::dev.off()
  
  
  # ---------------------------------------------------------------------------
  # TIFF
  # ---------------------------------------------------------------------------
  
  if (
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ) {
    
    ragg::agg_tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      background = "white"
    )
    
  } else {
    
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
  }
  
  print(
    plot
  )
  
  grDevices::dev.off()
  
  
  invisible(
    list(
      PDF = pdf_file,
      PNG = png_file,
      TIFF = tif_file
    )
  )
}


output_files <- export_plot(
  p,
  "Supplementary_Figure3_CA_recurrence_additional_LEfSe_marker_abundance_7KB",
  width = 8.6,
  height = 5.2
)


# ==============================================================================
# 14. Save source tables
# ==============================================================================

readr::write_csv(
  plot_dat,
  file.path(
    out_dir,
    "Supplementary_Figure3_CA_recurrence_additional_LEfSe_marker_abundance_plot_data_7KB.csv"
  )
)

readr::write_csv(
  stats_tbl,
  file.path(
    out_dir,
    "Supplementary_Figure3_CA_recurrence_additional_LEfSe_marker_abundance_statistics_7KB.csv"
  )
)

saveRDS(
  p,
  file.path(
    out_dir,
    "Supplementary_Figure3_CA_recurrence_additional_LEfSe_marker_abundance_7KB.rds"
  )
)


# ==============================================================================
# 15. Console summary
# ==============================================================================

message("")
message("==============================================================")
message("Supplementary Figure 3 completed.")
message("==============================================================")

message("")
message("Fixed marker order:")

for (current_marker in marker_order_plot) {
  message("  ", current_marker)
}

message("")
message("7KB statistics:")

for (i in seq_len(nrow(stats_tbl))) {
  
  message(
    "  ",
    as.character(
      stats_tbl$Marker_label_plot[i]
    ),
    ": P = ",
    signif(
      stats_tbl$p_value[i],
      6
    )
  )
}

message("")
message("PDF: ", output_files$PDF)
message("PNG: ", output_files$PNG)
message("TIFF: ", output_files$TIFF)
message("Output directory: ", out_dir)