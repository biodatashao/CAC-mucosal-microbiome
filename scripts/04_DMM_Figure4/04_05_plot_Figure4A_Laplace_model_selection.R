#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_05_plot_Figure4A_Laplace_model_selection.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Purpose:
## Generate Figure 4A for the final 7KB DMM rerun.
##
## These are the locked best-converged Laplace values from the completed
## K = 1:7 x 20-repeat primary DMM analysis.
##
## No model fitting is performed here.
## This script is plotting only.
############################################################

options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "scales"
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
  library(data.table)
  library(ggplot2)
  library(scales)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "Figure4_final_panels",
  "Figure4A_DMM_Laplace_model_selection"
)

source_dir <- file.path(
  output_dir,
  "source_tables"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  source_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Locked best-repeat results from primary 7KB DMM
# ==============================================================================

laplace_data <- data.frame(
  K = 1:7,
  
  Repeat = c(
    1,
    13,
    10,
    16,
    12,
    9,
    11
  ),
  
  Seed = c(
    20270716,
    20280728,
    20290725,
    20300731,
    20310727,
    20320724,
    20330726
  ),
  
  Laplace = c(
    90707,
    87151,
    86176,
    86185,
    86325,
    86635,
    87044
  ),
  
  AIC = c(
    90764,
    87548,
    87056,
    87531,
    88416,
    89315,
    90182
  ),
  
  BIC = c(
    91187,
    88395,
    88325,
    89225,
    90533,
    91856,
    93147
  ),
  
  stringsAsFactors = FALSE
)


# ==============================================================================
# 3. Audit optimum
# ==============================================================================

optimal_row <- laplace_data[
  which.min(
    laplace_data$Laplace
  ),
  ,
  drop = FALSE
]

optimal_k <- optimal_row$K
optimal_laplace <- optimal_row$Laplace

if (optimal_k != 3) {
  stop(
    paste0(
      "Optimal K audit failed. Expected K=3, observed K=",
      optimal_k,
      "."
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 4. Save exact plotting source
# ==============================================================================

data.table::fwrite(
  laplace_data,
  file.path(
    source_dir,
    "Figure4A_DMM_best_repeat_model_selection_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 5. Plot settings
#
# Reproduces the original final Figure 4A structure:
# - K = 1 to 7
# - Laplace approximation on y axis
# - connected black points
# - optimal K=3 highlighted by text annotation
#
# Final panel tag A and final typography are applied again in the
# full A-E assembly script.
# ==============================================================================

base_family <- "Helvetica"


panel_a <- ggplot2::ggplot(
  laplace_data,
  ggplot2::aes(
    x = K,
    y = Laplace
  )
) +
  
  ggplot2::geom_line(
    linewidth = 0.55,
    color = "black"
  ) +
  
  ggplot2::geom_point(
    size = 2.0,
    shape = 16,
    color = "black"
  ) +
  
  ggplot2::geom_point(
    data = optimal_row,
    size = 2.6,
    shape = 16,
    color = "black"
  ) +
  
  ggplot2::annotate(
    "text",
    x = optimal_k,
    y = optimal_laplace - 230,
    label = "K = 3",
    family = base_family,
    fontface = "bold",
    size = 3.0,
    hjust = 0.5,
    vjust = 1
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = 1:7,
    limits = c(
      0.8,
      7.2
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  ggplot2::scale_y_continuous(
    labels = scales::label_comma(
      accuracy = 1
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.08,
        0.05
      )
    )
  ) +
  
  ggplot2::labs(
    x = "Number of DMM components (K)",
    y = "Laplace approximation"
  ) +
  
  ggplot2::theme_classic(
    base_size = 9,
    base_family = base_family
  ) +
  
  ggplot2::theme(
    text = ggplot2::element_text(
      family = base_family,
      color = "black"
    ),
    
    axis.title = ggplot2::element_text(
      family = base_family,
      size = 9,
      face = "bold",
      color = "black"
    ),
    
    axis.text = ggplot2::element_text(
      family = base_family,
      size = 8,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      linewidth = 0.4,
      color = "black"
    ),
    
    axis.ticks = ggplot2::element_line(
      linewidth = 0.4,
      color = "black"
    ),
    
    plot.margin = ggplot2::margin(
      t = 6,
      r = 6,
      b = 6,
      l = 6,
      unit = "pt"
    )
  )


# ==============================================================================
# 6. Save RDS
#
# This is what the final Figure 4 assembly will read.
# ==============================================================================

rds_file <- file.path(
  output_dir,
  "Figure4A_DMM_Laplace_model_selection_7KB.rds"
)

saveRDS(
  panel_a,
  rds_file
)


# ==============================================================================
# 7. Export PDF
# ==============================================================================

pdf_file <- file.path(
  output_dir,
  "Figure4A_DMM_Laplace_model_selection_7KB.pdf"
)

grDevices::pdf(
  file = pdf_file,
  width = 4.1,
  height = 3.3,
  family = base_family,
  useDingbats = FALSE
)

print(
  panel_a
)

grDevices::dev.off()


# ==============================================================================
# 8. Export PNG
# ==============================================================================

png_file <- file.path(
  output_dir,
  "Figure4A_DMM_Laplace_model_selection_7KB.png"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_png(
    filename = png_file,
    width = 4.1,
    height = 3.3,
    units = "in",
    res = 600,
    background = "white"
  )
  
} else {
  
  grDevices::png(
    filename = png_file,
    width = 4.1,
    height = 3.3,
    units = "in",
    res = 600,
    type = "quartz",
    bg = "white"
  )
}

print(
  panel_a
)

grDevices::dev.off()


# ==============================================================================
# 9. Export TIFF
# ==============================================================================

tiff_file <- file.path(
  output_dir,
  "Figure4A_DMM_Laplace_model_selection_7KB.tif"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_tiff(
    filename = tiff_file,
    width = 4.1,
    height = 3.3,
    units = "in",
    res = 600,
    compression = "lzw",
    background = "white"
  )
  
} else {
  
  grDevices::tiff(
    filename = tiff_file,
    width = 4.1,
    height = 3.3,
    units = "in",
    res = 600,
    compression = "lzw",
    type = "quartz",
    bg = "white"
  )
}

print(
  panel_a
)

grDevices::dev.off()


# ==============================================================================
# 10. Console summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("Final 7KB Figure 4A completed successfully.\n")
cat("============================================================\n")

cat("\nBest model per K:\n")

print(
  laplace_data,
  row.names = FALSE
)

cat("\nOptimal model:\n")

cat(
  "K = ",
  optimal_k,
  "\n",
  sep = ""
)

cat(
  "Laplace = ",
  optimal_laplace,
  "\n",
  sep = ""
)

cat("\nRDS for final Figure 4 assembly:\n")

cat(
  rds_file,
  "\n"
)

cat("\nOutput directory:\n")

cat(
  output_dir,
  "\n"
)