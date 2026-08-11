#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_06_plot_Figure4B_state_distribution.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Purpose:
## Reproduce final manuscript Figure 4B using the 7KB K=3 DMM assignments.
##
## Panel:
## 100% stacked bar chart showing C1/C2/C3 proportions across:
## Polyp / UC remission / UC active / Dysplasia / CAC
##
## CAC C2 label:
## 17/23
## 73.9%
############################################################


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
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
      paste(
        missing_packages,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

assignment_file <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "DMM_progression127_primary_prevalence10_7KB",
  "tables",
  "DMM_sample_assignments_optimalK3_native_7KB.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "Figure4_final_panels",
  "Figure4B_DMM_state_distribution"
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
# 2. Validate input
# ==============================================================================

if (!file.exists(assignment_file)) {
  stop(
    paste0(
      "Assignment file not found:\n",
      assignment_file
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 3. Read assignments
# ==============================================================================

assignment_raw <- data.table::fread(
  assignment_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "SampleID",
  "DMM_cluster",
  "Progression5"
)

missing_columns <- setdiff(
  required_columns,
  colnames(
    assignment_raw
  )
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 4. Factors and labels
# ==============================================================================

group_order <- c(
  "Polyp",
  "UC_remission",
  "UC_active",
  "Dysplasia",
  "CA"
)

group_display <- c(
  Polyp = "Polyp",
  UC_remission = "UC remission",
  UC_active = "UC active",
  Dysplasia = "Dysplasia",
  CA = "CAC"
)

cluster_order <- c(
  "C1",
  "C2",
  "C3"
)

assignment <- assignment_raw %>%
  transmute(
    SampleID = as.character(
      SampleID
    ),
    
    Progression5 = factor(
      as.character(
        Progression5
      ),
      levels = group_order
    ),
    
    DMM_cluster = factor(
      as.character(
        DMM_cluster
      ),
      levels = cluster_order
    )
  )


# ==============================================================================
# 5. Strict sample-count audit
# ==============================================================================

expected_group_counts <- c(
  Polyp = 26,
  UC_remission = 36,
  UC_active = 25,
  Dysplasia = 17,
  CA = 23
)

observed_group_counts <- assignment %>%
  count(
    Progression5,
    name = "N"
  )

observed_group_named <- setNames(
  observed_group_counts$N,
  as.character(
    observed_group_counts$Progression5
  )
)

if (
  !identical(
    as.integer(
      observed_group_named[
        names(
          expected_group_counts
        )
      ]
    ),
    as.integer(
      expected_group_counts
    )
  )
) {
  stop(
    "Progression group count audit failed.",
    call. = FALSE
  )
}


expected_cluster_counts <- c(
  C1 = 55,
  C2 = 46,
  C3 = 26
)

observed_cluster_counts <- assignment %>%
  count(
    DMM_cluster,
    name = "N"
  )

observed_cluster_named <- setNames(
  observed_cluster_counts$N,
  as.character(
    observed_cluster_counts$DMM_cluster
  )
)

if (
  !identical(
    as.integer(
      observed_cluster_named[
        names(
          expected_cluster_counts
        )
      ]
    ),
    as.integer(
      expected_cluster_counts
    )
  )
) {
  stop(
    "DMM state count audit failed.",
    call. = FALSE
  )
}


# ==============================================================================
# 6. Build complete state-distribution table
# ==============================================================================

plot_table <- assignment %>%
  count(
    Progression5,
    DMM_cluster,
    name = "N"
  ) %>%
  tidyr::complete(
    Progression5 = factor(
      group_order,
      levels = group_order
    ),
    
    DMM_cluster = factor(
      cluster_order,
      levels = cluster_order
    ),
    
    fill = list(
      N = 0L
    )
  ) %>%
  group_by(
    Progression5
  ) %>%
  mutate(
    Total = sum(
      N
    ),
    
    Proportion =
      N /
      Total
  ) %>%
  ungroup() %>%
  mutate(
    Progression_display = factor(
      group_display[
        as.character(
          Progression5
        )
      ],
      levels = unname(
        group_display[
          group_order
        ]
      )
    )
  )


# ==============================================================================
# 7. Exact CAC C2 audit
# ==============================================================================

cac_c2 <- plot_table %>%
  filter(
    Progression5 == "CA",
    DMM_cluster == "C2"
  )

if (
  nrow(
    cac_c2
  ) != 1 ||
  cac_c2$N != 17 ||
  cac_c2$Total != 23
) {
  stop(
    paste0(
      "CAC C2 audit failed. Expected 17/23, observed ",
      paste0(
        cac_c2$N,
        "/",
        cac_c2$Total
      )
    ),
    call. = FALSE
  )
}

cac_c2_percent <- 100 *
  cac_c2$Proportion

if (
  abs(
    cac_c2_percent -
    73.9130434783
  ) > 1e-6
) {
  stop(
    "CAC C2 percentage audit failed.",
    call. = FALSE
  )
}


# ==============================================================================
# 8. Export source table
# ==============================================================================

data.table::fwrite(
  plot_table,
  file.path(
    source_dir,
    "Figure4B_DMM_state_distribution_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 9. State colors
#
# Keep consistent throughout Figure 4.
# ==============================================================================

state_colors <- c(
  "C1" = "#6C8EBF",
  "C2" = "#C95E57",
  "C3" = "#88B27C"
)


# ==============================================================================
# 10. Label position for CAC C2
# ==============================================================================

cac_c2_label_data <- plot_table %>%
  filter(
    Progression5 == "CA",
    DMM_cluster == "C2"
  ) %>%
  mutate(
    Label = paste0(
      N,
      "/",
      Total,
      "\n",
      sprintf(
        "%.1f%%",
        100 * Proportion
      )
    )
  )


# ==============================================================================
# 11. Final Figure 4B
# ==============================================================================

panel_b <- ggplot(
  plot_table,
  aes(
    x = Progression_display,
    y = Proportion,
    fill = DMM_cluster
  )
) +
  geom_col(
    width = 0.72
  ) +
  geom_text(
    data = cac_c2_label_data,
    aes(
      x = Progression_display,
      y = 0.56,
      label = Label
    ),
    inherit.aes = FALSE,
    family = "Helvetica",
    size = 2.7,
    lineheight = 0.95,
    color = "black"
  ) +
  scale_fill_manual(
    values = state_colors,
    breaks = cluster_order,
    drop = FALSE,
    name = "DMM state"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = c(
      0,
      0.25,
      0.50,
      0.75,
      1
    ),
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = c(
      0,
      0
    )
  ) +
  labs(
    x = NULL,
    y = "Proportion of samples"
  ) +
  theme_classic(
    base_size = 8,
    base_family = "Helvetica"
  ) +
  theme(
    text = element_text(
      family = "Helvetica",
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 8.5,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 8,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    
    axis.line = element_line(
      linewidth = 0.35,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.35,
      color = "black"
    ),
    
    legend.position = "bottom",
    
    legend.direction = "horizontal",
    
    legend.title = element_text(
      size = 8,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.width = grid::unit(
      4.0,
      "mm"
    ),
    
    legend.key.height = grid::unit(
      3.0,
      "mm"
    ),
    
    plot.margin = margin(
      4,
      5,
      3,
      5,
      unit = "pt"
    )
  )


# ==============================================================================
# 12. Save RDS
# ==============================================================================

saveRDS(
  panel_b,
  file.path(
    output_dir,
    "Figure4B_DMM_state_distribution_7KB.rds"
  )
)


# ==============================================================================
# 13. Export PDF
# ==============================================================================

grDevices::pdf(
  file = file.path(
    output_dir,
    "Figure4B_DMM_state_distribution_7KB.pdf"
  ),
  width = 5.5,
  height = 3.7,
  family = "Helvetica",
  useDingbats = FALSE
)

print(
  panel_b
)

grDevices::dev.off()


# ==============================================================================
# 14. Export PNG
# ==============================================================================

png_file <- file.path(
  output_dir,
  "Figure4B_DMM_state_distribution_7KB.png"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_png(
    filename = png_file,
    width = 5.5,
    height = 3.7,
    units = "in",
    res = 600,
    background = "white"
  )
  
  print(
    panel_b
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::png(
    filename = png_file,
    width = 5.5,
    height = 3.7,
    units = "in",
    res = 600,
    type = "quartz",
    bg = "white"
  )
  
  print(
    panel_b
  )
  
  grDevices::dev.off()
}


# ==============================================================================
# 15. Export TIFF
# ==============================================================================

tiff_file <- file.path(
  output_dir,
  "Figure4B_DMM_state_distribution_7KB.tif"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_tiff(
    filename = tiff_file,
    width = 5.5,
    height = 3.7,
    units = "in",
    res = 600,
    compression = "lzw",
    background = "white"
  )
  
  print(
    panel_b
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::tiff(
    filename = tiff_file,
    width = 5.5,
    height = 3.7,
    units = "in",
    res = 600,
    type = "quartz",
    compression = "lzw",
    bg = "white"
  )
  
  print(
    panel_b
  )
  
  grDevices::dev.off()
}


# ==============================================================================
# 16. Console summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("Final 7KB Figure 4B completed.\n")
cat("============================================================\n")

cat("\n")
cat("DMM state distribution by progression group:\n")

print(
  plot_table %>%
    select(
      Progression5,
      DMM_cluster,
      N,
      Total,
      Proportion
    )
)

cat("\n")
cat("CAC C2:\n")
cat(
  cac_c2$N,
  "/",
  cac_c2$Total,
  " = ",
  sprintf(
    "%.1f%%",
    100 * cac_c2$Proportion
  ),
  "\n",
  sep = ""
)

cat("\n")
cat("Output directory:\n")
cat(
  output_dir,
  "\n"
)
