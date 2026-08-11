#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_08_plot_Figure4D_ecological_metrics.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Purpose:
## Reproduce final manuscript Figure 4D using the 7KB K=3 DMM states.
##
## Final three panels:
##   Shannon diversity
##   Observed genera
##   Dominant genus proportion
##
## Statistics:
##   Overall Kruskal-Wallis test across C1/C2/C3.
##
## The original final panel displays only the overall P value.
## No pairwise significance brackets are added.
############################################################


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "patchwork"
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
  library(ggplot2)
  library(patchwork)
})


# ==============================================================================
# 1. Paths
# ==============================================================================

project_root <- PROJECT_ROOT

input_file <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "DMM_progression127_K3_characterization_7KB",
  "tables",
  "DMM_K3_genus_diversity_and_dominance_7KB.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "Figure4_final_panels",
  "Figure4D_DMM_ecological_characteristics"
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

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found:\n",
      input_file
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 3. Read source table
# ==============================================================================

dat <- data.table::fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)


required_columns <- c(
  "SampleID",
  "DMM_cluster",
  "Genus_Shannon",
  "Observed_genera",
  "Dominant_genus_relative_abundance"
)

missing_columns <- setdiff(
  required_columns,
  colnames(
    dat
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
# 4. Prepare state factor
# ==============================================================================

cluster_order <- c(
  "C1",
  "C2",
  "C3"
)

dat <- dat %>%
  mutate(
    SampleID = as.character(
      SampleID
    ),
    
    DMM_cluster = factor(
      as.character(
        DMM_cluster
      ),
      levels = cluster_order
    )
  )


if (
  anyNA(
    dat$DMM_cluster
  )
) {
  stop(
    "Unexpected DMM state label detected.",
    call. = FALSE
  )
}


# ==============================================================================
# 5. Strict state-count audit
# ==============================================================================

state_counts <- dat %>%
  count(
    DMM_cluster,
    name = "N"
  )

expected_state_counts <- c(
  C1 = 55,
  C2 = 46,
  C3 = 26
)

observed_state_counts <- setNames(
  state_counts$N,
  as.character(
    state_counts$DMM_cluster
  )
)

if (
  !identical(
    as.integer(
      observed_state_counts[
        names(
          expected_state_counts
        )
      ]
    ),
    as.integer(
      expected_state_counts
    )
  )
) {
  stop(
    "State-count audit failed. Expected C1=55, C2=46, C3=26.",
    call. = FALSE
  )
}


# ==============================================================================
# 6. Overall Kruskal-Wallis tests
# ==============================================================================

kw_shannon <- kruskal.test(
  Genus_Shannon ~ DMM_cluster,
  data = dat
)

kw_observed <- kruskal.test(
  Observed_genera ~ DMM_cluster,
  data = dat
)

kw_dominant <- kruskal.test(
  Dominant_genus_relative_abundance ~ DMM_cluster,
  data = dat
)


statistics_table <- data.frame(
  Metric = c(
    "Shannon diversity",
    "Observed genera",
    "Dominant genus proportion"
  ),
  
  Kruskal_Wallis_chi_squared = c(
    unname(
      kw_shannon$statistic
    ),
    
    unname(
      kw_observed$statistic
    ),
    
    unname(
      kw_dominant$statistic
    )
  ),
  
  Degrees_of_freedom = c(
    unname(
      kw_shannon$parameter
    ),
    
    unname(
      kw_observed$parameter
    ),
    
    unname(
      kw_dominant$parameter
    )
  ),
  
  P_value = c(
    kw_shannon$p.value,
    kw_observed$p.value,
    kw_dominant$p.value
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  statistics_table,
  file.path(
    source_dir,
    "Figure4D_overall_Kruskal_Wallis_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 7. State summary table
# ==============================================================================

summary_table <- dat %>%
  group_by(
    DMM_cluster
  ) %>%
  summarise(
    N = n(),
    
    Shannon_median = median(
      Genus_Shannon,
      na.rm = TRUE
    ),
    
    Shannon_IQR = IQR(
      Genus_Shannon,
      na.rm = TRUE
    ),
    
    Observed_genera_median = median(
      Observed_genera,
      na.rm = TRUE
    ),
    
    Observed_genera_IQR = IQR(
      Observed_genera,
      na.rm = TRUE
    ),
    
    Dominant_genus_proportion_median = median(
      Dominant_genus_relative_abundance,
      na.rm = TRUE
    ),
    
    Dominant_genus_proportion_IQR = IQR(
      Dominant_genus_relative_abundance,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


data.table::fwrite(
  summary_table,
  file.path(
    source_dir,
    "Figure4D_state_summary_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 8. Plot settings
#
# Follow final manuscript panel:
# - grey boxplots
# - individual sample points
# - C1/C2/C3 on x-axis
# - overall P value only
# ==============================================================================

base_family <- "Helvetica"

base_size <- 8

box_fill <- "grey85"

box_line <- "black"

point_alpha <- 0.55

point_size <- 1.3

jitter_width <- 0.12


format_overall_p <- function(p_value) {
  
  if (
    is.na(
      p_value
    )
  ) {
    return(
      "P = NA"
    )
  }
  
  if (
    p_value < 0.001
  ) {
    return(
      "P < 0.001"
    )
  }
  
  paste0(
    "P = ",
    sprintf(
      "%.3f",
      p_value
    )
  )
}


# ==============================================================================
# 9. Common theme
# ==============================================================================

theme_dmm_metric <- theme_classic(
  base_size = base_size,
  base_family = base_family
) +
  theme(
    text = element_text(
      family = base_family,
      color = "black"
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 8,
      face = "plain",
      color = "black"
    ),
    
    axis.text.x = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 7.5,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.35,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.35,
      color = "black"
    ),
    
    plot.margin = margin(
      5,
      4,
      4,
      4,
      unit = "pt"
    )
  )


# ==============================================================================
# 10. Panel D1: Shannon diversity
# ==============================================================================

set.seed(
  20260726
)

shannon_y_max <- max(
  dat$Genus_Shannon,
  na.rm = TRUE
)

shannon_y_min <- min(
  dat$Genus_Shannon,
  na.rm = TRUE
)

shannon_range <- shannon_y_max -
  shannon_y_min

if (
  shannon_range <= 0
) {
  shannon_range <- 1
}

shannon_annotation_y <-
  shannon_y_max +
  0.08 * shannon_range


p_shannon <- ggplot(
  dat,
  aes(
    x = DMM_cluster,
    y = Genus_Shannon
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    fill = box_fill,
    color = box_line,
    linewidth = 0.38
  ) +
  geom_jitter(
    width = jitter_width,
    height = 0,
    size = point_size,
    alpha = point_alpha
  ) +
  annotate(
    "text",
    x = 2,
    y = shannon_annotation_y,
    label = format_overall_p(
      kw_shannon$p.value
    ),
    family = base_family,
    size = 2.8
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.15
      )
    )
  ) +
  labs(
    y = "Shannon diversity"
  ) +
  theme_dmm_metric


# ==============================================================================
# 11. Panel D2: Observed genera
# ==============================================================================

set.seed(
  20260727
)

observed_y_max <- max(
  dat$Observed_genera,
  na.rm = TRUE
)

observed_y_min <- min(
  dat$Observed_genera,
  na.rm = TRUE
)

observed_range <- observed_y_max -
  observed_y_min

if (
  observed_range <= 0
) {
  observed_range <- 1
}

observed_annotation_y <-
  observed_y_max +
  0.08 * observed_range


p_observed <- ggplot(
  dat,
  aes(
    x = DMM_cluster,
    y = Observed_genera
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    fill = box_fill,
    color = box_line,
    linewidth = 0.38
  ) +
  geom_jitter(
    width = jitter_width,
    height = 0,
    size = point_size,
    alpha = point_alpha
  ) +
  annotate(
    "text",
    x = 2,
    y = observed_annotation_y,
    label = format_overall_p(
      kw_observed$p.value
    ),
    family = base_family,
    size = 2.8
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.04,
        0.15
      )
    )
  ) +
  labs(
    y = "Observed genera"
  ) +
  theme_dmm_metric


# ==============================================================================
# 12. Panel D3: Dominant genus proportion
# ==============================================================================

set.seed(
  20260728
)

dominant_y_max <- max(
  dat$Dominant_genus_relative_abundance,
  na.rm = TRUE
)

dominant_y_min <- min(
  dat$Dominant_genus_relative_abundance,
  na.rm = TRUE
)

dominant_range <- dominant_y_max -
  dominant_y_min

if (
  dominant_range <= 0
) {
  dominant_range <- 1
}

dominant_annotation_y <-
  dominant_y_max +
  0.08 * dominant_range


p_dominant <- ggplot(
  dat,
  aes(
    x = DMM_cluster,
    y = Dominant_genus_relative_abundance
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    fill = box_fill,
    color = box_line,
    linewidth = 0.38
  ) +
  geom_jitter(
    width = jitter_width,
    height = 0,
    size = point_size,
    alpha = point_alpha
  ) +
  annotate(
    "text",
    x = 2,
    y = dominant_annotation_y,
    label = format_overall_p(
      kw_dominant$p.value
    ),
    family = base_family,
    size = 2.8
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(
        0.04,
        0.15
      )
    )
  ) +
  labs(
    y = "Dominant genus proportion"
  ) +
  theme_dmm_metric


# ==============================================================================
# 13. Combine three panels
#
# Final Figure 4D arrangement:
# Shannon | Observed genera | Dominant genus proportion
# ==============================================================================

panel_d <- (
  p_shannon |
    p_observed |
    p_dominant
) +
  plot_layout(
    widths = c(
      1,
      1,
      1
    )
  )


# ==============================================================================
# 14. Export plotting source values
# ==============================================================================

plot_source <- dat %>%
  select(
    SampleID,
    DMM_cluster,
    Genus_Shannon,
    Observed_genera,
    Dominant_genus_relative_abundance
  )


data.table::fwrite(
  plot_source,
  file.path(
    source_dir,
    "Figure4D_plot_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 15. Save individual plot RDS objects
# ==============================================================================

saveRDS(
  p_shannon,
  file.path(
    output_dir,
    "Figure4D_Shannon_7KB.rds"
  )
)

saveRDS(
  p_observed,
  file.path(
    output_dir,
    "Figure4D_Observed_genera_7KB.rds"
  )
)

saveRDS(
  p_dominant,
  file.path(
    output_dir,
    "Figure4D_Dominant_genus_proportion_7KB.rds"
  )
)

saveRDS(
  panel_d,
  file.path(
    output_dir,
    "Figure4D_DMM_ecological_characteristics_7KB.rds"
  )
)


# ==============================================================================
# 16. Export PDF
#
# Base PDF only; no Cairo.
# ==============================================================================

pdf_file <- file.path(
  output_dir,
  "Figure4D_DMM_ecological_characteristics_7KB.pdf"
)

grDevices::pdf(
  file = pdf_file,
  width = 7.2,
  height = 3.0,
  family = "Helvetica",
  useDingbats = FALSE
)

print(
  panel_d
)

grDevices::dev.off()


# ==============================================================================
# 17. Export PNG
# ==============================================================================

png_file <- file.path(
  output_dir,
  "Figure4D_DMM_ecological_characteristics_7KB.png"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_png(
    filename = png_file,
    width = 7.2,
    height = 3.0,
    units = "in",
    res = 600,
    background = "white"
  )
  
  print(
    panel_d
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::png(
    filename = png_file,
    width = 7.2,
    height = 3.0,
    units = "in",
    res = 600,
    type = "quartz",
    bg = "white"
  )
  
  print(
    panel_d
  )
  
  grDevices::dev.off()
}


# ==============================================================================
# 18. Export TIFF
# ==============================================================================

tiff_file <- file.path(
  output_dir,
  "Figure4D_DMM_ecological_characteristics_7KB.tif"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_tiff(
    filename = tiff_file,
    width = 7.2,
    height = 3.0,
    units = "in",
    res = 600,
    compression = "lzw",
    background = "white"
  )
  
  print(
    panel_d
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::tiff(
    filename = tiff_file,
    width = 7.2,
    height = 3.0,
    units = "in",
    res = 600,
    type = "quartz",
    compression = "lzw",
    bg = "white"
  )
  
  print(
    panel_d
  )
  
  grDevices::dev.off()
}


# ==============================================================================
# 19. Console summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("Final 7KB Figure 4D completed.\n")
cat("============================================================\n")

cat("\n")
cat("State counts:\n")
print(
  state_counts
)

cat("\n")
cat("State summaries:\n")
print(
  tibble::as_tibble(
    summary_table
  )
)

cat("\n")
cat("Overall Kruskal-Wallis tests:\n")
print(
  tibble::as_tibble(
    statistics_table
  )
)

cat("\n")
cat("Final displayed P labels:\n")

cat(
  "Shannon diversity: ",
  format_overall_p(
    kw_shannon$p.value
  ),
  "\n",
  sep = ""
)

cat(
  "Observed genera: ",
  format_overall_p(
    kw_observed$p.value
  ),
  "\n",
  sep = ""
)

cat(
  "Dominant genus proportion: ",
  format_overall_p(
    kw_dominant$p.value
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