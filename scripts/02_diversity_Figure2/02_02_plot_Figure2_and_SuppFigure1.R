#!/usr/bin/env Rscript




############################################################
## 02_02_plot_Figure2_and_SuppFigure1.R
##
## Module 02 - Alpha and beta diversity (Figure 2, Supplementary Figure 1)
##
## Final 7-KB manuscript Figure 2 + Supplementary Figure 1
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
  library(stringr)
  library(ggplot2)
  library(cowplot)
})


############################################################
## 0. Paths
############################################################

base_dir <- file.path(PROJECT_ROOT, "output/analysis")


stat_dir <- file.path(
  base_dir,
  "01_Figure2_alpha_beta",
  "statistics"
)


out_dir <- file.path(
  base_dir,
  "01_Figure2_alpha_beta",
  "figures"
)


dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 1. Input files
############################################################

f_alpha_5g <- file.path(
  stat_dir,
  "Table_01_progression127_alpha_by_sample.csv"
)


f_kw_5g <- file.path(
  stat_dir,
  "Table_03_progression127_alpha_Kruskal_Wallis.csv"
)


f_pcoa_5g <- file.path(
  stat_dir,
  "Table_05_progression127_Bray_PCoA_coordinates.csv"
)


f_pcoa_var_5g <- file.path(
  stat_dir,
  "Table_06_progression127_Bray_PCoA_variance.csv"
)


f_perm_5g <- file.path(
  stat_dir,
  "Table_07_progression127_PERMANOVA.csv"
)


f_alpha_pair <- file.path(
  stat_dir,
  "Table_10_paired_CA23_nonCA23_alpha_by_sample.csv"
)


f_wil_pair <- file.path(
  stat_dir,
  "Table_12_paired_CA23_nonCA23_Wilcoxon.csv"
)


f_pcoa_pair <- file.path(
  stat_dir,
  "Table_13_paired_CA23_nonCA23_Bray_PCoA_coordinates.csv"
)


f_pcoa_var_pair <- file.path(
  stat_dir,
  "Table_14_paired_CA23_nonCA23_Bray_PCoA_variance.csv"
)


f_perm_pair <- file.path(
  stat_dir,
  "Table_15_paired_CA23_nonCA23_PERMANOVA.csv"
)


required_files <- c(
  f_alpha_5g,
  f_kw_5g,
  f_pcoa_5g,
  f_pcoa_var_5g,
  f_perm_5g,
  f_alpha_pair,
  f_wil_pair,
  f_pcoa_pair,
  f_pcoa_var_pair,
  f_perm_pair
)


missing_files <- required_files[
  !file.exists(required_files)
]


if (length(missing_files) > 0) {
  
  stop(
    paste0(
      "Missing required files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


############################################################
## 2. Global style
##
## EXACTLY preserve previous submission settings
############################################################

progression_levels_lab <- c(
  "Polyp",
  "UC remission",
  "UC active",
  "Dysplasia",
  "CAC"
)


progression_colors <- c(
  "Polyp"        = "#7DA9B7",
  "UC remission" = "#8DBA91",
  "UC active"    = "#E3A04F",
  "Dysplasia"    = "#C77DA8",
  "CAC"          = "#A32635"
)


pair_levels_lab <- c(
  "nonCAC",
  "CAC"
)


pair_colors <- c(
  "nonCAC" = "#A99E79",
  "CAC"    = "#A32635"
)


base_size <- 8

font_family <- "Helvetica"


fig_theme <- theme_classic(
  base_size = base_size,
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
      size = 8.5,
      face = "plain",
      hjust = 0.5,
      margin = margin(
        b = 4
      )
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 8
    ),
    
    legend.key.size = unit(
      3.5,
      "mm"
    ),
    
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
  )


############################################################
## 3. Helper functions
##
## Preserved from previous final plotting script
############################################################

format_p <- function(
    x,
    digits = 3
) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )[1]
  
  
  if (is.na(x)) {
    
    return(
      "NA"
    )
  }
  
  
  if (x < 0.001) {
    
    return(
      formatC(
        x,
        format = "e",
        digits = 2
      )
    )
  }
  
  
  formatC(
    x,
    format = "f",
    digits = digits
  )
}


format_p_fixed <- function(
    x,
    digits = 4
) {
  
  x <- suppressWarnings(
    as.numeric(x)
  )[1]
  
  
  if (is.na(x)) {
    
    return(
      "NA"
    )
  }
  
  
  formatC(
    x,
    format = "f",
    digits = digits
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
    
    
    print(
      plot
    )
    
    
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
    
    
    print(
      plot
    )
    
    
    grDevices::dev.off()
    
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
    
    
    print(
      plot
    )
    
    
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
    
    
    print(
      plot
    )
    
    
    grDevices::dev.off()
  }
}


make_alpha_box <- function(
    dat,
    group_col,
    metric,
    colors,
    ylab,
    title_text,
    x_angle = 40,
    show_pair_lines = FALSE
) {
  
  dat2 <- dat %>%
    mutate(
      .group = .data[[group_col]],
      .value = .data[[metric]]
    )
  
  
  p <- ggplot(
    dat2,
    aes(
      x = .group,
      y = .value,
      color = .group,
      fill = .group
    )
  )
  
  
  if (
    show_pair_lines &&
    "PairID" %in% colnames(dat2)
  ) {
    
    p <- p +
      geom_line(
        aes(
          group = PairID
        ),
        color = "grey72",
        linewidth = 0.35,
        alpha = 0.55
      )
  }
  
  
  p +
    geom_boxplot(
      width = 0.58,
      linewidth = 0.55,
      outlier.shape = NA,
      alpha = 0.25,
      color = "black"
    ) +
    
    geom_jitter(
      width = 0.11,
      size = 1.8,
      alpha = 0.82,
      stroke = 0
    ) +
    
    scale_color_manual(
      values = colors,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = colors,
      drop = FALSE
    ) +
    
    labs(
      x = NULL,
      y = ylab,
      title = title_text
    ) +
    
    fig_theme +
    
    theme(
      legend.position = "none",
      
      axis.text.x = element_text(
        angle = x_angle,
        hjust = 1,
        vjust = 1,
        size = 9
      )
    )
}


make_pcoa <- function(
    dat,
    group_col,
    colors,
    xlab,
    ylab,
    title_text,
    ellipse_level = 0.95,
    point_size = 1.6,
    xlim = NULL,
    ylim = NULL
) {
  
  dat2 <- dat %>%
    mutate(
      .group = .data[[group_col]]
    )
  
  
  p <- ggplot(
    dat2,
    aes(
      x = PCoA1,
      y = PCoA2,
      color = .group
    )
  ) +
    
    stat_ellipse(
      aes(
        group = .group
      ),
      type = "norm",
      level = ellipse_level,
      linewidth = 0.55,
      alpha = 0.75,
      show.legend = FALSE
    ) +
    
    geom_point(
      size = point_size,
      alpha = 0.82
    ) +
    
    scale_color_manual(
      values = colors,
      drop = FALSE
    ) +
    
    labs(
      x = xlab,
      y = ylab,
      title = title_text
    ) +
    
    fig_theme +
    
    theme(
      legend.position = "right",
      
      axis.title.x = element_text(
        size = 8.5,
        face = "bold",
        margin = margin(
          t = 1.5
        )
      ),
      
      axis.title.y = element_text(
        size = 8.5,
        face = "bold",
        margin = margin(
          r = 3
        )
      ),
      
      plot.margin = margin(
        5,
        5,
        3,
        5
      )
    )
  
  
  if (
    !is.null(xlim) ||
    !is.null(ylim)
  ) {
    
    p <- p +
      coord_cartesian(
        xlim = xlim,
        ylim = ylim
      )
  }
  
  
  p
}


############################################################
## 4. Read and format 7-KB data
############################################################

alpha_5g <- read_csv(
  f_alpha_5g,
  show_col_types = FALSE
) %>%
  mutate(
    Group_display = recode(
      Progression5,
      "Polyp" = "Polyp",
      "UC_remission" = "UC remission",
      "UC_active" = "UC active",
      "Dysplasia" = "Dysplasia",
      "CA" = "CAC"
    ),
    
    Group_display = factor(
      Group_display,
      levels = progression_levels_lab
    )
  )


kw_5g <- read_csv(
  f_kw_5g,
  show_col_types = FALSE
)


pcoa_5g <- read_csv(
  f_pcoa_5g,
  show_col_types = FALSE
) %>%
  mutate(
    Group_display = recode(
      Progression5,
      "Polyp" = "Polyp",
      "UC_remission" = "UC remission",
      "UC_active" = "UC active",
      "Dysplasia" = "Dysplasia",
      "CA" = "CAC"
    ),
    
    Group_display = factor(
      Group_display,
      levels = progression_levels_lab
    )
  )


pcoa_var_5g <- read_csv(
  f_pcoa_var_5g,
  show_col_types = FALSE
)


perm_5g <- read_csv(
  f_perm_5g,
  show_col_types = FALSE
)


alpha_pair <- read_csv(
  f_alpha_pair,
  show_col_types = FALSE
) %>%
  mutate(
    Group_display = factor(
      Group,
      levels = pair_levels_lab
    )
  )


wil_pair <- read_csv(
  f_wil_pair,
  show_col_types = FALSE
)


pcoa_pair <- read_csv(
  f_pcoa_pair,
  show_col_types = FALSE
) %>%
  mutate(
    Group_display = factor(
      Group,
      levels = pair_levels_lab
    )
  )


pcoa_var_pair <- read_csv(
  f_pcoa_var_pair,
  show_col_types = FALSE
)


perm_pair <- read_csv(
  f_perm_pair,
  show_col_types = FALSE
)


############################################################
## 5. Statistical labels
############################################################

kw_shannon_p <- kw_5g %>%
  filter(
    Metric == "Shannon"
  ) %>%
  pull(
    p_value_raw
  )


perm_5g_r2 <- perm_5g$R2[1]

perm_5g_p <- perm_5g$p_value[1]


label_A <- if (
  kw_shannon_p < 0.001
) {
  
  "P < 0.001"
  
} else {
  
  paste0(
    "P = ",
    format_p_fixed(
      kw_shannon_p,
      3
    )
  )
}


label_B <- paste0(
  "PERMANOVA R² = ",
  formatC(
    perm_5g_r2,
    format = "f",
    digits = 3
  ),
  ", P = ",
  formatC(
    perm_5g_p,
    format = "f",
    digits = 3
  )
)


shannon_pair_p <- wil_pair %>%
  filter(
    Metric == "Shannon"
  ) %>%
  pull(
    p_value_raw
  )


label_C <- paste0(
  "P = ",
  formatC(
    shannon_pair_p,
    format = "f",
    digits = 3
  )
)


perm_pair_r2 <- perm_pair$R2[1]

perm_pair_p <- perm_pair$p_value[1]


label_D <- paste0(
  "PERMANOVA R² = ",
  formatC(
    perm_pair_r2,
    format = "f",
    digits = 3
  ),
  ", P = ",
  formatC(
    perm_pair_p,
    format = "f",
    digits = 3
  )
)


############################################################
## 6. Supplementary labels
############################################################

kw_simpson_p <- kw_5g %>%
  filter(
    Metric == "Simpson"
  ) %>%
  pull(
    p_value_raw
  )


kw_observed_p <- kw_5g %>%
  filter(
    Metric == "Observed_ASVs"
  ) %>%
  pull(
    p_value_raw
  )


wil_simpson_p <- wil_pair %>%
  filter(
    Metric == "Simpson"
  ) %>%
  pull(
    p_value_raw
  )


wil_observed_p <- wil_pair %>%
  filter(
    Metric == "Observed_ASVs"
  ) %>%
  pull(
    p_value_raw
  )


s_label_A <- if (
  kw_simpson_p < 0.001
) {
  
  "P < 0.001"
  
} else {
  
  paste0(
    "P = ",
    formatC(
      kw_simpson_p,
      format = "f",
      digits = 3
    )
  )
}


s_label_B <- paste0(
  "P = ",
  formatC(
    kw_observed_p,
    format = "f",
    digits = 3
  )
)


s_label_C <- paste0(
  "P = ",
  formatC(
    wil_simpson_p,
    format = "f",
    digits = 3
  )
)


s_label_D <- paste0(
  "P = ",
  formatC(
    wil_observed_p,
    format = "f",
    digits = 3
  )
)


############################################################
## 7. PCoA percentages
##
## Read directly from 7-KB statistics.
############################################################

pct_5g <- c(
  pcoa_var_5g$Variance_percent[
    pcoa_var_5g$Axis == "PCoA1"
  ],
  
  pcoa_var_5g$Variance_percent[
    pcoa_var_5g$Axis == "PCoA2"
  ]
)


pct_pair <- c(
  pcoa_var_pair$Variance_percent[
    pcoa_var_pair$Axis == "PCoA1"
  ],
  
  pcoa_var_pair$Variance_percent[
    pcoa_var_pair$Axis == "PCoA2"
  ]
)


pct_5g <- round(
  pct_5g,
  1
)


pct_pair <- round(
  pct_pair,
  1
)


############################################################
## 8. Figure 2A
##
## EXACT plotting settings from old submission script
############################################################

p1A <- make_alpha_box(
  dat = alpha_5g,
  group_col = "Group_display",
  metric = "Shannon",
  colors = progression_colors,
  ylab = "Shannon Index",
  title_text = label_A,
  x_angle = 40
)


############################################################
## 9. Figure 2B
##
## 95% ellipse, as in old submission script
############################################################

p1B <- make_pcoa(
  dat = pcoa_5g,
  group_col = "Group_display",
  colors = progression_colors,
  xlab = paste0(
    "Bray-Curtis PCoA 1 (",
    pct_5g[1],
    "%)"
  ),
  ylab = paste0(
    "Bray-Curtis PCoA 2 (",
    pct_5g[2],
    "%)"
  ),
  title_text = label_B,
  ellipse_level = 0.95,
  point_size = 1.6
)


############################################################
## 10. Figure 2C
############################################################

p1C <- make_alpha_box(
  dat = alpha_pair,
  group_col = "Group_display",
  metric = "Shannon",
  colors = pair_colors,
  ylab = "Shannon Index",
  title_text = label_C,
  x_angle = 40,
  show_pair_lines = FALSE
)


############################################################
## 11. Figure 2D
##
## CRITICAL:
## Old final submission uses 68% ellipse, NOT 95%.
############################################################

p1D <- make_pcoa(
  dat = pcoa_pair,
  group_col = "Group_display",
  colors = pair_colors,
  xlab = paste0(
    "Bray-Curtis PCoA 1 (",
    pct_pair[1],
    "%)"
  ),
  ylab = paste0(
    "Bray-Curtis PCoA 2 (",
    pct_pair[2],
    "%)"
  ),
  title_text = label_D,
  ellipse_level = 0.95,
  point_size = 1.8
)


############################################################
## 12. Main Figure 2
##
## EXACT layout from previous final submission script
############################################################

fig2 <- cowplot::plot_grid(
  p1A,
  p1B,
  p1C,
  p1D,
  
  ncol = 2,
  
  labels = c(
    "A",
    "B",
    "C",
    "D"
  ),
  
  label_size = 14,
  
  label_fontface = "bold",
  
  label_fontfamily = font_family,
  
  align = "none",
  
  axis = "none",
  
  rel_widths = c(
    1,
    1.18
  ),
  
  rel_heights = c(
    1,
    1
  )
)


############################################################
## 13. Supplementary Figure 1
############################################################

s1A <- make_alpha_box(
  dat = alpha_5g,
  group_col = "Group_display",
  metric = "Simpson",
  colors = progression_colors,
  ylab = "Simpson diversity",
  title_text = s_label_A,
  x_angle = 40
)


s1B <- make_alpha_box(
  dat = alpha_5g,
  group_col = "Group_display",
  metric = "Observed_ASVs",
  colors = progression_colors,
  ylab = "Observed features",
  title_text = s_label_B,
  x_angle = 40
)


s1C <- make_alpha_box(
  dat = alpha_pair,
  group_col = "Group_display",
  metric = "Simpson",
  colors = pair_colors,
  ylab = "Simpson diversity",
  title_text = s_label_C,
  x_angle = 40,
  show_pair_lines = FALSE
)


s1D <- make_alpha_box(
  dat = alpha_pair,
  group_col = "Group_display",
  metric = "Observed_ASVs",
  colors = pair_colors,
  ylab = "Observed features",
  title_text = s_label_D,
  x_angle = 40,
  show_pair_lines = FALSE
)


supp_fig1 <- cowplot::plot_grid(
  s1A,
  s1B,
  s1C,
  s1D,
  
  ncol = 2,
  
  labels = c(
    "A",
    "B",
    "C",
    "D"
  ),
  
  label_size = 14,
  
  label_fontface = "bold",
  
  label_fontfamily = font_family,
  
  align = "none",
  
  axis = "none",
  
  rel_widths = c(
    1,
    1
  ),
  
  rel_heights = c(
    1,
    1
  )
)


############################################################
## 14. Export
##
## EXACT output dimensions from old final script
############################################################

export_plot(
  fig2,
  "Figure2_7KB_submission",
  width = 7.2,
  height = 6.4
)


export_plot(
  supp_fig1,
  "Supplementary_Figure1_7KB_submission",
  width = 7.2,
  height = 6.4
)


############################################################
## 15. Individual panels
############################################################

export_plot(
  p1A,
  "Figure2A_Shannon_5group",
  width = 3.2,
  height = 3.0
)


export_plot(
  p1B,
  "Figure2B_Bray_PCoA_5group",
  width = 3.8,
  height = 3.0
)


export_plot(
  p1C,
  "Figure2C_Shannon_CA_vs_nonCAC_paired",
  width = 3.2,
  height = 3.0
)


export_plot(
  p1D,
  "Figure2D_Bray_PCoA_CA_vs_nonCAC_paired",
  width = 3.8,
  height = 3.0
)


export_plot(
  s1A,
  "Supplementary_Figure1A_Simpson_5group",
  width = 3.2,
  height = 3.0
)


export_plot(
  s1B,
  "Supplementary_Figure1B_Observed_features_5group",
  width = 3.2,
  height = 3.0
)


export_plot(
  s1C,
  "Supplementary_Figure1C_Simpson_CA_vs_nonCAC_paired",
  width = 3.2,
  height = 3.0
)


export_plot(
  s1D,
  "Supplementary_Figure1D_Observed_features_CA_vs_nonCAC_paired",
  width = 3.2,
  height = 3.0
)


############################################################
## 16. Summary
############################################################

message(
  "\nDone. Outputs written to: ",
  out_dir
)


message(
  "\nFigure 2 labels:"
)


message(
  "  A: ",
  label_A
)


message(
  "  B: ",
  label_B
)


message(
  "  C: ",
  label_C
)


message(
  "  D: ",
  label_D
)


message(
  "\nPCoA percentages:"
)


message(
  "  5-group: ",
  pct_5g[1],
  "%, ",
  pct_5g[2],
  "%"
)


message(
  "  CA vs nonCAC paired: ",
  pct_pair[1],
  "%, ",
  pct_pair[2],
  "%"
)


message(
  "\nEllipse levels:"
)


message(
  "  Figure 2B = 95%"
)


message(
  "  Figure 2D = 68%  [same as old final submission]"
)