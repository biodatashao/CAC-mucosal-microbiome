#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 05_04_assemble_Figure7.R
##
## Module 05 - Clinical association and recurrence (Figure 7, Supplementary Figure 3)
##
## Final Figure 7 — 7KB
## A B C / D E F / G H I
############################################################

options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "cowplot",
  "survival",
  "pROC"
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
      "Missing packages: ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(survival)
  library(pROC)
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

core_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "01_core_analysis"
)

dhi_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "02_DHI_original_method"
)

alpha_dir <- file.path(core_dir, "alpha_diversity")
survival_dir <- file.path(core_dir, "survival")
marker_dir <- file.path(core_dir, "marker_abundance")

beta_dir <- file.path(
  dhi_dir,
  "Figure7D_genus_Bray"
)

rmi_dir <- file.path(
  dhi_dir,
  "Figure7HI_RMI_LOOCV"
)

output_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "03_final_Figure7_exact_original_style"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Files
# ==============================================================================

ca_file <- file.path(
  alpha_dir,
  "Figure7_CAC23_Shannon_clinical_source_7KB.tsv"
)

dys_file <- file.path(
  alpha_dir,
  "Figure7_Dysplasia17_Shannon_grade_source_7KB.tsv"
)

dfs_file <- file.path(
  survival_dir,
  "Figure7E_DFS_Shannon_source_7KB.tsv"
)

marker_file <- file.path(
  marker_dir,
  "Figure7_SuppFig3_CAC23_marker_abundance_source_7KB.tsv"
)

pcoa_file <- file.path(
  beta_dir,
  "Figure7D_genus_Bray_PCoA_coordinates_7KB.tsv"
)

pcoa_stats_file <- file.path(
  beta_dir,
  "Figure7D_genus_Bray_summary_7KB.tsv"
)

rmi_file <- file.path(
  rmi_dir,
  "Figure7HI_original_LOOCV_RMI_source_7KB.tsv"
)

rmi_stats_file <- file.path(
  rmi_dir,
  "Figure7HI_original_LOOCV_RMI_statistics_7KB.tsv"
)

required_files <- c(
  ca_file,
  dys_file,
  dfs_file,
  marker_file,
  pcoa_file,
  pcoa_stats_file,
  rmi_file,
  rmi_stats_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required files:\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 3. Read
# ==============================================================================

ca_dat <- fread(
  ca_file,
  data.table = FALSE
)

dys_dat <- fread(
  dys_file,
  data.table = FALSE
)

dfs_dat <- fread(
  dfs_file,
  data.table = FALSE
)

marker_dat <- fread(
  marker_file,
  data.table = FALSE
)

pcoa_dat <- fread(
  pcoa_file,
  data.table = FALSE
)

pcoa_stats <- fread(
  pcoa_stats_file,
  data.table = FALSE
)

rmi_dat <- fread(
  rmi_file,
  data.table = FALSE
)

rmi_stats <- fread(
  rmi_stats_file,
  data.table = FALSE
)


# ==============================================================================
# 4. Common style
# ==============================================================================

font_family <- "Helvetica"

recurrence_colors <- c(
  "No recurrence" = "#A99E79",
  "Recurrence" = "#A32635"
)

stage_colors <- c(
  "Stage I" = "#7DA9B7",
  "Stage II" = "#8DBA91",
  "Stage III" = "#E3A04F"
)

grade_colors <- c(
  "High-grade dysplasia" = "#8DBA91",
  "Low-grade dysplasia" = "#C77DA8"
)

shannon_colors <- c(
  "Low Shannon" = "#7DA9B7",
  "High Shannon" = "#A32635"
)


base_theme <- theme_classic(
  base_size = 8,
  base_family = font_family
) +
  theme(
    text = element_text(
      family = font_family,
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
      hjust = 0.5,
      margin = margin(b = 5)
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


format_p <- function(p) {
  
  if (!is.finite(p)) {
    return("P = NA")
  }
  
  if (p < 0.001) {
    return("P < 0.001")
  }
  
  if (p < 0.01) {
    return(
      paste0(
        "P = ",
        formatC(
          p,
          format = "f",
          digits = 4
        )
      )
    )
  }
  
  if (p < 0.1) {
    return(
      paste0(
        "P = ",
        formatC(
          p,
          format = "f",
          digits = 3
        )
      )
    )
  }
  
  paste0(
    "P = ",
    formatC(
      p,
      format = "f",
      digits = 2
    )
  )
}


make_box <- function(
    dat,
    x_col,
    y_col,
    colors,
    y_lab,
    title_text,
    test = "wilcox"
) {
  
  dat$.x <- factor(
    dat[, x_col],
    levels = names(colors)
  )
  
  dat$.y <- as.numeric(
    dat[, y_col]
  )
  
  keep <- !is.na(dat$.x) &
    is.finite(dat$.y)
  
  dat <- dat[
    keep,
    ,
    drop = FALSE
  ]
  
  if (test == "kruskal") {
    
    p_value <- kruskal.test(
      .y ~ .x,
      data = dat
    )$p.value
    
  } else {
    
    p_value <- suppressWarnings(
      wilcox.test(
        .y ~ .x,
        data = dat,
        exact = FALSE
      )$p.value
    )
  }
  
  ggplot(
    dat,
    aes(
      x = .x,
      y = .y,
      color = .x
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
      width = 0.11,
      height = 0,
      size = 1.8,
      alpha = 0.82,
      stroke = 0
    ) +
    scale_color_manual(
      values = colors,
      drop = FALSE
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = format_p(p_value),
      hjust = 1.05,
      vjust = 1.35,
      size = 3,
      fontface = "bold",
      family = font_family
    ) +
    labs(
      x = NULL,
      y = y_lab,
      title = title_text
    ) +
    base_theme +
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 35,
        hjust = 1,
        vjust = 1
      )
    )
}


get_beta <- function(metric_name) {
  
  value <- pcoa_stats$Value[
    pcoa_stats$Metric == metric_name
  ]
  
  if (length(value) != 1) {
    stop(
      paste0(
        "Cannot extract beta metric: ",
        metric_name
      ),
      call. = FALSE
    )
  }
  
  as.numeric(value)
}


get_rmi <- function(metric_name) {
  
  value <- rmi_stats$Value[
    rmi_stats$Metric == metric_name
  ]
  
  if (length(value) != 1) {
    stop(
      paste0(
        "Cannot extract RMI metric: ",
        metric_name
      ),
      call. = FALSE
    )
  }
  
  as.numeric(value)
}


# ==============================================================================
# 5. A
# ==============================================================================

ca_dat$Tumor_stage_plot <- factor(
  ca_dat$Tumor_stage_plot,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III"
  )
)

pA <- make_box(
  dat = ca_dat,
  x_col = "Tumor_stage_plot",
  y_col = "Shannon",
  colors = stage_colors,
  y_lab = "Shannon diversity",
  title_text = "Tumor Shannon diversity by tumor stage",
  test = "kruskal"
)


# ==============================================================================
# 6. B
# ==============================================================================

dys_dat$Dysplasia_grade_plot <- factor(
  dys_dat$Dysplasia_grade_plot,
  levels = c(
    "High-grade dysplasia",
    "Low-grade dysplasia"
  )
)

pB <- make_box(
  dat = dys_dat,
  x_col = "Dysplasia_grade_plot",
  y_col = "Shannon",
  colors = grade_colors,
  y_lab = "Shannon diversity",
  title_text = "Dysplasia Shannon diversity by grade",
  test = "wilcox"
)


# ==============================================================================
# 7. C
# ==============================================================================

ca_dat$Recurrence_plot <- factor(
  ca_dat$Recurrence_plot,
  levels = c(
    "No recurrence",
    "Recurrence"
  )
)

pC <- make_box(
  dat = ca_dat,
  x_col = "Recurrence_plot",
  y_col = "Shannon",
  colors = recurrence_colors,
  y_lab = "Shannon diversity",
  title_text = "Tumor Shannon diversity by recurrence",
  test = "wilcox"
)


# ==============================================================================
# 8. D — genus-level Bray-Curtis PCoA
# ==============================================================================

pcoa_dat$Recurrence_plot <- factor(
  pcoa_dat$Recurrence_plot,
  levels = c(
    "No recurrence",
    "Recurrence"
  )
)

pc1 <- get_beta(
  "PCoA1_variance_percent"
)

pc2 <- get_beta(
  "PCoA2_variance_percent"
)

perm_r2 <- get_beta(
  "PERMANOVA_R2_recurrence"
)

perm_p <- get_beta(
  "PERMANOVA_P_recurrence"
)


pD <- ggplot(
  pcoa_dat,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Recurrence_plot
  )
) +
  stat_ellipse(
    aes(
      group = Recurrence_plot
    ),
    type = "norm",
    level = 0.95,
    linewidth = 0.55,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  geom_point(
    size = 1.8,
    alpha = 0.82
  ) +
  scale_color_manual(
    values = recurrence_colors,
    breaks = c(
      "No recurrence",
      "Recurrence"
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(
      -0.4,
      0.0,
      0.4,
      0.8
    )
  ) +
  scale_y_continuous(
    breaks = c(
      -0.4,
      0.0,
      0.4
    )
  ) +
  labs(
    title = paste0(
      "PERMANOVA R² = ",
      formatC(
        perm_r2,
        format = "f",
        digits = 3
      ),
      ", P = ",
      formatC(
        perm_p,
        format = "f",
        digits = 3
      )
    ),
    x = paste0(
      "Bray-Curtis PCoA 1 (",
      formatC(
        pc1,
        format = "f",
        digits = 1
      ),
      "%)"
    ),
    y = paste0(
      "Bray-Curtis PCoA 2 (",
      formatC(
        pc2,
        format = "f",
        digits = 1
      ),
      "%)"
    ),
    color = NULL
  ) +
  base_theme +
  theme(
    plot.title = element_text(
      size = 9,
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right",
    legend.title = element_blank()
  )


# ==============================================================================
# 9. E — DFS
# ==============================================================================

dfs_dat$Shannon_group <- factor(
  dfs_dat$Shannon_group,
  levels = c(
    "Low Shannon",
    "High Shannon"
  )
)

dfs_dat$DFS_time_months <- as.numeric(
  dfs_dat$DFS_time_months
)

dfs_dat$DFS_event <- as.numeric(
  dfs_dat$DFS_event
)


km_fit <- survfit(
  Surv(
    DFS_time_months,
    DFS_event
  ) ~ Shannon_group,
  data = dfs_dat
)

km_diff <- survdiff(
  Surv(
    DFS_time_months,
    DFS_event
  ) ~ Shannon_group,
  data = dfs_dat
)

logrank_p <- pchisq(
  km_diff$chisq,
  df = length(km_diff$n) - 1,
  lower.tail = FALSE
)


km_sum <- summary(
  km_fit
)

km_dat <- data.frame(
  time = km_sum$time,
  surv = km_sum$surv,
  n.censor = km_sum$n.censor,
  strata = as.character(
    km_sum$strata
  ),
  stringsAsFactors = FALSE
)

km_dat$Shannon_group <- sub(
  "^Shannon_group=",
  "",
  km_dat$strata
)

start_dat <- data.frame(
  time = c(0, 0),
  surv = c(1, 1),
  n.censor = c(0, 0),
  strata = c(
    "Shannon_group=Low Shannon",
    "Shannon_group=High Shannon"
  ),
  Shannon_group = c(
    "Low Shannon",
    "High Shannon"
  ),
  stringsAsFactors = FALSE
)

km_dat <- rbind(
  start_dat,
  km_dat
)

km_dat$Shannon_group <- factor(
  km_dat$Shannon_group,
  levels = c(
    "Low Shannon",
    "High Shannon"
  )
)

km_dat <- km_dat[
  order(
    km_dat$Shannon_group,
    km_dat$time
  ),
  ,
  drop = FALSE
]

censor_dat <- km_dat[
  km_dat$n.censor > 0,
  ,
  drop = FALSE
]


pE <- ggplot(
  km_dat,
  aes(
    x = time,
    y = surv,
    color = Shannon_group,
    group = Shannon_group
  )
) +
  geom_step(
    linewidth = 0.75
  ) +
  geom_point(
    data = censor_dat,
    shape = 3,
    size = 1.8,
    stroke = 0.7
  ) +
  scale_color_manual(
    values = shannon_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1.02
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "Log-rank ",
      format_p(
        logrank_p
      )
    ),
    hjust = 1.05,
    vjust = 1.35,
    size = 3,
    fontface = "bold",
    family = font_family
  ) +
  labs(
    x = "DFS time (months)",
    y = "Disease-free survival",
    title = "DFS according to tumor Shannon diversity",
    color = "Shannon_group"
  ) +
  base_theme +
  theme(
    legend.position = c(
      0.70,
      0.18
    ),
    legend.background = element_blank(),
    legend.key = element_blank()
  )


# ==============================================================================
# 10. F — Lactococcus
# ==============================================================================

marker_dat$Recurrence_plot <- factor(
  marker_dat$Recurrence_plot,
  levels = c(
    "No recurrence",
    "Recurrence"
  )
)

lacto_dat <- marker_dat[
  marker_dat$Genus == "Lactococcus",
  ,
  drop = FALSE
]

pF <- make_box(
  dat = lacto_dat,
  x_col = "Recurrence_plot",
  y_col = "Relative_abundance_percent",
  colors = recurrence_colors,
  y_lab = "Relative abundance (%)",
  title_text = expression(
    italic(Lactococcus)~"abundance according to recurrence"
  ),
  test = "wilcox"
)


# ==============================================================================
# 11. G — UCG-005
# ==============================================================================

ucg_dat <- marker_dat[
  marker_dat$Genus == "UCG_005",
  ,
  drop = FALSE
]

pG <- make_box(
  dat = ucg_dat,
  x_col = "Recurrence_plot",
  y_col = "Relative_abundance_percent",
  colors = recurrence_colors,
  y_lab = "Relative abundance (%)",
  title_text = "UCG-005 abundance according to recurrence",
  test = "wilcox"
)


# ==============================================================================
# 12. H — original strict LOOCV RMI
# ==============================================================================

rmi_dat$Recurrence_plot <- factor(
  rmi_dat$Recurrence_plot,
  levels = c(
    "No recurrence",
    "Recurrence"
  )
)

rmi_p <- get_rmi(
  "RMI_Wilcoxon_P"
)


pH <- ggplot(
  rmi_dat,
  aes(
    x = Recurrence_plot,
    y = RMI_two_genus,
    color = Recurrence_plot
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
    width = 0.11,
    height = 0,
    size = 1.8,
    alpha = 0.82,
    stroke = 0
  ) +
  scale_color_manual(
    values = recurrence_colors,
    drop = FALSE
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = format_p(
      rmi_p
    ),
    hjust = 1.05,
    vjust = 1.35,
    size = 3,
    fontface = "bold",
    family = font_family
  ) +
  labs(
    x = NULL,
    y = "Two-genus microbial index"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
    )
  )


# ==============================================================================
# 13. I — ROC
#
# IMPORTANT:
# Directly reconstruct from the 23 out-of-fold RMI scores.
# Do NOT use previously exported ROC-coordinate table.
# ==============================================================================

roc_object <- pROC::roc(
  response = rmi_dat$Recurrence_binary,
  predictor = rmi_dat$RMI_two_genus,
  levels = c(
    0,
    1
  ),
  direction = "<",
  quiet = TRUE
)


auc_value <- get_rmi(
  "LOOCV_AUC"
)

auc_lower <- get_rmi(
  "LOOCV_AUC_CI_lower"
)

auc_upper <- get_rmi(
  "LOOCV_AUC_CI_upper"
)


roc_dat <- data.frame(
  FPR = 1 - as.numeric(
    roc_object$specificities
  ),
  Sensitivity = as.numeric(
    roc_object$sensitivities
  ),
  stringsAsFactors = FALSE
)


# Order the staircase strictly from left to right.
roc_dat <- roc_dat[
  order(
    roc_dat$FPR,
    roc_dat$Sensitivity
  ),
  ,
  drop = FALSE
]


# For identical FPR values retain the highest achieved sensitivity.
roc_dat <- aggregate(
  Sensitivity ~ FPR,
  data = roc_dat,
  FUN = max
)


roc_dat <- roc_dat[
  order(
    roc_dat$FPR
  ),
  ,
  drop = FALSE
]


# Enforce monotonic sensitivity.
roc_dat$Sensitivity <- cummax(
  roc_dat$Sensitivity
)


# Ensure standard ROC endpoints are present.
roc_dat <- rbind(
  data.frame(
    FPR = 0,
    Sensitivity = 0
  ),
  roc_dat,
  data.frame(
    FPR = 1,
    Sensitivity = 1
  )
)


roc_dat <- aggregate(
  Sensitivity ~ FPR,
  data = roc_dat,
  FUN = max
)

roc_dat <- roc_dat[
  order(
    roc_dat$FPR
  ),
  ,
  drop = FALSE
]

roc_dat$Sensitivity <- cummax(
  roc_dat$Sensitivity
)


roc_label <- paste0(
  "LOOCV AUC = ",
  formatC(
    auc_value,
    format = "f",
    digits = 3
  ),
  "\n95% CI ",
  formatC(
    auc_lower,
    format = "f",
    digits = 3
  ),
  "-",
  formatC(
    auc_upper,
    format = "f",
    digits = 3
  )
)


pI <- ggplot(
  roc_dat,
  aes(
    x = FPR,
    y = Sensitivity
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey60"
  ) +
  geom_step(
    direction = "hv",
    linewidth = 0.8,
    color = "#A32635"
  ) +
  annotate(
    "text",
    x = 0.97,
    y = 0.05,
    label = roc_label,
    hjust = 1,
    vjust = 0,
    size = 2.5,
    family = font_family,
    color = "black"
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand = c(
      0,
      0
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand = c(
      0,
      0
    )
  ) +
  coord_equal() +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  base_theme +
  theme(
    legend.position = "none"
  )


# ==============================================================================
# 14. Final assembly
# ==============================================================================

shared_theme <- theme(
  plot.title = element_text(
    family = font_family,
    size = 9,
    face = "bold",
    hjust = 0.5,
    margin = margin(
      b = 5,
      unit = "pt"
    )
  ),
  axis.title = element_text(
    family = font_family,
    size = 8.5,
    face = "bold"
  ),
  axis.text = element_text(
    family = font_family,
    size = 7.5,
    color = "black"
  ),
  legend.title = element_text(
    family = font_family,
    size = 7.5
  ),
  legend.text = element_text(
    family = font_family,
    size = 7
  ),
  plot.margin = margin(
    t = 7,
    r = 7,
    b = 7,
    l = 7,
    unit = "pt"
  )
)


pA <- pA + shared_theme
pB <- pB + shared_theme
pC <- pC + shared_theme
pD <- pD + shared_theme
pE <- pE + shared_theme
pF <- pF + shared_theme
pG <- pG + shared_theme
pH <- pH + shared_theme
pI <- pI + shared_theme


aligned_panels <- cowplot::align_plots(
  pA,
  pB,
  pC,
  pD,
  pE,
  pF,
  pG,
  pH,
  pI,
  align = "hv",
  axis = "tblr"
)


figure7_final <- cowplot::plot_grid(
  plotlist = aligned_panels,
  ncol = 3,
  nrow = 3,
  labels = LETTERS[1:9],
  label_size = 14,
  label_fontface = "bold",
  label_fontfamily = font_family,
  label_x = 0.005,
  label_y = 0.995,
  hjust = 0,
  vjust = 1,
  align = "hv",
  axis = "tblr"
)


# ==============================================================================
# 15. Export
# ==============================================================================

output_stem <- file.path(
  output_dir,
  "Figure7_clinical_recurrence_9panels_7KB_final"
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

rds_file <- paste0(
  output_stem,
  ".rds"
)


figure_width <- 11.5
figure_height <- 10.5
figure_dpi <- 600


grDevices::pdf(
  pdf_file,
  width = figure_width,
  height = figure_height,
  family = font_family,
  useDingbats = FALSE
)

print(
  figure7_final
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
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    background = "white"
  )
  
} else {
  
  grDevices::png(
    png_file,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    bg = "white",
    type = "quartz"
  )
}

print(
  figure7_final
)

grDevices::dev.off()


if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_tiff(
    tiff_file,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    compression = "lzw",
    background = "white"
  )
  
} else {
  
  grDevices::tiff(
    tiff_file,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    compression = "lzw",
    bg = "white",
    type = "quartz"
  )
}

print(
  figure7_final
)

grDevices::dev.off()


saveRDS(
  figure7_final,
  rds_file
)


# ==============================================================================
# 16. Summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("FINAL FIGURE 7 COMPLETED\n")
cat("============================================================\n")

cat(
  "D PERMANOVA R2 = ",
  perm_r2,
  ", P = ",
  perm_p,
  "\n",
  sep = ""
)

cat(
  "D PCoA1 = ",
  pc1,
  "%; PCoA2 = ",
  pc2,
  "%\n",
  sep = ""
)

cat(
  "H RMI P = ",
  rmi_p,
  "\n",
  sep = ""
)

cat(
  "I AUC = ",
  auc_value,
  " (95% CI ",
  auc_lower,
  "-",
  auc_upper,
  ")\n",
  sep = ""
)

cat("\nPDF:\n")
cat(pdf_file, "\n")

cat("\nPNG:\n")
cat(png_file, "\n")

cat("\nTIFF:\n")
cat(tiff_file, "\n")