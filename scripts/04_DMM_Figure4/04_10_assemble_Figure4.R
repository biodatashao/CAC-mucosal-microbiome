#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_10_assemble_Figure4.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Final 7KB Figure 4 — ALL IN ONE
##
## Reproduces the original final Figure 4 plotting workflow using:
##
## A. DMM Laplace model selection
## B. DMM-state distribution across progression groups
## C. Lavelle-style cluster-mean heatmap
## D. Shannon / observed genera / dominant genus proportion
## E. Paired nonCAC -> CAC DMM-state correspondence
##
## Final layout:
##
## A | B | E
## C | C | D
##
## IMPORTANT
############################################################
# - Uses NEW 7KB results only.
# - Does NOT use the previously generated temporary panel RDS files.
# - Restores the original manuscript colour palette and plotting style.
# - Panel C uses the original final 29-genus display set:
#     C1 = 6
#     C2 = 13
#     C3 = 10
# - Panel D is VERTICAL, as in the original final manuscript.
# - Panel E uses discern = TRUE, as in the original final script.
# - No Cairo / X11 device.
# ==============================================================================


options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "scales",
  "dplyr",
  "tidyr",
  "ggalluvial"
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
  library(patchwork)
  library(scales)
  library(dplyr)
  library(tidyr)
  library(ggalluvial)
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

dmm_root <- file.path(
  rerun_root,
  "03_Figure4_DMM"
)

primary_dir <- file.path(
  dmm_root,
  "DMM_progression127_primary_prevalence10_7KB"
)

characterization_dir <- file.path(
  dmm_root,
  "DMM_progression127_K3_characterization_7KB"
)

heatmap_source_dir <- file.path(
  dmm_root,
  "Figure4C_Lavelle_cluster_mean_heatmap_7KB",
  "source_tables"
)

projection_dir <- file.path(
  dmm_root,
  "projection_CA23_nonCA23_fixed_K3_7KB",
  "tables"
)

output_dir <- file.path(
  dmm_root,
  "Figure4_final_all_in_one_7KB"
)

source_dir <- file.path(
  output_dir,
  "source_tables"
)

panel_dir <- file.path(
  output_dir,
  "panels"
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

dir.create(
  panel_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 2. Input files
# ==============================================================================

assignment_file <- file.path(
  primary_dir,
  "tables",
  "DMM_sample_assignments_optimalK3_native_7KB.tsv"
)

metric_file <- file.path(
  characterization_dir,
  "tables",
  "DMM_K3_genus_diversity_and_dominance_7KB.tsv"
)

heatmap_all_file <- file.path(
  heatmap_source_dir,
  "all_genera_cluster_summary_long_7KB.tsv"
)

heatmap_selected_file <- file.path(
  heatmap_source_dir,
  "Figure4C_selected_qualified_genera_7KB.tsv"
)

paired_transition_file <- file.path(
  projection_dir,
  "DMM_fixed_K3_paired_state_transition_source_7KB.tsv"
)


required_files <- c(
  assignment_file,
  metric_file,
  heatmap_all_file,
  paired_transition_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing required file(s):\n",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 3. Style — EXACT ORIGINAL FINAL FIGURE 4 PALETTE
# ==============================================================================

base_family <- "Helvetica"
base_size <- 8

cluster_levels <- c(
  "C1",
  "C2",
  "C3"
)

pathology_levels <- c(
  "Polyp",
  "UC remission",
  "UC active",
  "Dysplasia",
  "CAC"
)

cluster_colors <- c(
  C1 = "#78A89C",
  C2 = "#B74F5B",
  C3 = "#D5A24C"
)

heatmap_low <- "#5579B6"
heatmap_mid <- "#F7F7F7"
heatmap_high <- "#C95E57"


theme_manuscript <- ggplot2::theme_classic(
  base_size = base_size,
  base_family = base_family
) +
  ggplot2::theme(
    text = ggplot2::element_text(
      family = base_family,
      color = "black"
    ),
    
    axis.title = ggplot2::element_text(
      size = 9,
      face = "bold",
      color = "black"
    ),
    
    axis.text = ggplot2::element_text(
      size = 8,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    axis.ticks = ggplot2::element_line(
      linewidth = 0.40,
      color = "black"
    ),
    
    legend.title = ggplot2::element_text(
      size = 8.5,
      face = "bold"
    ),
    
    legend.text = ggplot2::element_text(
      size = 7.5
    ),
    
    plot.title = ggplot2::element_text(
      size = 9,
      face = "bold",
      hjust = 0
    ),
    
    plot.subtitle = ggplot2::element_text(
      size = 7.5,
      color = "black"
    )
  )


# ==============================================================================
# 4. Helpers
# ==============================================================================

section <- function(text) {
  
  cat("\n")
  cat(
    paste(
      rep("=", 80),
      collapse = ""
    ),
    "\n"
  )
  
  cat(
    text,
    "\n"
  )
  
  cat(
    paste(
      rep("=", 80),
      collapse = ""
    ),
    "\n"
  )
}


format_p <- function(p) {
  
  if (is.na(p)) {
    return("P = NA")
  }
  
  if (p < 0.001) {
    return("P < 0.001")
  }
  
  paste0(
    "P = ",
    formatC(
      p,
      format = "f",
      digits = 3
    )
  )
}


find_column <- function(data, patterns, description) {
  
  column_names <- colnames(data)
  
  for (pattern in patterns) {
    
    hit <- column_names[
      grepl(
        pattern,
        column_names,
        ignore.case = TRUE
      )
    ]
    
    if (length(hit) > 0) {
      return(hit[1])
    }
  }
  
  stop(
    paste0(
      "Could not identify column: ",
      description,
      "\nColumns available:\n",
      paste(
        column_names,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


clean_cluster <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  out <- rep(
    NA_character_,
    length(x)
  )
  
  out[
    grepl(
      "C1|cluster.?1|component.?1|state.?1|^1$",
      x,
      ignore.case = TRUE
    )
  ] <- "C1"
  
  out[
    grepl(
      "C2|cluster.?2|component.?2|state.?2|^2$",
      x,
      ignore.case = TRUE
    )
  ] <- "C2"
  
  out[
    grepl(
      "C3|cluster.?3|component.?3|state.?3|^3$",
      x,
      ignore.case = TRUE
    )
  ] <- "C3"
  
  out
}


row_zscore <- function(x) {
  
  x <- as.numeric(x)
  
  if (
    length(x) != 3 ||
    any(!is.finite(x))
  ) {
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }
  
  current_sd <- stats::sd(x)
  
  if (
    is.na(current_sd) ||
    current_sd == 0
  ) {
    return(
      rep(
        0,
        length(x)
      )
    )
  }
  
  z <- (
    x -
      mean(x)
  ) /
    current_sd
  
  z[
    z > 1.5
  ] <- 1.5
  
  z[
    z < -1.5
  ] <- -1.5
  
  z
}


# ==============================================================================
# 5. PANEL A — Laplace model selection
#
# Locked values from the completed 7KB primary DMM rerun.
# Plotting parameters are locked to the submitted Figure 4.
# ==============================================================================

section(
  "PANEL A: LAPLACE MODEL SELECTION"
)

laplace_table <- data.frame(
  K = 1:7,
  
  Laplace = c(
    90707,
    87151,
    86176,
    86185,
    86325,
    86635,
    87044
  ),
  
  Selected = 1:7 == 3,
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  laplace_table,
  file.path(
    source_dir,
    "Figure4A_Laplace_values_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


panel_a <- ggplot2::ggplot(
  laplace_table,
  ggplot2::aes(
    x = K,
    y = Laplace
  )
) +
  
  ggplot2::geom_line(
    linewidth = 0.65,
    color = "#4A4A4A"
  ) +
  
  ggplot2::geom_point(
    data = laplace_table[
      !laplace_table$Selected,
      ,
      drop = FALSE
    ],
    
    shape = 21,
    size = 2.5,
    stroke = 0.45,
    fill = "white",
    color = "black"
  ) +
  
  ggplot2::geom_point(
    data = laplace_table[
      laplace_table$Selected,
      ,
      drop = FALSE
    ],
    
    shape = 21,
    size = 2.8,
    stroke = 0.45,
    fill = cluster_colors["C2"],
    color = "black"
  ) +
  
  ggplot2::annotate(
    geom = "text",
    x = 3,
    y = laplace_table$Laplace[3],
    label = "K = 3",
    hjust = 0.5,
    vjust = -1.15,
    family = base_family,
    fontface = "bold",
    size = 2.7
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = 1:7
  ) +
  
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 1,
      big.mark = ","
    ),
    
    expand = ggplot2::expansion(
      mult = c(
        0.08,
        0.08
      )
    )
  ) +
  
  ggplot2::labs(
    x = "Number of DMM components (K)",
    y = "Laplace approximation"
  ) +
  
  theme_manuscript


# ==============================================================================
# 6. PANEL B — progression-group DMM state distribution
# ==============================================================================

section(
  "PANEL B: DMM STATE DISTRIBUTION"
)

assignment <- data.table::fread(
  assignment_file,
  data.table = FALSE,
  check.names = FALSE
)


required_assignment_columns <- c(
  "SampleID",
  "DMM_cluster",
  "Progression5"
)

missing_assignment_columns <- setdiff(
  required_assignment_columns,
  colnames(assignment)
)

if (length(missing_assignment_columns) > 0) {
  stop(
    paste0(
      "Assignment table missing: ",
      paste(
        missing_assignment_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


assignment$Cluster <- clean_cluster(
  assignment$DMM_cluster
)

assignment$Pathology <- as.character(
  assignment$Progression5
)

assignment$Pathology[
  assignment$Pathology == "UC_remission"
] <- "UC remission"

assignment$Pathology[
  assignment$Pathology == "UC_active"
] <- "UC active"

assignment$Pathology[
  assignment$Pathology %in%
    c(
      "CA",
      "CAC"
    )
] <- "CAC"


assignment$Pathology <- factor(
  assignment$Pathology,
  levels = pathology_levels
)

assignment$Cluster <- factor(
  assignment$Cluster,
  levels = cluster_levels
)


if (
  anyNA(
    assignment$Pathology
  ) ||
  anyNA(
    assignment$Cluster
  )
) {
  stop(
    "Unexpected progression or DMM-state label in assignment table.",
    call. = FALSE
  )
}


distribution_table <- as.data.frame(
  table(
    assignment$Pathology,
    assignment$Cluster
  ),
  stringsAsFactors = FALSE
)

colnames(
  distribution_table
) <- c(
  "Pathology",
  "Cluster",
  "Count"
)

distribution_table$Total <- ave(
  distribution_table$Count,
  distribution_table$Pathology,
  FUN = sum
)

distribution_table$Proportion <-
  distribution_table$Count /
  distribution_table$Total


distribution_table$Pathology <- factor(
  distribution_table$Pathology,
  levels = pathology_levels
)

distribution_table$Cluster <- factor(
  distribution_table$Cluster,
  levels = cluster_levels
)


set.seed(
  20260715
)

fisher_result <- stats::fisher.test(
  table(
    assignment$Pathology,
    assignment$Cluster
  ),
  simulate.p.value = TRUE,
  B = 100000
)

fisher_p <- fisher_result$p.value


cac_c2 <- distribution_table[
  distribution_table$Pathology == "CAC" &
    distribution_table$Cluster == "C2",
  ,
  drop = FALSE
]


cac_c2_label <- paste0(
  cac_c2$Count,
  "/",
  cac_c2$Total,
  "\n",
  formatC(
    100 *
      cac_c2$Proportion,
    format = "f",
    digits = 1
  ),
  "%"
)


panel_b <- ggplot2::ggplot(
  distribution_table,
  ggplot2::aes(
    x = Pathology,
    y = Proportion,
    fill = Cluster
  )
) +
  
  ggplot2::geom_col(
    width = 0.70,
    color = "white",
    linewidth = 0.25
  ) +
  
  ggplot2::annotate(
    geom = "text",
    x = 5,
    y = 0.50,
    label = cac_c2_label,
    color = "white",
    family = base_family,
    fontface = "bold",
    lineheight = 0.9,
    size = 2.6
  ) +
  
  ggplot2::scale_fill_manual(
    values = cluster_colors,
    breaks = cluster_levels,
    drop = FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      1
    ),
    
    breaks = seq(
      0,
      1,
      0.25
    ),
    
    labels = scales::label_percent(
      accuracy = 1
    ),
    
    expand = c(
      0,
      0
    )
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Proportion of samples",
    fill = "DMM state"
  ) +
  
  theme_manuscript +
  
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 7.5
    ),
    
    legend.position = "bottom"
  ) +
  
  ggplot2::guides(
    fill = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE,
      title.position = "left"
    )
  )


data.table::fwrite(
  distribution_table,
  file.path(
    source_dir,
    "Figure4B_cluster_distribution_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 7. PANEL C — final 29-genus heatmap
#
# Original final display set:
#
# C1 = fixed 6 representative genera
# C2 = 13 genera
# C3 = 10 genera
#
# Here all 29 genera are extracted directly from the NEW 7KB all-genus
# cluster-mean table. This prevents the two original representative C1 genera
# from disappearing merely because they missed the upstream qualification
# threshold in this rerun.
# ==============================================================================

section(
  "PANEL C: FINAL 29-GENUS HEATMAP"
)


heatmap_all <- data.table::fread(
  heatmap_all_file,
  data.table = FALSE,
  check.names = FALSE
)


genus_col <- find_column(
  heatmap_all,
  c(
    "^Genus$",
    "genus"
  ),
  "genus"
)


cluster_col <- find_column(
  heatmap_all,
  c(
    "^DMM_cluster$",
    "^Cluster$",
    "assigned_cluster",
    "cluster"
  ),
  "DMM cluster"
)


mean_col <- find_column(
  heatmap_all,
  c(
    "^Mean_relative_abundance$",
    "mean.*relative.*abundance",
    "cluster.*mean",
    "^Mean_RA$",
    "^Mean$"
  ),
  "cluster mean relative abundance"
)


heatmap_all_simple <- data.frame(
  Genus = as.character(
    heatmap_all[, genus_col]
  ),
  
  Cluster = clean_cluster(
    heatmap_all[, cluster_col]
  ),
  
  Mean_abundance = as.numeric(
    heatmap_all[, mean_col]
  ),
  
  stringsAsFactors = FALSE
)


heatmap_all_simple <- heatmap_all_simple[
  !is.na(
    heatmap_all_simple$Genus
  ) &
    heatmap_all_simple$Genus != "" &
    heatmap_all_simple$Cluster %in%
    cluster_levels,
  ,
  drop = FALSE
]


# ------------------------------------------------------------------------------
# Exact final genus order from the original manuscript
# ------------------------------------------------------------------------------

c1_genera <- c(
  "Enterococcus",
  "Neobacillus",
  "Clostridium",
  "Streptococcus",
  "Actinomyces",
  "Lachnospiraceae_NK3A20_group"
)


c2_genera <- c(
  "UCG_005",
  "Extibacter",
  "Mesobacillus",
  "Anaerococcus",
  "Romboutsia",
  "UCG_009",
  "Blautia",
  "Caldimonas",
  "Thomasclavelia",
  "Butyribacter",
  "Peptoniphilus",
  "Faecalibaculum",
  "Lachnospiraceae_UCG_010"
)


c3_genera <- c(
  "Lactobacillus",
  "Christensenellaceae_R_7_group",
  "Lachnospiraceae_NK4A136_group",
  "Akkermansia",
  "Listeria",
  "Fannyhessea",
  "Adlercreutzia",
  "Weissella",
  "Dubosiella",
  "Prevotellaceae_UCG_001"
)


# ------------------------------------------------------------------------------
# Name alias handling
#
# Original final plot uses hyphenated display names while the 7KB analysis
# source may use underscore taxonomy labels.
# ------------------------------------------------------------------------------

canonicalize_genus <- function(x) {
  
  x <- as.character(x)
  
  x[
    x == "UCG-005"
  ] <- "UCG_005"
  
  x[
    x == "UCG-009"
  ] <- "UCG_009"
  
  x[
    x == "Christensenellaceae_R-7_group"
  ] <- "Christensenellaceae_R_7_group"
  
  x[
    x == "Lachnospiraceae_UCG-010"
  ] <- "Lachnospiraceae_UCG_010"
  
  x[
    x == "Prevotellaceae_UCG-001"
  ] <- "Prevotellaceae_UCG_001"
  
  x
}


heatmap_all_simple$Genus <- canonicalize_genus(
  heatmap_all_simple$Genus
)


final_genus_order <- c(
  c1_genera,
  c2_genera,
  c3_genera
)


missing_final_genera <- setdiff(
  final_genus_order,
  unique(
    heatmap_all_simple$Genus
  )
)


if (length(missing_final_genera) > 0) {
  stop(
    paste0(
      "The following final Figure 4C genera are absent from the NEW 7KB ",
      "all-genus cluster-mean source:\n",
      paste(
        missing_final_genera,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# Build complete 29 x 3 matrix
# ------------------------------------------------------------------------------

heatmap_selected <- heatmap_all_simple[
  heatmap_all_simple$Genus %in%
    final_genus_order,
  ,
  drop = FALSE
]


heatmap_selected <- stats::aggregate(
  Mean_abundance ~ Genus + Cluster,
  data = heatmap_selected,
  FUN = mean
)


heatmap_complete <- tidyr::complete(
  heatmap_selected,
  Genus = final_genus_order,
  Cluster = cluster_levels
)


if (
  anyNA(
    heatmap_complete$Mean_abundance
  )
) {
  
  bad_rows <- heatmap_complete[
    is.na(
      heatmap_complete$Mean_abundance
    ),
    ,
    drop = FALSE
  ]
  
  print(
    bad_rows
  )
  
  stop(
    "Figure 4C 29 x 3 abundance matrix is incomplete.",
    call. = FALSE
  )
}


# ------------------------------------------------------------------------------
# Row-scale each genus exactly across C1/C2/C3
# ------------------------------------------------------------------------------

heatmap_wide <- tidyr::pivot_wider(
  heatmap_complete,
  names_from = Cluster,
  values_from = Mean_abundance
)


heatmap_z <- heatmap_wide


for (i in seq_len(nrow(heatmap_z))) {
  
  current_values <- c(
    heatmap_wide$C1[i],
    heatmap_wide$C2[i],
    heatmap_wide$C3[i]
  )
  
  current_z <- row_zscore(
    current_values
  )
  
  heatmap_z$C1[i] <- current_z[1]
  heatmap_z$C2[i] <- current_z[2]
  heatmap_z$C3[i] <- current_z[3]
}


heatmap_long <- tidyr::pivot_longer(
  heatmap_z,
  cols = c(
    "C1",
    "C2",
    "C3"
  ),
  names_to = "Cluster",
  values_to = "Value"
)


heatmap_long$Cluster <- factor(
  heatmap_long$Cluster,
  levels = cluster_levels
)


heatmap_long$Genus <- factor(
  heatmap_long$Genus,
  levels = rev(
    final_genus_order
  )
)


# Original boundaries:
# 6 C1 + 13 C2 + 10 C3

separator_y <- c(
  length(final_genus_order) - 6 + 0.5,
  length(final_genus_order) - 19 + 0.5
)


# ------------------------------------------------------------------------------
# Display taxonomy labels exactly as original final figure
# ------------------------------------------------------------------------------

display_name_map <- c(
  "UCG_005" = "UCG-005",
  "UCG_009" = "UCG-009",
  "Christensenellaceae_R_7_group" =
    "Christensenellaceae_R-7_group",
  "Lachnospiraceae_UCG_010" =
    "Lachnospiraceae_UCG-010",
  "Prevotellaceae_UCG_001" =
    "Prevotellaceae_UCG-001"
)


display_genus <- function(x) {
  
  x <- as.character(x)
  
  replace_idx <- x %in%
    names(
      display_name_map
    )
  
  x[
    replace_idx
  ] <- unname(
    display_name_map[
      x[
        replace_idx
      ]
    ]
  )
  
  x
}


panel_c <- ggplot2::ggplot(
  heatmap_long,
  ggplot2::aes(
    x = Cluster,
    y = Genus,
    fill = Value
  )
) +
  
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.30
  ) +
  
  ggplot2::geom_hline(
    yintercept = separator_y,
    linewidth = 0.55,
    color = "black"
  ) +
  
  ggplot2::scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    name = "Row-scaled\nmean abundance"
  ) +
  
  ggplot2::scale_x_discrete(
    position = "top",
    drop = FALSE
  ) +
  
  ggplot2::scale_y_discrete(
    labels = display_genus
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = NULL
  ) +
  
  ggplot2::theme_minimal(
    base_size = base_size,
    base_family = base_family
  ) +
  
  ggplot2::theme(
    text = ggplot2::element_text(
      family = base_family,
      color = "black"
    ),
    
    panel.grid = ggplot2::element_blank(),
    
    axis.text.x = ggplot2::element_text(
      size = 8,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 7.0,
      color = "black"
    ),
    
    axis.ticks = ggplot2::element_blank(),
    
    legend.position = "right",
    
    legend.title = ggplot2::element_text(
      size = 7.5,
      face = "bold"
    ),
    
    legend.text = ggplot2::element_text(
      size = 7
    ),
    
    plot.margin = ggplot2::margin(
      3,
      4,
      3,
      3
    )
  )


data.table::fwrite(
  heatmap_long,
  file.path(
    source_dir,
    "Figure4C_final_29_genus_heatmap_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 8. PANEL D — ecological metrics
#
# Original final format:
# Shannon
# ↓
# Observed genera
# ↓
# Dominant genus proportion
# ==============================================================================

section(
  "PANEL D: ECOLOGICAL METRICS"
)


metric_table <- data.table::fread(
  metric_file,
  data.table = FALSE,
  check.names = FALSE
)


required_metric_columns <- c(
  "SampleID",
  "DMM_cluster",
  "Genus_Shannon",
  "Observed_genera",
  "Dominant_genus_relative_abundance"
)


missing_metric_columns <- setdiff(
  required_metric_columns,
  colnames(metric_table)
)


if (length(missing_metric_columns) > 0) {
  stop(
    paste0(
      "Metric table missing: ",
      paste(
        missing_metric_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


metric_table$Cluster <- factor(
  clean_cluster(
    metric_table$DMM_cluster
  ),
  levels = cluster_levels
)


metric_table$Shannon <- as.numeric(
  metric_table$Genus_Shannon
)

metric_table$Dominant_genus_proportion <- as.numeric(
  metric_table$Dominant_genus_relative_abundance
)


test_table <- data.frame(
  Metric = c(
    "Shannon",
    "Observed genera",
    "Dominant genus proportion"
  ),
  
  P_value = c(
    stats::kruskal.test(
      Shannon ~ Cluster,
      data = metric_table
    )$p.value,
    
    stats::kruskal.test(
      Observed_genera ~ Cluster,
      data = metric_table
    )$p.value,
    
    stats::kruskal.test(
      Dominant_genus_proportion ~ Cluster,
      data = metric_table
    )$p.value
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  test_table,
  file.path(
    source_dir,
    "Figure4D_Kruskal_Wallis_tests_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


make_metric_plot <- function(
    metric_column,
    y_label,
    p_value,
    percent_axis = FALSE
) {
  
  plot_object <- ggplot2::ggplot(
    metric_table,
    ggplot2::aes(
      x = Cluster,
      y = .data[[metric_column]],
      fill = Cluster,
      color = Cluster
    )
  ) +
    
    ggplot2::geom_boxplot(
      width = 0.52,
      linewidth = 0.42,
      outlier.shape = NA,
      alpha = 0.78,
      color = "#333333"
    ) +
    
    ggplot2::geom_jitter(
      width = 0.10,
      height = 0,
      shape = 21,
      size = 0.90,
      stroke = 0.18,
      alpha = 0.62
    ) +
    
    ggplot2::scale_fill_manual(
      values = cluster_colors,
      guide = "none",
      drop = FALSE
    ) +
    
    ggplot2::scale_color_manual(
      values = cluster_colors,
      guide = "none",
      drop = FALSE
    ) +
    
    ggplot2::labs(
      x = NULL,
      y = y_label,
      subtitle = format_p(
        p_value
      )
    ) +
    
    theme_manuscript +
    
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        size = 6.8,
        face = "bold",
        color = "black"
      ),
      
      axis.text.y = ggplot2::element_text(
        size = 6.5,
        color = "black"
      ),
      
      axis.title.y = ggplot2::element_text(
        size = 7.2,
        face = "bold",
        color = "black"
      ),
      
      plot.subtitle = ggplot2::element_text(
        size = 6.2,
        hjust = 0.5,
        margin = ggplot2::margin(
          b = 0
        )
      ),
      
      plot.margin = ggplot2::margin(
        3,
        3,
        3,
        3
      )
    )
  
  
  if (percent_axis) {
    
    plot_object <- plot_object +
      
      ggplot2::scale_y_continuous(
        labels = scales::label_percent(
          accuracy = 1
        ),
        
        expand = ggplot2::expansion(
          mult = c(
            0.04,
            0.12
          )
        )
      )
    
  } else {
    
    plot_object <- plot_object +
      
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(
          mult = c(
            0.04,
            0.12
          )
        )
      )
  }
  
  
  plot_object
}


panel_d1 <- make_metric_plot(
  metric_column = "Shannon",
  y_label = "Shannon diversity",
  p_value = test_table$P_value[
    test_table$Metric == "Shannon"
  ]
)


panel_d2 <- make_metric_plot(
  metric_column = "Observed_genera",
  y_label = "Observed genera",
  p_value = test_table$P_value[
    test_table$Metric == "Observed genera"
  ]
)


panel_d3 <- make_metric_plot(
  metric_column = "Dominant_genus_proportion",
  y_label = "Dominant genus proportion",
  p_value = test_table$P_value[
    test_table$Metric ==
      "Dominant genus proportion"
  ],
  percent_axis = TRUE
)


panel_d <- panel_d1 /
  panel_d2 /
  panel_d3


panel_d <- panel_d +
  patchwork::plot_layout(
    heights = c(
      1,
      1,
      1
    )
  )


# ==============================================================================
# 9. PANEL E — paired nonCAC -> CAC state correspondence
#
# Uses original Figure4E plotting parameters:
# alpha = 0.78
# knot.pos = 0.45
# linewidth = 0.32
# discern = TRUE
# ==============================================================================

section(
  "PANEL E: PAIRED DMM STATE CORRESPONDENCE"
)


paired_data <- data.table::fread(
  paired_transition_file,
  data.table = FALSE,
  check.names = FALSE
)


required_pair_columns <- c(
  "PairID",
  "nonCA_state",
  "CA_state"
)


missing_pair_columns <- setdiff(
  required_pair_columns,
  colnames(paired_data)
)


if (length(missing_pair_columns) > 0) {
  stop(
    paste0(
      "Paired transition table missing: ",
      paste(
        missing_pair_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


paired_data$nonCA_state <- clean_cluster(
  paired_data$nonCA_state
)

paired_data$CA_state <- clean_cluster(
  paired_data$CA_state
)


if (nrow(paired_data) != 23) {
  stop(
    paste0(
      "Expected 23 paired subjects, found ",
      nrow(paired_data),
      "."
    ),
    call. = FALSE
  )
}


alluvial_data <- paired_data %>%
  dplyr::count(
    nonCA_state,
    CA_state,
    name = "Number_of_pairs"
  )


alluvial_data$nonCA_state <- factor(
  alluvial_data$nonCA_state,
  levels = cluster_levels
)

alluvial_data$CA_state <- factor(
  alluvial_data$CA_state,
  levels = cluster_levels
)


n_pairs <- sum(
  alluvial_data$Number_of_pairs
)


if (n_pairs != 23) {
  stop(
    "Paired alluvial total is not 23.",
    call. = FALSE
  )
}


panel_e <- ggplot2::ggplot(
  alluvial_data,
  ggplot2::aes(
    axis1 = nonCA_state,
    axis2 = CA_state,
    y = Number_of_pairs
  )
) +
  
  ggalluvial::geom_alluvium(
    ggplot2::aes(
      fill = nonCA_state
    ),
    
    width = 0.16,
    alpha = 0.78,
    knot.pos = 0.45,
    linewidth = 0.32,
    color = "white",
    discern = TRUE
  ) +
  
  ggalluvial::geom_stratum(
    width = 0.18,
    fill = "white",
    color = "#333333",
    linewidth = 0.50,
    discern = TRUE
  ) +
  
  ggplot2::geom_text(
    stat = "stratum",
    
    ggplot2::aes(
      label = sub(
        "\\.[0-9]+$",
        "",
        after_stat(stratum)
      )
    ),
    
    size = 3.8,
    fontface = "bold",
    family = base_family,
    discern = TRUE
  ) +
  
  ggplot2::scale_x_discrete(
    limits = c(
      "Paired nonCAC",
      "CAC"
    ),
    
    expand = c(
      0.10,
      0.10
    )
  ) +
  
  ggplot2::scale_y_continuous(
    breaks = seq(
      0,
      n_pairs,
      by = 5
    ),
    
    expand = ggplot2::expansion(
      mult = c(
        0,
        0.025
      )
    )
  ) +
  
  ggplot2::scale_fill_manual(
    values = cluster_colors,
    breaks = cluster_levels,
    name = "DMM state"
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Number of paired patients"
  ) +
  
  ggplot2::theme_classic(
    base_size = 10,
    base_family = base_family
  ) +
  
  ggplot2::theme(
    text = ggplot2::element_text(
      family = base_family,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      color = "black",
      linewidth = 0.45
    ),
    
    axis.ticks = ggplot2::element_line(
      color = "black",
      linewidth = 0.40
    ),
    
    axis.text.x = ggplot2::element_text(
      size = 9.2,
      face = "bold",
      color = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 8.3,
      color = "black"
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 9,
      face = "bold"
    ),
    
    legend.position = "none",
    
    plot.margin = ggplot2::margin(
      t = 4,
      r = 3,
      b = 3,
      l = 3,
      unit = "pt"
    )
  )


# ==============================================================================
# 10. Save standalone panels
# ==============================================================================

saveRDS(
  panel_a,
  file.path(
    panel_dir,
    "Figure4A_7KB.rds"
  )
)

saveRDS(
  panel_b,
  file.path(
    panel_dir,
    "Figure4B_7KB.rds"
  )
)

saveRDS(
  panel_c,
  file.path(
    panel_dir,
    "Figure4C_7KB.rds"
  )
)

saveRDS(
  panel_d,
  file.path(
    panel_dir,
    "Figure4D_7KB.rds"
  )
)

saveRDS(
  panel_e,
  file.path(
    panel_dir,
    "Figure4E_7KB.rds"
  )
)


# ==============================================================================
# 11. FINAL ASSEMBLY
#
# Original final assembly:
#
# top:
# A | B | E
#
# bottom:
# C     | D
#
# widths top:    0.78, 1.22, 1.45
# widths bottom: 1.72, 1.00
# heights:       0.78, 1.52
# ==============================================================================

section(
  "FINAL FIGURE 4 ASSEMBLY"
)


shared_margin <- ggplot2::theme(
  plot.margin = ggplot2::margin(
    t = 8,
    r = 7,
    b = 6,
    l = 7,
    unit = "pt"
  )
)


# ------------------------------------------------------------------------------
# A
# ------------------------------------------------------------------------------

panel_a_final <- panel_a +
  
  ggplot2::labs(
    tag = "A"
  ) +
  
  shared_margin +
  
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = base_family,
      size = 13,
      face = "bold",
      color = "black"
    ),
    
    axis.title = ggplot2::element_text(
      family = base_family,
      size = 8.5,
      face = "bold"
    ),
    
    axis.text = ggplot2::element_text(
      family = base_family,
      size = 7.5,
      color = "black"
    )
  )


# ------------------------------------------------------------------------------
# B
# ------------------------------------------------------------------------------

panel_b_final <- panel_b +
  
  ggplot2::labs(
    tag = "B"
  ) +
  
  shared_margin +
  
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = base_family,
      size = 13,
      face = "bold",
      color = "black"
    ),
    
    axis.title = ggplot2::element_text(
      family = base_family,
      size = 8.5,
      face = "bold"
    ),
    
    axis.text.x = ggplot2::element_text(
      family = base_family,
      size = 7,
      color = "black",
      angle = 30,
      hjust = 1
    ),
    
    axis.text.y = ggplot2::element_text(
      family = base_family,
      size = 7.5,
      color = "black"
    ),
    
    legend.title = ggplot2::element_text(
      family = base_family,
      size = 8,
      face = "bold"
    ),
    
    legend.text = ggplot2::element_text(
      family = base_family,
      size = 7
    )
  )


# ------------------------------------------------------------------------------
# C
# ------------------------------------------------------------------------------

panel_c_final <- panel_c +
  
  ggplot2::labs(
    tag = "C"
  ) +
  
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = base_family,
      size = 13,
      face = "bold",
      color = "black"
    ),
    
    axis.text.x = ggplot2::element_text(
      family = base_family,
      size = 7.5,
      color = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      family = base_family,
      size = 7.2,
      color = "black"
    ),
    
    legend.title = ggplot2::element_text(
      family = base_family,
      size = 8,
      face = "bold"
    ),
    
    legend.text = ggplot2::element_text(
      family = base_family,
      size = 7
    ),
    
    plot.margin = ggplot2::margin(
      t = 8,
      r = 9,
      b = 6,
      l = 7,
      unit = "pt"
    )
  )


# ------------------------------------------------------------------------------
# D
#
# Keep internal three-plot vertical patchwork.
# ------------------------------------------------------------------------------

panel_d_final <- panel_d +
  
  patchwork::plot_annotation(
    tag_prefix = "D",
    
    theme = ggplot2::theme(
      plot.tag = ggplot2::element_text(
        family = base_family,
        size = 13,
        face = "bold",
        color = "black"
      ),
      
      plot.tag.position = c(
        0.005,
        0.995
      ),
      
      plot.margin = ggplot2::margin(
        t = 7,
        r = 7,
        b = 5,
        l = 7,
        unit = "pt"
      )
    )
  )


# ------------------------------------------------------------------------------
# E
# ------------------------------------------------------------------------------

panel_e_final <- panel_e +
  
  ggplot2::labs(
    tag = "E"
  ) +
  
  shared_margin +
  
  ggplot2::theme(
    plot.tag = ggplot2::element_text(
      family = base_family,
      size = 13,
      face = "bold",
      color = "black"
    ),
    
    axis.title = ggplot2::element_text(
      family = base_family,
      size = 8.5,
      face = "bold"
    ),
    
    axis.text.x = ggplot2::element_text(
      family = base_family,
      size = 7.8,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    
    axis.text.y = ggplot2::element_text(
      family = base_family,
      size = 7.5,
      color = "black"
    ),
    
    legend.position = "none"
  )


# ==============================================================================
# 12. Layout
# ==============================================================================

top_row <- panel_a_final |
  panel_b_final |
  panel_e_final


top_row <- top_row +
  
  patchwork::plot_layout(
    widths = c(
      0.78,
      1.22,
      1.45
    )
  )


bottom_row <- panel_c_final |
  panel_d_final


bottom_row <- bottom_row +
  
  patchwork::plot_layout(
    widths = c(
      1.72,
      1.00
    )
  )


figure4_final <- top_row /
  bottom_row


figure4_final <- figure4_final +
  
  patchwork::plot_layout(
    heights = c(
      0.78,
      1.52
    )
  ) +
  
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.background = ggplot2::element_rect(
        fill = "white",
        color = NA
      ),
      
      plot.margin = ggplot2::margin(
        t = 7,
        r = 7,
        b = 7,
        l = 7,
        unit = "pt"
      )
    )
  )


# ==============================================================================
# 13. Final output
# ==============================================================================

output_stem <- file.path(
  output_dir,
  "Figure4_DMM_states_with_paired_E_FINAL_7KB"
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


figure_width <- 12.5
figure_height <- 9.2
figure_dpi <- 600


# ==============================================================================
# 14. PDF
# ==============================================================================

grDevices::pdf(
  file = pdf_file,
  width = figure_width,
  height = figure_height,
  family = base_family,
  useDingbats = FALSE,
  onefile = TRUE
)

print(
  figure4_final
)

grDevices::dev.off()


# ==============================================================================
# 15. PNG
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
  
} else {
  
  grDevices::png(
    filename = png_file,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    bg = "white",
    type = "quartz"
  )
}


print(
  figure4_final
)

grDevices::dev.off()


# ==============================================================================
# 16. TIFF
# ==============================================================================

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_tiff(
    filename = tiff_file,
    width = figure_width,
    height = figure_height,
    units = "in",
    res = figure_dpi,
    compression = "lzw",
    background = "white"
  )
  
} else {
  
  grDevices::tiff(
    filename = tiff_file,
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
  figure4_final
)

grDevices::dev.off()


# ==============================================================================
# 17. RDS
# ==============================================================================

saveRDS(
  figure4_final,
  rds_file
)


# ==============================================================================
# 18. Final audit
# ==============================================================================

section(
  "FINAL AUDIT"
)


cat(
  "Panel A optimal K: 3\n"
)

cat(
  "Panel A Laplace K3: 86176\n"
)

cat(
  "\nPanel B CAC C2:\n"
)

cat(
  cac_c2$Count,
  "/",
  cac_c2$Total,
  " = ",
  sprintf(
    "%.1f%%",
    100 *
      cac_c2$Proportion
  ),
  "\n",
  sep = ""
)


cat(
  "\nPanel C genera:\n"
)

cat(
  "C1 = ",
  length(c1_genera),
  "\n",
  sep = ""
)

cat(
  "C2 = ",
  length(c2_genera),
  "\n",
  sep = ""
)

cat(
  "C3 = ",
  length(c3_genera),
  "\n",
  sep = ""
)

cat(
  "Total = ",
  length(final_genus_order),
  "\n",
  sep = ""
)


cat(
  "\nPanel D P values:\n"
)

print(
  test_table,
  row.names = FALSE
)


cat(
  "\nPanel E transition matrix:\n"
)

print(
  table(
    factor(
      paired_data$nonCA_state,
      levels = cluster_levels
    ),
    
    factor(
      paired_data$CA_state,
      levels = cluster_levels
    )
  )
)


cat(
  "\nFinal outputs:\n"
)

cat(
  pdf_file,
  "\n"
)

cat(
  png_file,
  "\n"
)

cat(
  tiff_file,
  "\n"
)

cat(
  rds_file,
  "\n"
)

cat(
  "\nFigure 4 all-in-one generation completed.\n"
)