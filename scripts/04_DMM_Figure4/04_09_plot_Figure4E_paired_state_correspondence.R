#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_09_plot_Figure4E_paired_state_correspondence.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Purpose:
## Generate Figure 4E for the fixed paired nonCAC23 -> CAC23 DMM states.
##
## Validated transition matrix from script 14:
##
##                    CAC
## nonCAC        C1   C2   C3
## C1             2    4    2
## C2             0   11    1
## C3             0    2    1
##
## Total paired patients = 23
##
## IMPORTANT:
## ggalluvial is used in ALLUVIA/WIDE form:
##   axis1 = nonCAC
##   axis2 = CAC
##   y     = N
##
## Do NOT add an alluvium aesthetic here.
############################################################


options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "ggplot2",
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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
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
  "projection_CA23_nonCA23_fixed_K3_7KB",
  "tables",
  "DMM_fixed_K3_paired_state_transition_source_7KB.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "Figure4_final_panels",
  "Figure4E_paired_DMM_state_correspondence"
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
# 2. Input check
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
# 3. Read paired transition source
# ==============================================================================

dat <- data.table::fread(
  input_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "PairID",
  "nonCA_SampleID",
  "nonCA_state",
  "CA_SampleID",
  "CA_state"
)

missing_columns <- setdiff(
  required_columns,
  colnames(dat)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 4. Strict paired audit
# ==============================================================================

if (nrow(dat) != 23) {
  stop(
    paste0(
      "Expected exactly 23 paired patients, found ",
      nrow(dat),
      "."
    ),
    call. = FALSE
  )
}

if (anyDuplicated(dat$PairID) > 0) {
  stop(
    "Duplicated PairID detected.",
    call. = FALSE
  )
}

if (anyDuplicated(dat$nonCA_SampleID) > 0) {
  stop(
    "Duplicated nonCAC SampleID detected.",
    call. = FALSE
  )
}

if (anyDuplicated(dat$CA_SampleID) > 0) {
  stop(
    "Duplicated CAC SampleID detected.",
    call. = FALSE
  )
}


# ==============================================================================
# 5. State labels
# ==============================================================================

state_order <- c(
  "C1",
  "C2",
  "C3"
)

dat$nonCA_state <- as.character(
  dat$nonCA_state
)

dat$CA_state <- as.character(
  dat$CA_state
)

if (!all(dat$nonCA_state %in% state_order)) {
  stop(
    "Unexpected nonCAC DMM state label detected.",
    call. = FALSE
  )
}

if (!all(dat$CA_state %in% state_order)) {
  stop(
    "Unexpected CAC DMM state label detected.",
    call. = FALSE
  )
}


# ==============================================================================
# 6. Observed transition matrix
# ==============================================================================

transition_matrix <- table(
  factor(
    dat$nonCA_state,
    levels = state_order
  ),
  factor(
    dat$CA_state,
    levels = state_order
  )
)

dimnames(transition_matrix) <- list(
  nonCAC = state_order,
  CAC = state_order
)

cat("\n")
cat("Observed paired transition matrix:\n")

print(
  transition_matrix
)


# ==============================================================================
# 7. Validate transition matrix
# ==============================================================================

expected_transition_matrix <- matrix(
  c(
    2L, 4L, 2L,
    0L, 11L, 1L,
    0L, 2L, 1L
  ),
  nrow = 3,
  ncol = 3,
  byrow = TRUE,
  dimnames = list(
    nonCAC = state_order,
    CAC = state_order
  )
)

observed_numeric <- matrix(
  as.integer(
    transition_matrix
  ),
  nrow = 3,
  ncol = 3,
  dimnames = list(
    nonCAC = state_order,
    CAC = state_order
  )
)

if (
  !all(
    observed_numeric ==
    expected_transition_matrix
  )
) {
  
  cat("\n")
  cat("Expected transition matrix:\n")
  
  print(
    expected_transition_matrix
  )
  
  cat("\n")
  cat("Observed transition matrix:\n")
  
  print(
    observed_numeric
  )
  
  stop(
    "Observed transition values differ from the validated 7KB projection.",
    call. = FALSE
  )
}

cat("\n")
cat("Transition-matrix validation: PASS\n")


# ==============================================================================
# 8. State-total audit
# ==============================================================================

nonca_state_totals <- table(
  factor(
    dat$nonCA_state,
    levels = state_order
  )
)

ca_state_totals <- table(
  factor(
    dat$CA_state,
    levels = state_order
  )
)

expected_nonca_totals <- c(
  C1 = 8L,
  C2 = 12L,
  C3 = 3L
)

expected_ca_totals <- c(
  C1 = 2L,
  C2 = 17L,
  C3 = 4L
)

if (
  !all(
    as.integer(nonca_state_totals) ==
    as.integer(expected_nonca_totals)
  )
) {
  stop(
    "nonCAC state-total validation failed.",
    call. = FALSE
  )
}

if (
  !all(
    as.integer(ca_state_totals) ==
    as.integer(expected_ca_totals)
  )
) {
  stop(
    "CAC state-total validation failed.",
    call. = FALSE
  )
}

cat("\n")
cat("State-total validation: PASS\n")


# ==============================================================================
# 9. Export transition matrix
# ==============================================================================

transition_matrix_df <- data.frame(
  nonCAC_state = state_order,
  C1 = observed_numeric[, "C1"],
  C2 = observed_numeric[, "C2"],
  C3 = observed_numeric[, "C3"],
  stringsAsFactors = FALSE
)

data.table::fwrite(
  transition_matrix_df,
  file.path(
    source_dir,
    "Figure4E_transition_matrix_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 10. Exact McNemar test
#
# Binary comparison:
# C2 versus non-C2
#
# nonCAC non-C2 -> CAC C2:
# C1 -> C2 = 4
# C3 -> C2 = 2
# total       = 6
#
# nonCAC C2 -> CAC non-C2:
# C2 -> C1 = 0
# C2 -> C3 = 1
# total       = 1
# ==============================================================================

nonca_is_c2 <- dat$nonCA_state == "C2"

ca_is_c2 <- dat$CA_state == "C2"

nonca_non_c2_to_ca_c2 <- sum(
  !nonca_is_c2 &
    ca_is_c2
)

nonca_c2_to_ca_non_c2 <- sum(
  nonca_is_c2 &
    !ca_is_c2
)

discordant_pairs <-
  nonca_non_c2_to_ca_c2 +
  nonca_c2_to_ca_non_c2

if (discordant_pairs > 0) {
  
  smaller_discordant_count <- min(
    nonca_non_c2_to_ca_c2,
    nonca_c2_to_ca_non_c2
  )
  
  exact_mcnemar_p <- 2 *
    stats::pbinom(
      smaller_discordant_count,
      size = discordant_pairs,
      prob = 0.5
    )
  
  exact_mcnemar_p <- min(
    exact_mcnemar_p,
    1
  )
  
} else {
  
  exact_mcnemar_p <- 1
}

mcnemar_summary <- data.frame(
  Contrast = "C2 versus non-C2",
  nonCAC_nonC2_to_CAC_C2 = nonca_non_c2_to_ca_c2,
  nonCAC_C2_to_CAC_nonC2 = nonca_c2_to_ca_non_c2,
  Discordant_pairs = discordant_pairs,
  Exact_McNemar_P = exact_mcnemar_p,
  stringsAsFactors = FALSE
)

data.table::fwrite(
  mcnemar_summary,
  file.path(
    source_dir,
    "Figure4E_exact_McNemar_C2_vs_nonC2_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 11. Build ALLUVIA/WIDE transition-count table
#
# One row = one nonCAC -> CAC transition type.
#
# IMPORTANT:
# This is the correct format for:
#
# aes(
#   axis1 = nonCAC,
#   axis2 = CAC,
#   y = N
# )
#
# No alluvium aesthetic is used.
# ==============================================================================

flow_data <- dat %>%
  count(
    nonCA_state,
    CA_state,
    name = "N"
  ) %>%
  tidyr::complete(
    nonCA_state = state_order,
    CA_state = state_order,
    fill = list(
      N = 0L
    )
  ) %>%
  filter(
    N > 0
  )

flow_data$nonCA_state <- factor(
  flow_data$nonCA_state,
  levels = state_order
)

flow_data$CA_state <- factor(
  flow_data$CA_state,
  levels = state_order
)

flow_data <- flow_data %>%
  arrange(
    nonCA_state,
    CA_state
  )

flow_data$Flow_state <- flow_data$nonCA_state


data.table::fwrite(
  flow_data,
  file.path(
    source_dir,
    "Figure4E_flow_counts_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 12. Check flow total
# ==============================================================================

if (sum(flow_data$N) != 23) {
  stop(
    paste0(
      "Alluvial flow total should be 23, found ",
      sum(flow_data$N),
      "."
    ),
    call. = FALSE
  )
}

cat("\n")
cat("Alluvial flow-count validation: PASS\n")


# ==============================================================================
# 13. State colors
#
# Temporary standalone-panel palette.
# Final full Figure 4 assembly will enforce the original manuscript style.
# ==============================================================================

state_colors <- c(
  "C1" = "#6C8EBF",
  "C2" = "#C95E57",
  "C3" = "#88B27C"
)


# ==============================================================================
# 14. Figure 4E
#
# Correct ggalluvial ALLUVIA/WIDE syntax:
#
# axis1 = nonCAC state
# axis2 = CAC state
# y     = patient count
#
# Do NOT specify alluvium = ...
# ==============================================================================

panel_e <- ggplot2::ggplot(
  flow_data,
  ggplot2::aes(
    axis1 = nonCA_state,
    axis2 = CA_state,
    y = N
  )
) +
  
  ggalluvial::geom_alluvium(
    ggplot2::aes(
      fill = Flow_state
    ),
    width = 0.16,
    alpha = 0.68,
    knot.pos = 0.5
  ) +
  
  ggalluvial::geom_stratum(
    width = 0.16,
    fill = "white",
    color = "black",
    linewidth = 0.40
  ) +
  
  ggplot2::geom_text(
    stat = "stratum",
    ggplot2::aes(
      label = after_stat(stratum)
    ),
    family = "Helvetica",
    size = 2.8,
    color = "black"
  ) +
  
  ggplot2::scale_x_discrete(
    limits = c(
      "Paired nonCAC",
      "CAC"
    ),
    expand = c(
      0.18,
      0.18
    )
  ) +
  
  ggplot2::scale_fill_manual(
    values = state_colors,
    breaks = state_order,
    drop = FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    limits = c(
      0,
      23
    ),
    breaks = c(
      0,
      5,
      10,
      15,
      20
    ),
    expand = c(
      0,
      0
    )
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Number of paired patients"
  ) +
  
  ggplot2::guides(
    fill = "none"
  ) +
  
  ggplot2::theme_classic(
    base_size = 8,
    base_family = "Helvetica"
  ) +
  
  ggplot2::theme(
    text = ggplot2::element_text(
      family = "Helvetica",
      color = "black"
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 8,
      color = "black"
    ),
    
    axis.text.y = ggplot2::element_text(
      size = 7.5,
      color = "black"
    ),
    
    axis.text.x = ggplot2::element_text(
      size = 8,
      color = "black"
    ),
    
    axis.line = ggplot2::element_line(
      linewidth = 0.35,
      color = "black"
    ),
    
    axis.ticks = ggplot2::element_line(
      linewidth = 0.35,
      color = "black"
    ),
    
    plot.margin = ggplot2::margin(
      5,
      5,
      4,
      5,
      unit = "pt"
    )
  )


# ==============================================================================
# 15. Save RDS
# ==============================================================================

saveRDS(
  panel_e,
  file.path(
    output_dir,
    "Figure4E_paired_DMM_state_correspondence_7KB.rds"
  )
)


# ==============================================================================
# 16. Export PDF
#
# Base PDF only.
# No Cairo.
# ==============================================================================

pdf_file <- file.path(
  output_dir,
  "Figure4E_paired_DMM_state_correspondence_7KB.pdf"
)

grDevices::pdf(
  file = pdf_file,
  width = 4.2,
  height = 4.0,
  family = "Helvetica",
  useDingbats = FALSE
)

print(
  panel_e
)

grDevices::dev.off()


# ==============================================================================
# 17. Export PNG
# ==============================================================================

png_file <- file.path(
  output_dir,
  "Figure4E_paired_DMM_state_correspondence_7KB.png"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_png(
    filename = png_file,
    width = 4.2,
    height = 4.0,
    units = "in",
    res = 600,
    background = "white"
  )
  
  print(
    panel_e
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::png(
    filename = png_file,
    width = 4.2,
    height = 4.0,
    units = "in",
    res = 600,
    type = "quartz",
    bg = "white"
  )
  
  print(
    panel_e
  )
  
  grDevices::dev.off()
}


# ==============================================================================
# 18. Export TIFF
# ==============================================================================

tiff_file <- file.path(
  output_dir,
  "Figure4E_paired_DMM_state_correspondence_7KB.tif"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_tiff(
    filename = tiff_file,
    width = 4.2,
    height = 4.0,
    units = "in",
    res = 600,
    compression = "lzw",
    background = "white"
  )
  
  print(
    panel_e
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::tiff(
    filename = tiff_file,
    width = 4.2,
    height = 4.0,
    units = "in",
    res = 600,
    type = "quartz",
    compression = "lzw",
    bg = "white"
  )
  
  print(
    panel_e
  )
  
  grDevices::dev.off()
}


# ==============================================================================
# 19. Console summary
# ==============================================================================

cat("\n")
cat("============================================================\n")
cat("Final 7KB Figure 4E completed successfully.\n")
cat("============================================================\n")

cat("\n")
cat("Validated paired transition matrix:\n")

print(
  transition_matrix
)

cat("\n")
cat("Paired state totals:\n")

cat(
  "nonCAC: C1=",
  as.integer(nonca_state_totals["C1"]),
  ", C2=",
  as.integer(nonca_state_totals["C2"]),
  ", C3=",
  as.integer(nonca_state_totals["C3"]),
  "\n",
  sep = ""
)

cat(
  "CAC: C1=",
  as.integer(ca_state_totals["C1"]),
  ", C2=",
  as.integer(ca_state_totals["C2"]),
  ", C3=",
  as.integer(ca_state_totals["C3"]),
  "\n",
  sep = ""
)

cat("\n")
cat("Transition counts used for alluvial plot:\n")

print(
  flow_data
)

cat("\n")
cat("Exact McNemar test, C2 versus non-C2:\n")

print(
  mcnemar_summary,
  row.names = FALSE
)

cat("\n")
cat("Output directory:\n")

cat(
  output_dir,
  "\n"
)