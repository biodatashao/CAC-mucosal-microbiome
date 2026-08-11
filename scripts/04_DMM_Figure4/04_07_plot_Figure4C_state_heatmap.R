#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 04_07_plot_Figure4C_state_heatmap.R
##
## Module 04 - Dirichlet multinomial community states (Figure 4)
##
## Script:
##
## Purpose:
## Generate final 7KB Figure 4C using the exact final-manuscript
## selection and plotting logic from:
##
## IMPORTANT:
## - C1 uses the predefined six-genus candidate list, BUT only
##   genera already present in the qualified selected-genera
##   source table are retained, exactly as in the original script.
## - C2 keeps first 13 selected genera.
## - C3 keeps first 10 selected genera.
## - No genus is manually restored if it fails the 7KB upstream
##   selection thresholds.
############################################################


# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "data.table",
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
      paste(
        missing_packages,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# ============================================================
# 1. Paths
# ============================================================

project_root <- PROJECT_ROOT

source_root <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "Figure4C_Lavelle_cluster_mean_heatmap_7KB",
  "source_tables"
)

heatmap_file <- file.path(
  source_root,
  "Figure3B_row_normalized_cluster_mean_source.tsv"
)

selected_genera_file <- file.path(
  source_root,
  "Figure3B_selected_top13_per_cluster.tsv"
)

output_dir <- file.path(
  project_root,
  "output",
  "analysis",
  "03_Figure4_DMM",
  "Figure4_final_panels",
  "Figure4C_DMM_cluster_mean_heatmap"
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


# ============================================================
# 2. Validate input files
# ============================================================

required_files <- c(
  heatmap_file,
  selected_genera_file
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Required file(s) not found:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}


# ============================================================
# 3. Original final-manuscript style
# ============================================================

base_family <- "Helvetica"

base_size <- 8

cluster_levels <- c(
  "C1",
  "C2",
  "C3"
)

heatmap_low <- "#5579B6"

heatmap_mid <- "#F7F7F7"

heatmap_high <- "#C95E57"


# ============================================================
# 4. Helper
# ============================================================

clean_cluster <- function(x) {
  
  x <- trimws(
    as.character(x)
  )
  
  x_upper <- toupper(x)
  
  output <- rep(
    NA_character_,
    length(x_upper)
  )
  
  output[
    grepl(
      "C1|CLUSTER.?1|COMPONENT.?1|STATE.?1|^1$",
      x_upper
    )
  ] <- "C1"
  
  output[
    grepl(
      "C2|CLUSTER.?2|COMPONENT.?2|STATE.?2|^2$",
      x_upper
    )
  ] <- "C2"
  
  output[
    grepl(
      "C3|CLUSTER.?3|COMPONENT.?3|STATE.?3|^3$",
      x_upper
    )
  ] <- "C3"
  
  output
}


# ============================================================
# 5. Read source data
# ============================================================

heatmap_raw <- data.table::fread(
  heatmap_file,
  data.table = FALSE,
  check.names = FALSE
)

selected_raw <- data.table::fread(
  selected_genera_file,
  data.table = FALSE,
  check.names = FALSE
)


# ============================================================
# 6. Validate required columns
# ============================================================

required_heatmap_columns <- c(
  "Genus",
  "DMM_cluster",
  "Row_normalized_mean"
)

missing_heatmap_columns <- setdiff(
  required_heatmap_columns,
  colnames(
    heatmap_raw
  )
)

if (length(missing_heatmap_columns) > 0) {
  stop(
    paste0(
      "Missing heatmap columns: ",
      paste(
        missing_heatmap_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


required_selected_columns <- c(
  "Genus",
  "Assigned_cluster",
  "Display_order"
)

missing_selected_columns <- setdiff(
  required_selected_columns,
  colnames(
    selected_raw
  )
)

if (length(missing_selected_columns) > 0) {
  stop(
    paste0(
      "Missing selected-genera columns: ",
      paste(
        missing_selected_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# ============================================================
# 7. Reproduce original selected-order table
# ============================================================

selected_order <- selected_raw[
  ,
  c(
    "Genus",
    "Assigned_cluster",
    "Display_order"
  ),
  drop = FALSE
]

selected_order[["Assigned_cluster"]] <- clean_cluster(
  selected_order[["Assigned_cluster"]]
)

selected_order <- selected_order[
  selected_order[["Assigned_cluster"]] %in%
    cluster_levels,
  ,
  drop = FALSE
]

selected_order <- selected_order[
  order(
    match(
      selected_order[["Assigned_cluster"]],
      cluster_levels
    ),
    selected_order[["Display_order"]]
  ),
  ,
  drop = FALSE
]


# ============================================================
# 8. Exact final-manuscript genus selection
#
# C1:
# predefined representative genera, but ONLY when present in
# the upstream selected-genera table.
#
# C2:
# first 13 selected genera.
#
# C3:
# first 10 selected genera.
# ============================================================

c1_keep_genera <- c(
  "Enterococcus",
  "Neobacillus",
  "Clostridium",
  "Streptococcus",
  "Actinomyces",
  "Lachnospiraceae_NK3A20_group"
)


# ------------------------------------------------------------
# C1
# ------------------------------------------------------------

selected_c1 <- selected_order[
  selected_order[["Assigned_cluster"]] == "C1" &
    selected_order[["Genus"]] %in%
    c1_keep_genera,
  ,
  drop = FALSE
]

missing_c1_genera <- setdiff(
  c1_keep_genera,
  selected_c1[["Genus"]]
)

selected_c1[["Manual_order"]] <- match(
  selected_c1[["Genus"]],
  c1_keep_genera
)

selected_c1 <- selected_c1[
  order(
    selected_c1[["Manual_order"]]
  ),
  ,
  drop = FALSE
]

selected_c1[["Manual_order"]] <- NULL


# ------------------------------------------------------------
# C2
# ------------------------------------------------------------

selected_c2 <- selected_order[
  selected_order[["Assigned_cluster"]] == "C2",
  ,
  drop = FALSE
]

selected_c2 <- head(
  selected_c2,
  13
)


# ------------------------------------------------------------
# C3
# ------------------------------------------------------------

selected_c3 <- selected_order[
  selected_order[["Assigned_cluster"]] == "C3",
  ,
  drop = FALSE
]

selected_c3 <- head(
  selected_c3,
  10
)


# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------

selected_order_main <- rbind(
  selected_c1,
  selected_c2,
  selected_c3
)

selected_order_main[["Display_order_main"]] <- seq_len(
  nrow(
    selected_order_main
  )
)


# ============================================================
# 9. Strict audit
# ============================================================

main_cluster_counts <- table(
  factor(
    selected_order_main[["Assigned_cluster"]],
    levels = cluster_levels
  )
)

cat("\n")
cat("Main-figure heatmap genera by state:\n")

print(
  main_cluster_counts
)

cat("\n")
cat("C1 predefined genera requested:\n")

print(
  c1_keep_genera
)

cat("\n")
cat("C1 genera retained after 7KB upstream qualification:\n")

print(
  selected_c1[["Genus"]]
)

cat("\n")
cat("C1 genera absent from the qualified 7KB selection table:\n")

if (length(missing_c1_genera) == 0) {
  
  cat("None\n")
  
} else {
  
  print(
    missing_c1_genera
  )
}


if (
  nrow(
    selected_c2
  ) != 13
) {
  warning(
    paste0(
      "C2 contains ",
      nrow(
        selected_c2
      ),
      " genera rather than 13."
    )
  )
}


if (
  nrow(
    selected_c3
  ) != 10
) {
  warning(
    paste0(
      "C3 contains ",
      nrow(
        selected_c3
      ),
      " genera rather than 10."
    )
  )
}


# ============================================================
# 10. Prepare heatmap long table
#
# Exact original logic.
# ============================================================

heatmap_long <- data.frame(
  Genus = as.character(
    heatmap_raw[["Genus"]]
  ),
  
  Cluster = clean_cluster(
    heatmap_raw[["DMM_cluster"]]
  ),
  
  Value = as.numeric(
    heatmap_raw[["Row_normalized_mean"]]
  ),
  
  stringsAsFactors = FALSE
)


heatmap_long <- merge(
  heatmap_long,
  selected_order_main,
  by = "Genus",
  all = FALSE,
  sort = FALSE
)


heatmap_long[["Cluster"]] <- factor(
  heatmap_long[["Cluster"]],
  levels = cluster_levels
)


genus_order <- selected_order_main[["Genus"]]


heatmap_long[["Genus"]] <- factor(
  heatmap_long[["Genus"]],
  levels = rev(
    genus_order
  )
)


# ============================================================
# 11. Audit dimensions
# ============================================================

expected_heatmap_rows <- length(
  genus_order
) * 3L

if (
  nrow(
    heatmap_long
  ) !=
  expected_heatmap_rows
) {
  
  stop(
    paste0(
      "Heatmap source audit failed: expected ",
      expected_heatmap_rows,
      " rows for ",
      length(
        genus_order
      ),
      " genera × 3 states, but found ",
      nrow(
        heatmap_long
      ),
      "."
    ),
    call. = FALSE
  )
}


# ============================================================
# 12. State separator positions
#
# Exact original logic.
# ============================================================

cluster_counts <- table(
  factor(
    selected_order_main[["Assigned_cluster"]],
    levels = cluster_levels
  )
)


cluster_boundaries <- cumsum(
  cluster_counts
)


separator_y <- (
  length(
    genus_order
  ) -
    cluster_boundaries[1:2] +
    0.5
)


# ============================================================
# 13. Final Figure 4C
#
# Exact style from final manuscript script.
# ============================================================

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


# ============================================================
# 14. Export final source tables
# ============================================================

data.table::fwrite(
  heatmap_long,
  file.path(
    source_dir,
    "Figure4C_heatmap_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  selected_order_main,
  file.path(
    source_dir,
    "Figure4C_selected_genera_final_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


c1_audit <- data.frame(
  Genus = c1_keep_genera,
  
  Present_in_qualified_7KB_table =
    c1_keep_genera %in%
    selected_c1[["Genus"]],
  
  stringsAsFactors = FALSE
)

data.table::fwrite(
  c1_audit,
  file.path(
    source_dir,
    "Figure4C_C1_predefined_genera_audit_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ============================================================
# 15. Save RDS
# ============================================================

saveRDS(
  panel_c,
  file.path(
    output_dir,
    "Figure4C_DMM_cluster_mean_heatmap_7KB.rds"
  )
)


# ============================================================
# 16. Export PDF
#
# Base PDF only. No Cairo/X11.
# ============================================================

pdf_file <- file.path(
  output_dir,
  "Figure4C_DMM_cluster_mean_heatmap_7KB.pdf"
)

grDevices::pdf(
  file = pdf_file,
  width = 5.2,
  height = 7.0,
  family = "Helvetica",
  useDingbats = FALSE
)

print(
  panel_c
)

grDevices::dev.off()


# ============================================================
# 17. Export PNG
# ============================================================

png_file <- file.path(
  output_dir,
  "Figure4C_DMM_cluster_mean_heatmap_7KB.png"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_png(
    filename = png_file,
    width = 5.2,
    height = 7.0,
    units = "in",
    res = 600,
    background = "white"
  )
  
  print(
    panel_c
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::png(
    filename = png_file,
    width = 5.2,
    height = 7.0,
    units = "in",
    res = 600,
    type = "quartz",
    bg = "white"
  )
  
  print(
    panel_c
  )
  
  grDevices::dev.off()
}


# ============================================================
# 18. Export TIFF
# ============================================================

tiff_file <- file.path(
  output_dir,
  "Figure4C_DMM_cluster_mean_heatmap_7KB.tif"
)

if (
  requireNamespace(
    "ragg",
    quietly = TRUE
  )
) {
  
  ragg::agg_tiff(
    filename = tiff_file,
    width = 5.2,
    height = 7.0,
    units = "in",
    res = 600,
    compression = "lzw",
    background = "white"
  )
  
  print(
    panel_c
  )
  
  grDevices::dev.off()
  
} else {
  
  grDevices::tiff(
    filename = tiff_file,
    width = 5.2,
    height = 7.0,
    units = "in",
    res = 600,
    type = "quartz",
    compression = "lzw",
    bg = "white"
  )
  
  print(
    panel_c
  )
  
  grDevices::dev.off()
}


# ============================================================
# 19. Completion summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Final 7KB Figure 4C completed.\n")
cat("============================================================\n")

cat("\n")
cat("Final state-specific genus counts:\n")

print(
  main_cluster_counts
)

cat("\n")
cat(
  "Total genera in final heatmap: ",
  length(
    genus_order
  ),
  "\n",
  sep = ""
)

cat("\n")
cat("C1 retained genera:\n")

print(
  selected_c1[["Genus"]]
)

cat("\n")
cat("C1 requested but not retained after 7KB qualification:\n")

if (length(missing_c1_genera) == 0) {
  
  cat("None\n")
  
} else {
  
  print(
    missing_c1_genera
  )
}

cat("\n")
cat("C2 retained genera:\n")

print(
  selected_c2[["Genus"]]
)

cat("\n")
cat("C3 retained genera:\n")

print(
  selected_c3[["Genus"]]
)

cat("\n")
cat("Output directory:\n")
cat(
  output_dir,
  "\n"
)