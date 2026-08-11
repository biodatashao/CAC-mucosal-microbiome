#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 07_01_JonckheereTerpstra_trend_top20_genera.R
##
## Module 07 - Sensitivity analyses and negative controls
############################################################



options(stringsAsFactors = FALSE)
options(width = 220)

if (!requireNamespace("clinfun", quietly = TRUE)) {
  stop("Package 'clinfun' is required.")
}

suppressPackageStartupMessages({
  library(data.table)
})

# ==============================================================================
# Paths
# ==============================================================================

IN_FILE <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "02_Figure3_taxa_LEfSe/JT_trend_top20_7KB/",
  "JT_trend_top20_genera_7KB_sample_abundance.tsv"
)

OUT_DIR <- paste0(
  file.path(PROJECT_ROOT, "output/analysis/"),
  "02_Figure3_taxa_LEfSe/JT_trend_top20_7KB"
)

OUT_FILE <- file.path(
  OUT_DIR,
  "JT_trend_top20_genera_7KB_results_FINAL.tsv"
)

if (!file.exists(IN_FILE)) {
  stop(
    paste0(
      "Input file not found:\n",
      IN_FILE
    )
  )
}

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# Read previously generated 7KB sample-level genus abundance
# ==============================================================================

dat <- fread(
  IN_FILE,
  data.table = FALSE
)

required_cols <- c(
  "Genus",
  "SampleID",
  "Relative_abundance_percent",
  "Group"
)

missing_cols <- setdiff(
  required_cols,
  colnames(dat)
)

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}


# ==============================================================================
# Normalize group labels
# ==============================================================================

dat$Group <- as.character(dat$Group)

dat$Group[
  dat$Group == "UC_remission"
] <- "UC remission"

dat$Group[
  dat$Group == "UC_active"
] <- "UC active"

dat$Group[
  dat$Group == "CA"
] <- "CAC"

ordered_groups <- c(
  "UC remission",
  "UC active",
  "Dysplasia",
  "CAC"
)

dat <- dat[
  dat$Group %in% ordered_groups,
  ,
  drop = FALSE
]

dat$Group_order <- match(
  dat$Group,
  ordered_groups
)

dat$Relative_abundance_percent <- as.numeric(
  dat$Relative_abundance_percent
)


# ==============================================================================
# Basic checks
# ==============================================================================

if (any(is.na(dat$Group_order))) {
  stop("Unexpected group label detected.")
}

if (any(!is.finite(dat$Relative_abundance_percent))) {
  stop("Non-finite abundance value detected.")
}

genera <- unique(
  as.character(dat$Genus)
)

cat("\n============================================================\n")
cat("JT AUDIT INPUT\n")
cat("============================================================\n")

for (g in ordered_groups) {
  
  n_samples <- length(
    unique(
      dat$SampleID[
        dat$Group == g
      ]
    )
  )
  
  cat(
    sprintf(
      "%-15s n = %d\n",
      g,
      n_samples
    )
  )
}

cat(
  sprintf(
    "\nNumber of genera tested = %d\n",
    length(genera)
  )
)


# ==============================================================================
# Jonckheere-Terpstra
#
# We test BOTH directions for every genus.
# The direction is then defined by whichever one-sided test is smaller.
#
# This avoids choosing the alternative based on the same data before testing.
# ==============================================================================

results_list <- vector(
  "list",
  length(genera)
)

for (i in seq_along(genera)) {
  
  genus_i <- genera[i]
  
  tmp <- dat[
    dat$Genus == genus_i,
    ,
    drop = FALSE
  ]
  
  tmp <- tmp[
    is.finite(tmp$Relative_abundance_percent) &
      !is.na(tmp$Group_order),
    ,
    drop = FALSE
  ]
  
  jt_inc <- suppressWarnings(
    clinfun::jonckheere.test(
      x = tmp$Relative_abundance_percent,
      g = tmp$Group_order,
      alternative = "increasing"
    )
  )
  
  jt_dec <- suppressWarnings(
    clinfun::jonckheere.test(
      x = tmp$Relative_abundance_percent,
      g = tmp$Group_order,
      alternative = "decreasing"
    )
  )
  
  p_inc <- as.numeric(
    jt_inc$p.value
  )
  
  p_dec <- as.numeric(
    jt_dec$p.value
  )
  
  if (
    !is.finite(p_inc) ||
    !is.finite(p_dec)
  ) {
    stop(
      paste0(
        "Non-finite JT P value for genus: ",
        genus_i
      )
    )
  }
  
  if (p_inc <= p_dec) {
    
    direction <- "increasing"
    p_value <- p_inc
    
  } else {
    
    direction <- "decreasing"
    p_value <- p_dec
  }
  
  group_medians <- rep(
    NA_real_,
    length(ordered_groups)
  )
  
  names(group_medians) <- ordered_groups
  
  for (group_i in ordered_groups) {
    
    x_i <- tmp$Relative_abundance_percent[
      tmp$Group == group_i
    ]
    
    group_medians[group_i] <- median(
      x_i,
      na.rm = TRUE
    )
  }
  
  results_list[[i]] <- data.frame(
    Genus = genus_i,
    Direction = direction,
    P_value = p_value,
    UC_remission_median = group_medians["UC remission"],
    UC_active_median = group_medians["UC active"],
    Dysplasia_median = group_medians["Dysplasia"],
    CAC_median = group_medians["CAC"],
    stringsAsFactors = FALSE
  )
}

results <- do.call(
  rbind,
  results_list
)

results$Q_value_BH <- p.adjust(
  results$P_value,
  method = "BH"
)

results <- results[
  order(
    results$Q_value_BH,
    results$P_value
  ),
  ,
  drop = FALSE
]

rownames(results) <- NULL


# ==============================================================================
# Save
# ==============================================================================

fwrite(
  results,
  OUT_FILE,
  sep = "\t"
)


# ==============================================================================
# Helper for console output
# ==============================================================================

show_genus <- function(target) {
  
  z <- results[
    results$Genus == target,
    ,
    drop = FALSE
  ]
  
  if (nrow(z) != 1) {
    
    cat(
      target,
      ": NOT FOUND\n",
      sep = ""
    )
    
    return(
      invisible(NULL)
    )
  }
  
  cat(
    sprintf(
      "%s\n",
      target
    )
  )
  
  cat(
    sprintf(
      "  medians: %.4f -> %.4f -> %.4f -> %.4f%%\n",
      z$UC_remission_median,
      z$UC_active_median,
      z$Dysplasia_median,
      z$CAC_median
    )
  )
  
  cat(
    sprintf(
      "  direction = %s\n",
      z$Direction
    )
  )
  
  cat(
    sprintf(
      "  P = %.8g\n",
      z$P_value
    )
  )
  
  cat(
    sprintf(
      "  BH q = %.8g\n",
      z$Q_value_BH
    )
  )
}


# ==============================================================================
# Final manuscript audit
# ==============================================================================

cat("\n============================================================\n")
cat("FINAL FIGURE 3 TREND AUDIT\n")
cat("============================================================\n\n")

show_genus(
  "UCG 005"
)

cat("\n")

show_genus(
  "Streptococcus"
)

cat("\n")

show_genus(
  "Bacteroides"
)


cat("\n============================================================\n")
cat("SIGNIFICANT JT TRENDS: BH q < 0.05\n")
cat("============================================================\n")

sig <- results[
  results$Q_value_BH < 0.05,
  ,
  drop = FALSE
]

if (nrow(sig) == 0) {
  
  cat("NONE\n")
  
} else {
  
  for (i in seq_len(nrow(sig))) {
    
    cat(
      sprintf(
        "%-35s %-10s P = %.6g   q = %.6g\n",
        sig$Genus[i],
        sig$Direction[i],
        sig$P_value[i],
        sig$Q_value_BH[i]
      )
    )
  }
}


cat("\n============================================================\n")
cat("OUTPUT\n")
cat("============================================================\n")

cat(
  OUT_FILE,
  "\n"
)

cat("\nDone.\n")