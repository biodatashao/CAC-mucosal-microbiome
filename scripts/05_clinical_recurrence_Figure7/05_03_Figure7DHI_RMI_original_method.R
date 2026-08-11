#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 05_03_Figure7DHI_RMI_original_method.R
##
## Module 05 - Clinical association and recurrence (Figure 7, Supplementary Figure 3)
##
## Recalculate Figure 7D / H / I using the ORIGINAL manuscript methodology.
##
## Figure 7D
##   Genus-level relative abundance
##   -> remove genera that are zero in all CAC23 samples
##   -> Bray-Curtis
##   -> PCoA
##   -> PERMANOVA, 9999 permutations
##   -> betadisper
##   -> AJCC stage-adjusted PERMANOVA
##
## Figure 7H/I
##   log10(relative abundance + 1e-6)
##   -> strict fold-wise standardization
##   -> RMI = z(UCG-005) - z(Lactococcus)
##   -> out-of-fold RMI score for every CAC sample
##   -> Wilcoxon test for recurrence
##   -> ROC directly from out-of-fold RMI scores
##   -> stratified bootstrap 95% CI, 5000 replicates
############################################################

options(stringsAsFactors = FALSE)
options(width = 220)


# ==============================================================================
# 0. Packages
# ==============================================================================

required_packages <- c(
  "data.table",
  "stringr",
  "vegan",
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
      "Missing required package(s): ",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(vegan)
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

progression_dir <- file.path(
  rerun_root,
  "00_clean_data",
  "progression127"
)

asv_count_file <- file.path(
  progression_dir,
  "asv_count_7KB_progression127.tsv"
)

taxonomy_file <- file.path(
  progression_dir,
  "taxonomy_7KB_progression127.tsv"
)

metadata_file <- file.path(
  progression_dir,
  "metadata_7KB_progression127.tsv"
)

clinical_file <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "00_clinical_metadata_audit",
  "CAC23_clinical_metadata_audit_7KB.tsv"
)

output_dir <- file.path(
  rerun_root,
  "04_Figure7_clinical_recurrence",
  "02_DHI_original_method"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

beta_dir <- file.path(
  output_dir,
  "Figure7D_genus_Bray"
)

rmi_dir <- file.path(
  output_dir,
  "Figure7HI_RMI_LOOCV"
)

dir.create(beta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rmi_dir, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 2. Input checks
# ==============================================================================

required_files <- c(
  asv_count_file,
  taxonomy_file,
  metadata_file,
  clinical_file
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
# 3. Helpers
# ==============================================================================

section <- function(title) {
  
  cat("\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat(title, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
}


extract_genus <- function(x) {
  
  x <- as.character(x)
  
  genus <- stringr::str_extract(
    x,
    "g__[^;]+"
  )
  
  genus <- stringr::str_replace(
    genus,
    "^g__",
    ""
  )
  
  genus <- stringr::str_trim(
    genus
  )
  
  genus[
    is.na(genus) |
      genus == ""
  ] <- NA_character_
  
  genus
}


canonicalize_genus <- function(x) {
  
  x <- as.character(x)
  
  x[x == "UCG-005"] <- "UCG_005"
  x[x == "UCG-009"] <- "UCG_009"
  
  x[
    x == "Christensenellaceae_R-7_group"
  ] <- "Christensenellaceae_R_7_group"
  
  x[
    x == "f__Eubacterium__eligens_group"
  ] <- "Eubacterium_eligens_group"
  
  x[
    x == "[Eubacterium]_eligens_group"
  ] <- "Eubacterium_eligens_group"
  
  x[
    x == "[Eubacterium] eligens group"
  ] <- "Eubacterium_eligens_group"
  
  x
}


bad_genus <- function(x) {
  
  x0 <- tolower(
    as.character(x)
  )
  
  is.na(x0) |
    x0 == "" |
    x0 == "na" |
    x0 == "unassigned" |
    x0 == "unclassified" |
    x0 == "uncultured" |
    x0 == "unknown" |
    stringr::str_detect(x0, "uncultured") |
    stringr::str_detect(x0, "unclassified") |
    stringr::str_detect(x0, "unknown")
}


training_z <- function(
    training_values,
    target_values
) {
  
  mu <- mean(
    training_values,
    na.rm = TRUE
  )
  
  sigma <- stats::sd(
    training_values,
    na.rm = TRUE
  )
  
  if (
    !is.finite(sigma) ||
    sigma <= 0
  ) {
    stop(
      "Invalid training SD during LOOCV.",
      call. = FALSE
    )
  }
  
  (
    target_values -
      mu
  ) /
    sigma
}


# ==============================================================================
# 4. Read data
# ==============================================================================

section(
  "READ DATA"
)

asv_table <- data.table::fread(
  asv_count_file,
  data.table = FALSE,
  check.names = FALSE
)

taxonomy <- data.table::fread(
  taxonomy_file,
  data.table = FALSE,
  check.names = FALSE
)

metadata <- data.table::fread(
  metadata_file,
  data.table = FALSE,
  check.names = FALSE
)

clinical <- data.table::fread(
  clinical_file,
  data.table = FALSE,
  check.names = FALSE
)


if (!"ASV" %in% colnames(asv_table)) {
  stop(
    "ASV count table lacks ASV column.",
    call. = FALSE
  )
}

if (
  !all(
    c("ASV", "Taxonomy") %in%
    colnames(taxonomy)
  )
) {
  stop(
    "Taxonomy table must contain ASV and Taxonomy.",
    call. = FALSE
  )
}

if (
  !all(
    c("SampleID", "Progression5") %in%
    colnames(metadata)
  )
) {
  stop(
    "Metadata must contain SampleID and Progression5.",
    call. = FALSE
  )
}

if (
  !all(
    c(
      "SampleID",
      "Recurrence_plot",
      "Tumor_stage_plot"
    ) %in%
    colnames(clinical)
  )
) {
  stop(
    "Clinical audit table lacks required columns.",
    call. = FALSE
  )
}


metadata$SampleID <- as.character(
  metadata$SampleID
)

metadata$Progression5 <- as.character(
  metadata$Progression5
)

clinical$SampleID <- as.character(
  clinical$SampleID
)


cac_samples <- metadata$SampleID[
  metadata$Progression5 %in%
    c("CAC", "CA")
]


if (length(cac_samples) != 23) {
  stop(
    paste0(
      "Expected CAC23; found ",
      length(cac_samples),
      "."
    ),
    call. = FALSE
  )
}


clinical <- clinical[
  match(
    cac_samples,
    clinical$SampleID
  ),
  ,
  drop = FALSE
]


if (
  anyNA(
    clinical$SampleID
  )
) {
  stop(
    "CAC23 clinical alignment failed.",
    call. = FALSE
  )
}


cat(
  "CAC samples: ",
  length(cac_samples),
  "\n",
  sep = ""
)

cat("\nRecurrence distribution:\n")
print(
  table(
    clinical$Recurrence_plot
  )
)

cat("\nTumor stage distribution:\n")
print(
  table(
    clinical$Tumor_stage_plot
  )
)


# ==============================================================================
# 5. Build CAC23 ASV matrix
# ==============================================================================

section(
  "BUILD CAC23 ASV MATRIX"
)

missing_samples <- setdiff(
  cac_samples,
  colnames(asv_table)
)

if (length(missing_samples) > 0) {
  stop(
    paste0(
      "CAC samples absent from count matrix:\n",
      paste(missing_samples, collapse = "\n")
    ),
    call. = FALSE
  )
}


asv_matrix <- as.matrix(
  asv_table[
    ,
    cac_samples,
    drop = FALSE
  ]
)

storage.mode(
  asv_matrix
) <- "numeric"

rownames(
  asv_matrix
) <- as.character(
  asv_table$ASV
)


# ==============================================================================
# 6. Align taxonomy and collapse to genus
# ==============================================================================

section(
  "COLLAPSE TO GENUS"
)

taxonomy_index <- match(
  rownames(asv_matrix),
  as.character(taxonomy$ASV)
)

if (anyNA(taxonomy_index)) {
  stop(
    "Some ASVs have no taxonomy match.",
    call. = FALSE
  )
}


genus <- extract_genus(
  taxonomy$Taxonomy[
    taxonomy_index
  ]
)

genus <- canonicalize_genus(
  genus
)

genus[
  bad_genus(genus)
] <- NA_character_


valid_rows <- !is.na(
  genus
)


genus_count <- rowsum(
  asv_matrix[
    valid_rows,
    ,
    drop = FALSE
  ],
  group = genus[
    valid_rows
  ],
  reorder = FALSE
)


cat(
  "Genus rows before CAC zero-filter: ",
  nrow(genus_count),
  "\n",
  sep = ""
)


# ==============================================================================
# 7. Remove genera that are zero in all CAC23
# ==============================================================================

nonzero_genus <- rowSums(
  genus_count
) > 0


genus_count <- genus_count[
  nonzero_genus,
  ,
  drop = FALSE
]


cat(
  "Genus rows after CAC zero-filter: ",
  nrow(genus_count),
  "\n",
  sep = ""
)


# ==============================================================================
# 8. Genus-level relative abundance
# ==============================================================================

genus_relative <- sweep(
  genus_count,
  2,
  colSums(genus_count),
  "/"
)


if (
  any(
    !is.finite(
      genus_relative
    )
  )
) {
  stop(
    "Non-finite genus relative abundance.",
    call. = FALSE
  )
}


# save genus matrix used for D/H/I

genus_relative_out <- data.frame(
  Genus = rownames(genus_relative),
  genus_relative,
  check.names = FALSE
)


data.table::fwrite(
  genus_relative_out,
  file.path(
    output_dir,
    "CAC23_genus_relative_abundance_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


# ==============================================================================
# 9. Figure 7D: Bray-Curtis
# ==============================================================================

section(
  "FIGURE 7D: GENUS-LEVEL BRAY-CURTIS"
)


sample_genus_matrix <- t(
  genus_relative
)


recurrence <- factor(
  clinical$Recurrence_plot,
  levels = c(
    "No recurrence",
    "Recurrence"
  )
)

stage <- factor(
  clinical$Tumor_stage_plot,
  levels = c(
    "Stage I",
    "Stage II",
    "Stage III"
  )
)


rownames(
  sample_genus_matrix
) <- clinical$SampleID


bray <- vegan::vegdist(
  sample_genus_matrix,
  method = "bray"
)


# ==============================================================================
# 10. PCoA
# ==============================================================================

pcoa <- stats::cmdscale(
  bray,
  eig = TRUE,
  k = 2,
  add = TRUE
)


coords <- as.data.frame(
  pcoa$points
)

colnames(coords) <- c(
  "PCoA1",
  "PCoA2"
)

coords$SampleID <- rownames(
  coords
)


coords$Recurrence_plot <- as.character(
  recurrence[
    match(
      coords$SampleID,
      clinical$SampleID
    )
  ]
)

coords$Tumor_stage_plot <- as.character(
  stage[
    match(
      coords$SampleID,
      clinical$SampleID
    )
  ]
)


positive_eigenvalues <- pcoa$eig[
  pcoa$eig > 0
]


pc1_percent <- 100 *
  pcoa$eig[1] /
  sum(
    positive_eigenvalues
  )


pc2_percent <- 100 *
  pcoa$eig[2] /
  sum(
    positive_eigenvalues
  )


# ==============================================================================
# 11. Primary PERMANOVA
# ==============================================================================

set.seed(
  20260716
)

permanova <- vegan::adonis2(
  bray ~ recurrence,
  permutations = 9999,
  by = "margin"
)


permanova_r2 <- permanova$R2[1]
permanova_p <- permanova$`Pr(>F)`[1]


# ==============================================================================
# 12. Betadisper
# ==============================================================================

dispersion <- vegan::betadisper(
  bray,
  group = recurrence,
  type = "median",
  bias.adjust = TRUE
)


set.seed(
  20260716
)

dispersion_test <- vegan::permutest(
  dispersion,
  permutations = 9999
)


betadisper_p <- dispersion_test$tab$`Pr(>F)`[1]


# ==============================================================================
# 13. AJCC-adjusted PERMANOVA
# ==============================================================================

set.seed(
  20260716
)

permanova_adjusted <- vegan::adonis2(
  bray ~ stage + recurrence,
  permutations = 9999,
  by = "margin"
)


adjusted_table <- as.data.frame(
  permanova_adjusted
)

adjusted_table$Term <- rownames(
  adjusted_table
)

rownames(
  adjusted_table
) <- NULL


# ==============================================================================
# 14. Save Figure 7D results
# ==============================================================================

beta_summary <- data.frame(
  Metric = c(
    "PCoA1_variance_percent",
    "PCoA2_variance_percent",
    "PERMANOVA_R2_recurrence",
    "PERMANOVA_P_recurrence",
    "Betadisper_P",
    "Permutations"
  ),
  
  Value = c(
    pc1_percent,
    pc2_percent,
    permanova_r2,
    permanova_p,
    betadisper_p,
    9999
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  coords,
  file.path(
    beta_dir,
    "Figure7D_genus_Bray_PCoA_coordinates_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  beta_summary,
  file.path(
    beta_dir,
    "Figure7D_genus_Bray_summary_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  adjusted_table,
  file.path(
    beta_dir,
    "Figure7D_genus_Bray_AJCC_adjusted_PERMANOVA_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


cat("\nPrimary genus-level Bray-Curtis:\n")

print(
  beta_summary,
  row.names = FALSE
)


cat("\nAJCC-adjusted PERMANOVA:\n")

print(
  adjusted_table,
  row.names = FALSE
)


# ==============================================================================
# 15. Figure 7H/I: identify target genera
# ==============================================================================

section(
  "FIGURE 7H/I: ORIGINAL RMI METHOD"
)


required_genera <- c(
  "Lactococcus",
  "UCG_005"
)


missing_genera <- setdiff(
  required_genera,
  rownames(
    genus_relative
  )
)


if (length(missing_genera) > 0) {
  stop(
    paste0(
      "Missing required genera:\n",
      paste(missing_genera, collapse = "\n")
    ),
    call. = FALSE
  )
}


lactococcus <- as.numeric(
  genus_relative[
    "Lactococcus",
    clinical$SampleID
  ]
)


ucg005 <- as.numeric(
  genus_relative[
    "UCG_005",
    clinical$SampleID
  ]
)


# ==============================================================================
# 16. Original log transform
# ==============================================================================

pseudocount <- 1e-6


lactococcus_log <- log10(
  lactococcus +
    pseudocount
)


ucg005_log <- log10(
  ucg005 +
    pseudocount
)


outcome <- as.integer(
  recurrence ==
    "Recurrence"
)


# ==============================================================================
# 17. Strict LOOCV RMI
#
# No logistic regression.
#
# Each held-out sample receives an RMI calculated using ONLY training mean/SD.
# ==============================================================================

n <- length(
  outcome
)


loocv_rmi <- rep(
  NA_real_,
  n
)


for (i in seq_len(n)) {
  
  training_index <- setdiff(
    seq_len(n),
    i
  )
  
  
  lacto_z_test <- training_z(
    training_values = lactococcus_log[
      training_index
    ],
    target_values = lactococcus_log[i]
  )
  
  
  ucg_z_test <- training_z(
    training_values = ucg005_log[
      training_index
    ],
    target_values = ucg005_log[i]
  )
  
  
  loocv_rmi[i] <-
    ucg_z_test -
    lacto_z_test
}


if (
  any(
    !is.finite(
      loocv_rmi
    )
  )
) {
  stop(
    "Invalid out-of-fold RMI score.",
    call. = FALSE
  )
}


# ==============================================================================
# 18. Figure 7H recurrence comparison
# ==============================================================================

rmi_source <- data.frame(
  SampleID = clinical$SampleID,
  Recurrence_plot = recurrence,
  Recurrence_binary = outcome,
  Lactococcus_relative_abundance = lactococcus,
  UCG_005_relative_abundance = ucg005,
  Lactococcus_log10 = lactococcus_log,
  UCG_005_log10 = ucg005_log,
  RMI_two_genus = loocv_rmi,
  stringsAsFactors = FALSE
)


rmi_wilcox <- stats::wilcox.test(
  RMI_two_genus ~ Recurrence_plot,
  data = rmi_source,
  exact = FALSE
)


rmi_p <- rmi_wilcox$p.value


# ==============================================================================
# 19. Figure 7I ROC
#
# Direct ROC on out-of-fold RMI scores.
# ==============================================================================

roc_object <- pROC::roc(
  response = outcome,
  predictor = loocv_rmi,
  levels = c(
    0,
    1
  ),
  direction = "<",
  quiet = TRUE
)


auc_value <- as.numeric(
  pROC::auc(
    roc_object
  )
)


# ==============================================================================
# 20. Original bootstrap CI
#
# Stratified bootstrap, 5000 replicates.
# ==============================================================================

set.seed(
  20260716
)


auc_ci <- pROC::ci.auc(
  roc_object,
  method = "bootstrap",
  boot.n = 5000,
  boot.stratified = TRUE
)


auc_ci <- as.numeric(
  auc_ci
)


# ==============================================================================
# 21. ROC coordinates
# ==============================================================================

roc_coordinates <- pROC::coords(
  roc_object,
  x = "all",
  ret = c(
    "specificity",
    "sensitivity",
    "threshold"
  ),
  transpose = FALSE
)


roc_coordinates <- as.data.frame(
  roc_coordinates
)

roc_coordinates$One_minus_specificity <-
  1 -
  roc_coordinates$specificity


# ==============================================================================
# 22. Save Figure 7H/I
# ==============================================================================

rmi_summary <- data.frame(
  Metric = c(
    "RMI_Wilcoxon_P",
    "LOOCV_AUC",
    "LOOCV_AUC_CI_lower",
    "LOOCV_AUC_CI_median",
    "LOOCV_AUC_CI_upper",
    "Bootstrap_replicates",
    "Pseudocount"
  ),
  
  Value = c(
    rmi_p,
    auc_value,
    auc_ci[1],
    auc_ci[2],
    auc_ci[3],
    5000,
    pseudocount
  ),
  
  stringsAsFactors = FALSE
)


data.table::fwrite(
  rmi_source,
  file.path(
    rmi_dir,
    "Figure7HI_original_LOOCV_RMI_source_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  rmi_summary,
  file.path(
    rmi_dir,
    "Figure7HI_original_LOOCV_RMI_statistics_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


data.table::fwrite(
  roc_coordinates,
  file.path(
    rmi_dir,
    "Figure7I_original_LOOCV_ROC_coordinates_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


saveRDS(
  roc_object,
  file.path(
    rmi_dir,
    "Figure7I_original_LOOCV_ROC_7KB.rds"
  )
)


cat("\nFigure 7H/I results:\n")

print(
  rmi_summary,
  row.names = FALSE
)


# ==============================================================================
# 23. Final summary
# ==============================================================================

section(
  "FINAL D/H/I SUMMARY"
)


final_summary <- data.frame(
  Result = c(
    "Figure7D_genus_PERMANOVA_R2",
    "Figure7D_genus_PERMANOVA_P",
    "Figure7D_genus_Betadisper_P",
    "Figure7D_PCoA1_percent",
    "Figure7D_PCoA2_percent",
    "Figure7H_RMI_Wilcoxon_P",
    "Figure7I_LOOCV_AUC",
    "Figure7I_AUC_CI_lower",
    "Figure7I_AUC_CI_upper"
  ),
  
  Value = c(
    permanova_r2,
    permanova_p,
    betadisper_p,
    pc1_percent,
    pc2_percent,
    rmi_p,
    auc_value,
    auc_ci[1],
    auc_ci[3]
  ),
  
  stringsAsFactors = FALSE
)


print(
  final_summary,
  row.names = FALSE
)


data.table::fwrite(
  final_summary,
  file.path(
    output_dir,
    "Figure7_DHI_original_method_summary_7KB.tsv"
  ),
  sep = "\t",
  quote = FALSE
)


cat("\nOutput directory:\n")
cat(
  output_dir,
  "\n"
)