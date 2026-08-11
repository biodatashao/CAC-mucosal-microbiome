


############################################################
## 06_03_LEfSe_inflammation3_genus.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## Purpose:
##   Genus-level LEfSe analysis for inflammation spectrum
##
##   Polyp
##       vs
##   UC_remission
##       vs
##   UC_active
##
## Dataset:
##   FFPE 7KB rerun progression127
############################################################
# 0. Clean environment
############################

rm(list = ls())

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################
# 1. Paths
############################

PROJECT_DIR <- PROJECT_ROOT

RUN_DIR <- file.path(
  PROJECT_DIR,
  "output/analysis"
)


DATA_DIR <- file.path(
  RUN_DIR,
  "00_clean_data/progression127"
)


OUTPUT_DIR <- file.path(
  RUN_DIR,
  "06_LEfSe_UC_inflammation3_7KB_GENUS_ONLY"
)


dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# LEfSe environment
############################

CONDA_SH <- "/opt/miniconda3/etc/profile.d/conda.sh"

CONDA_ENV <- "lefse"


LEFSE_INPUT_FILE <- file.path(
  OUTPUT_DIR,
  "LEfSe_input_inflammation3_7KB_GENUS_ONLY.tsv"
)


LEFSE_FORMATTED_FILE <- file.path(
  OUTPUT_DIR,
  "LEfSe_formatted_inflammation3_7KB_GENUS_ONLY.in"
)


LEFSE_RESULT_FILE <- file.path(
  OUTPUT_DIR,
  "LEfSe_results_inflammation3_7KB_GENUS_ONLY.res"
)


############################
# 2. Parameters
############################

GROUP_ORDER <- c(
  "Polyp",
  "UC_remission",
  "UC_active"
)


EXPECTED_N <- c(
  Polyp = 26,
  UC_remission = 36,
  UC_active = 25
)


MIN_PREVALENCE_N <- 5

MIN_MEAN_RA <- 1e-4

LDA_CUTOFF <- 2.0


############################
# 3. Input files
############################

ASV_FILE <- file.path(
  DATA_DIR,
  "asv_count_7KB_progression127.tsv"
)


TAX_FILE <- file.path(
  DATA_DIR,
  "taxonomy_7KB_progression127.tsv"
)


META_FILE <- file.path(
  DATA_DIR,
  "metadata_7KB_progression127.tsv"
)


############################
# 4. Read data
############################

asv <- read.table(
  ASV_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)


taxonomy <- read.table(
  TAX_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


metadata <- read.table(
  META_FILE,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)



############################
# 5. Select inflammation spectrum groups
############################

if (!"Progression5" %in% colnames(metadata)) {

  stop(
    "metadata does not contain Progression5 column"
  )

}


metadata_sub <- metadata[
  metadata$Progression5 %in% GROUP_ORDER,
  ,
  drop = FALSE
]


metadata_sub$Group <- factor(
  metadata_sub$Progression5,
  levels = GROUP_ORDER
)


cat("\nSelected groups:\n")
print(table(metadata_sub$Group))


for (g in GROUP_ORDER) {

  n_now <- sum(metadata_sub$Group == g)

  if (n_now != EXPECTED_N[g]) {

    stop(
      paste0(
        "Sample number mismatch for ",
        g,
        ": expected ",
        EXPECTED_N[g],
        ", got ",
        n_now
      )
    )

  }

}



############################
# 6. Subset ASV table
############################

sample_ids <- rownames(metadata_sub)


asv <- asv[, sample_ids]



############################
# 7. Strict genus extraction
############################

tax_string <- taxonomy$Taxonomy


# Extract genus level from taxonomy string
genus <- sapply(
  strsplit(
    tax_string,
    ";"
  ),
  function(x) {
    
    g <- x[grepl(
      "^g__",
      x
    )]
    
    if (length(g) == 0) {
      return(NA)
    }
    
    return(g)
    
  }
)


genus_clean <- gsub(
  "^g__",
  "",
  genus
)


bad_pattern <- paste(
  c(
    "uncultured",
    "unclassified",
    "norank",
    "unknown",
    "ambiguous"
  ),
  collapse = "|"
)


keep_genus <- !grepl(
  bad_pattern,
  genus_clean,
  ignore.case = TRUE
)


keep_genus <- keep_genus &
  !is.na(genus_clean) &
  genus_clean != ""


asv_genus <- asv[keep_genus, ]

genus_clean <- genus_clean[keep_genus]



############################
# 8. Aggregate ASV to genus
############################

genus_table <- rowsum(
  as.matrix(asv_genus),
  group = genus_clean
)



############################
# 9. Relative abundance
############################

genus_RA <- sweep(
  genus_table,
  2,
  colSums(genus_table),
  FUN = "/"
)



############################
# 10. Filtering
############################

prevalence <- rowSums(
  genus_RA > 0
)


mean_RA <- rowMeans(
  genus_RA
)


keep_filter <- (
  prevalence >= MIN_PREVALENCE_N
) &
(
  mean_RA >= MIN_MEAN_RA
)


genus_RA_filtered <- genus_RA[
  keep_filter,
  ,
  drop = FALSE
]


cat(
  "\nNumber of genera after filtering:",
  nrow(genus_RA_filtered),
  "\n"
)



############################
# 11. Generate LEfSe input
############################

lefse_input <- rbind(
  Class = as.character(metadata_sub$Group),
  genus_RA_filtered
)


LEFSE_INPUT_FILE <- file.path(
  OUTPUT_DIR,
  "LEfSe_input_inflammation3_7KB_GENUS_ONLY.tsv"
)


write.table(
  lefse_input,
  file = LEFSE_INPUT_FILE,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)



############################
# 12. Generate LEfSe input
############################

lefse_input <- rbind(
  Class = as.character(metadata_sub$Group),
  genus_RA_filtered
)


write.table(
  lefse_input,
  file = LEFSE_INPUT_FILE,
  sep = "\t",
  quote = FALSE,
  col.names = FALSE
)



############################################################
## Run LEfSe
##
## EXACTLY same conda strategy as original script
############################################################

cat(
  "\nRunning LEfSe inside R with clean conda environment...\n"
)


unlink(
  c(
    LEFSE_FORMATTED_FILE,
    LEFSE_RESULT_FILE
  )
)


run_cmd <- paste(
  "unset R_HOME R_LIBS R_LIBS_USER R_LIBS_SITE R_PROFILE R_ENVIRON &&",
  "source", shQuote(CONDA_SH), "&&",
  "conda activate", CONDA_ENV, "&&",
  "export R_HOME=$CONDA_PREFIX/lib/R &&",
  "export R_DEFAULT_PACKAGES='utils,stats,graphics,grDevices,datasets,methods' &&",
  "echo CONDA_PREFIX=$CONDA_PREFIX &&",
  "echo R_HOME=$R_HOME &&",
  "which R &&",
  "R --version | head -n 1 &&",
  "which lefse_format_input.py &&",
  "which lefse_run.py &&",
  "lefse_format_input.py",
  shQuote(
    LEFSE_INPUT_FILE
  ),
  shQuote(
    LEFSE_FORMATTED_FILE
  ),
  "-c 1 -o 1000000 &&",
  "lefse_run.py",
  shQuote(
    LEFSE_FORMATTED_FILE
  ),
  shQuote(
    LEFSE_RESULT_FILE
  ),
  "-l",
  as.character(
    LDA_CUTOFF
  )
)


lefse_log <- system2(
  command = "bash",
  args = c(
    "-lc",
    shQuote(
      run_cmd
    )
  ),
  stdout = TRUE,
  stderr = TRUE
)


cat(
  "\nLEfSe log:\n"
)


cat(
  paste(
    lefse_log,
    collapse = "\n"
  )
)



############################
# 13. Check LEfSe output
############################

if (!file.exists(LEFSE_RESULT_FILE)) {
  
  stop(
    "LEfSe result file was not generated"
  )
  
}


############################
# 14. Clean marker table
############################

marker <- read.table(
  LEFSE_RESULT_FILE,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  fill = TRUE
)


colnames(marker)[1:4] <- c(
  "Class",
  "Feature",
  "LDA",
  "p_value"
)


marker_clean <- marker[
  marker$LDA >= LDA_CUTOFF,
  ,
  drop = FALSE
]


CLEAN_FILE <- file.path(
  OUTPUT_DIR,
  "LEfSe_clean_markers_inflammation3_7KB_GENUS_ONLY.csv"
)


write.csv(
  marker_clean,
  CLEAN_FILE,
  row.names = FALSE
)



############################
# 15. Summary
############################

cat("\n==============================\n")
cat("Inflammation3 LEfSe completed\n")
cat("==============================\n")


cat("\nSample numbers:\n")
print(table(metadata_sub$Group))


cat("\nFiltered genera:")
cat(nrow(genus_RA_filtered), "\n")


cat("\nMarkers:")
cat(nrow(marker_clean), "\n")


cat("\nOutput:\n")
cat(OUTPUT_DIR, "\n")


cat("\nMarker file:\n")
cat(CLEAN_FILE, "\n")


############################################################
# END
############################################################