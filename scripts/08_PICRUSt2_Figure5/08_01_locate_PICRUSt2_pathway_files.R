


############################################################
## 08_01_locate_PICRUSt2_pathway_files.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
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
  library(tibble)
  library(stringr)
})

############################################################
## Paths
############################################################

PROJECT_DIR <- PROJECT_ROOT

META_FILE <- file.path(
  OUTPUT_ROOT,
  "00_clean_data",
  "progression127",
  "metadata_7KB_progression127.tsv"
)

OUT_DIR <- file.path(
  OUTPUT_ROOT,
  "08_PICRUSt2_pathways",
  "PICRUSt2_pathway_file_check_7KB"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_CANDIDATES <- file.path(
  OUT_DIR,
  "candidate_PICRUSt2_pathway_files.tsv"
)

OUT_SAMPLE_MATCH <- file.path(
  OUT_DIR,
  "PICRUSt2_pathway_file_sample_match_progression127.tsv"
)

############################################################
## Helper functions
############################################################

read_pathway_header <- function(file) {
  con <- if (str_detect(file, "\\.gz$")) {
    gzfile(file, open = "rt")
  } else {
    file(file, open = "rt")
  }
  
  on.exit(close(con), add = TRUE)
  
  header <- readLines(con, n = 1, warn = FALSE)
  str_split(header, "\t", simplify = TRUE) %>% as.character()
}

guess_file_type <- function(path) {
  x <- basename(path)
  
  case_when(
    str_detect(x, "path_abun_unstrat_descrip") ~ "PICRUSt2 MetaCyc pathway abundance with descriptions",
    str_detect(x, "path_abun_unstrat") ~ "PICRUSt2 MetaCyc pathway abundance",
    str_detect(x, "pred_metagenome_unstrat") ~ "PICRUSt2 predicted metagenome",
    str_detect(x, "pathways") ~ "pathway-related",
    TRUE ~ "unknown"
  )
}

############################################################
## Read metadata
############################################################

if (!file.exists(META_FILE)) {
  stop("Missing metadata file:\n", META_FILE)
}

meta_df <- read_tsv(META_FILE, show_col_types = FALSE)

GROUP_ORDER <- c("Polyp", "UC_remission", "UC_active", "Dysplasia", "CA")

meta_use <- meta_df %>%
  filter(Progression5 %in% GROUP_ORDER) %>%
  mutate(Progression5 = factor(Progression5, levels = GROUP_ORDER))

progression127_samples <- meta_use$SampleID

cat("\nCurrent progression127 sample counts:\n")
print(meta_use %>% count(Progression5, name = "n"), n = Inf, width = Inf)

############################################################
## Find candidate files
############################################################

all_files <- list.files(
  PROJECT_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)

candidate_files <- all_files[
  str_detect(
    basename(all_files),
    "path_abun_unstrat|pathways_out|pred_metagenome_unstrat|metacyc|MetaCyc|pathway"
  )
]

candidate_files <- candidate_files[
  str_detect(candidate_files, "\\.(tsv|tsv.gz|txt|txt.gz|csv|csv.gz)$")
]

candidate_tbl <- tibble(
  File = candidate_files,
  File_name = basename(candidate_files),
  File_type_guess = guess_file_type(candidate_files),
  Size_MB = round(file.info(candidate_files)$size / 1024 / 1024, 3)
) %>%
  arrange(desc(Size_MB), File)

write_tsv(candidate_tbl, OUT_CANDIDATES)

cat("\nCandidate pathway-related files:\n")
print(candidate_tbl, n = Inf, width = Inf)

############################################################
## Check sample matching for candidate tables
############################################################

match_tbl <- bind_rows(
  lapply(candidate_files, function(f) {
    header <- tryCatch(
      read_pathway_header(f),
      error = function(e) character(0)
    )
    
    if (length(header) == 0) {
      return(tibble(
        File = f,
        File_name = basename(f),
        First_column = NA_character_,
        N_columns = NA_integer_,
        N_sample_columns = NA_integer_,
        N_progression127_matched = NA_integer_,
        N_progression127_missing = NA_integer_,
        N_extra_sample_columns = NA_integer_,
        Match_fraction_progression127 = NA_real_,
        Example_matched_samples = NA_character_,
        Example_missing_samples = NA_character_
      ))
    }
    
    first_col <- header[1]
    sample_cols <- header[-1]
    
    matched <- intersect(progression127_samples, sample_cols)
    missing <- setdiff(progression127_samples, sample_cols)
    extra <- setdiff(sample_cols, progression127_samples)
    
    tibble(
      File = f,
      File_name = basename(f),
      First_column = first_col,
      N_columns = length(header),
      N_sample_columns = length(sample_cols),
      N_progression127_matched = length(matched),
      N_progression127_missing = length(missing),
      N_extra_sample_columns = length(extra),
      Match_fraction_progression127 = length(matched) / length(progression127_samples),
      Example_matched_samples = paste(head(matched, 10), collapse = ", "),
      Example_missing_samples = paste(head(missing, 10), collapse = ", ")
    )
  })
) %>%
  arrange(desc(N_progression127_matched), N_progression127_missing, File_name)

write_tsv(match_tbl, OUT_SAMPLE_MATCH)

cat("\nSample matching check:\n")
print(match_tbl, n = Inf, width = Inf)

cat("\nFiles written:\n")
cat("Candidate files: ", OUT_CANDIDATES, "\n")
cat("Sample matching: ", OUT_SAMPLE_MATCH, "\n")

cat("\nDone.\n")