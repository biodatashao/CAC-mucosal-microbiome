


############################################################
## 02_01_alpha_beta_statistics.R
##
## Module 02 - Alpha and beta diversity (Figure 2, Supplementary Figure 1)
##
## Final 7-KB alpha / beta statistics
##
## Inputs:
##   1. progression127 7-KB clean dataset
##   2. paired CA23 vs nonCA23 7-KB clean dataset
##
## Outputs support:
##   - Figure 2
##   - Supplementary Figure 1
##   - Supplementary Table S4
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
  library(tibble)
  library(vegan)
})


############################################################
## 0. Paths
############################################################

PROJECT_DIR <- PROJECT_ROOT


FINAL_ROOT <- file.path(
  PROJECT_DIR,
  "output",
  "analysis"
)


PROGRESSION_DIR <- file.path(
  FINAL_ROOT,
  "00_clean_data",
  "progression127"
)


PAIRED_DIR <- file.path(
  FINAL_ROOT,
  "00_clean_data",
  "CA23_nonCA23_paired"
)


OUT_DIR <- file.path(
  FINAL_ROOT,
  "01_Figure2_alpha_beta",
  "statistics"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
## 1. Input files
############################################################

PROG_COUNT_FILE <- file.path(
  PROGRESSION_DIR,
  "asv_count_7KB_progression127.tsv"
)


PROG_META_FILE <- file.path(
  PROGRESSION_DIR,
  "metadata_7KB_progression127.tsv"
)


PAIRED_COUNT_FILE <- file.path(
  PAIRED_DIR,
  "asv_count_7KB_CA23_nonCA23_paired.tsv"
)


PAIRED_META_FILE <- file.path(
  PAIRED_DIR,
  "metadata_7KB_CA23_nonCA23_paired.tsv"
)


required_files <- c(
  PROG_COUNT_FILE,
  PROG_META_FILE,
  PAIRED_COUNT_FILE,
  PAIRED_META_FILE
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
## 2. Analysis settings
############################################################

GROUP_ORDER <- c(
  "Polyp",
  "UC_remission",
  "UC_active",
  "Dysplasia",
  "CA"
)


SEED_MAIN <- 20260628

PROGRESSION_PERMUTATIONS <- 999

PAIRED_PERMUTATIONS <- 9999


############################################################
## 3. Helper: read count matrix
############################################################

read_count_matrix <- function(path) {
  
  df <- read_tsv(
    path,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  
  
  names(df) <- trimws(
    names(df)
  )
  
  
  asv_col <- names(df)[1]
  
  
  mat <- df %>%
    column_to_rownames(
      asv_col
    ) %>%
    as.matrix()
  
  
  storage.mode(
    mat
  ) <- "numeric"
  
  
  if (anyNA(mat)) {
    
    stop(
      paste0(
        "NA detected in count matrix:\n",
        path
      )
    )
  }
  
  
  if (any(mat < 0)) {
    
    stop(
      paste0(
        "Negative counts detected:\n",
        path
      )
    )
  }
  
  
  mat
}


############################################################
## 4. Helper: alpha diversity
############################################################

calculate_alpha <- function(sample_by_asv) {
  
  tibble(
    SampleID = rownames(
      sample_by_asv
    ),
    
    Sequencing_depth = rowSums(
      sample_by_asv
    ),
    
    Observed_ASVs = rowSums(
      sample_by_asv > 0
    ),
    
    Shannon = vegan::diversity(
      sample_by_asv,
      index = "shannon"
    ),
    
    Simpson = vegan::diversity(
      sample_by_asv,
      index = "simpson"
    )
  )
}


############################################################
## 5. Helper: extract adonis2 result
############################################################

extract_adonis <- function(x, label) {
  
  df <- as.data.frame(
    x
  )
  
  
  tibble(
    Analysis = label,
    
    Df = df$Df[1],
    
    SumOfSqs = df$SumOfSqs[1],
    
    R2 = df$R2[1],
    
    F = df$F[1],
    
    p_value = df$`Pr(>F)`[1]
  )
}


############################################################
## 6. Read progression127 count table
############################################################

prog_asv_by_sample <- read_count_matrix(
  PROG_COUNT_FILE
)


prog_meta <- read_tsv(
  PROG_META_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


required_prog_cols <- c(
  "SampleID",
  "Progression5"
)


missing_prog_cols <- setdiff(
  required_prog_cols,
  names(prog_meta)
)


if (length(missing_prog_cols) > 0) {
  
  stop(
    paste0(
      "Missing progression metadata columns: ",
      paste(
        missing_prog_cols,
        collapse = ", "
      )
    )
  )
}


prog_meta <- prog_meta %>%
  mutate(
    SampleID = as.character(
      SampleID
    ),
    
    Progression5 = factor(
      Progression5,
      levels = GROUP_ORDER
    )
  )


############################################################
## 7. Verify progression cohort
############################################################

if (nrow(prog_meta) != 127) {
  
  stop(
    paste0(
      "Expected 127 progression samples; found ",
      nrow(prog_meta)
    )
  )
}


if (anyNA(
  prog_meta$Progression5
)) {
  
  stop(
    "NA detected in Progression5."
  )
}


missing_prog_samples <- setdiff(
  prog_meta$SampleID,
  colnames(prog_asv_by_sample)
)


if (length(missing_prog_samples) > 0) {
  
  stop(
    paste0(
      "Progression samples missing from count table:\n",
      paste(
        missing_prog_samples,
        collapse = ", "
      )
    )
  )
}


prog_asv_by_sample <- prog_asv_by_sample[
  ,
  prog_meta$SampleID,
  drop = FALSE
]


stopifnot(
  identical(
    colnames(prog_asv_by_sample),
    prog_meta$SampleID
  )
)


prog_sample_by_asv <- t(
  prog_asv_by_sample
)


############################################################
## 8. Progression alpha diversity
############################################################

prog_alpha <- calculate_alpha(
  prog_sample_by_asv
) %>%
  left_join(
    prog_meta %>%
      select(
        SampleID,
        Progression5
      ),
    by = "SampleID"
  )


write_csv(
  prog_alpha,
  file.path(
    OUT_DIR,
    "Table_01_progression127_alpha_by_sample.csv"
  )
)


############################################################
## 9. Progression alpha summary
############################################################

prog_alpha_summary <- prog_alpha %>%
  group_by(
    Progression5
  ) %>%
  summarise(
    n = n(),
    
    Median_depth = median(
      Sequencing_depth
    ),
    
    Q1_depth = quantile(
      Sequencing_depth,
      0.25
    ),
    
    Q3_depth = quantile(
      Sequencing_depth,
      0.75
    ),
    
    Median_Observed = median(
      Observed_ASVs
    ),
    
    Q1_Observed = quantile(
      Observed_ASVs,
      0.25
    ),
    
    Q3_Observed = quantile(
      Observed_ASVs,
      0.75
    ),
    
    Median_Shannon = median(
      Shannon
    ),
    
    Q1_Shannon = quantile(
      Shannon,
      0.25
    ),
    
    Q3_Shannon = quantile(
      Shannon,
      0.75
    ),
    
    Median_Simpson = median(
      Simpson
    ),
    
    Q1_Simpson = quantile(
      Simpson,
      0.25
    ),
    
    Q3_Simpson = quantile(
      Simpson,
      0.75
    ),
    
    .groups = "drop"
  )


write_csv(
  prog_alpha_summary,
  file.path(
    OUT_DIR,
    "Table_02_progression127_alpha_summary.csv"
  )
)


############################################################
## 10. Progression global alpha tests
##
## Raw P:
##   primary reporting
##
## q_value_BH:
##   BH adjustment across the 3 predefined alpha metrics
############################################################

kw_observed <- kruskal.test(
  Observed_ASVs ~ Progression5,
  data = prog_alpha
)


kw_shannon <- kruskal.test(
  Shannon ~ Progression5,
  data = prog_alpha
)


kw_simpson <- kruskal.test(
  Simpson ~ Progression5,
  data = prog_alpha
)


prog_kw <- tibble(
  Metric = c(
    "Observed_ASVs",
    "Shannon",
    "Simpson"
  ),
  
  statistic = c(
    unname(
      kw_observed$statistic
    ),
    
    unname(
      kw_shannon$statistic
    ),
    
    unname(
      kw_simpson$statistic
    )
  ),
  
  df = c(
    unname(
      kw_observed$parameter
    ),
    
    unname(
      kw_shannon$parameter
    ),
    
    unname(
      kw_simpson$parameter
    )
  ),
  
  p_value_raw = c(
    kw_observed$p.value,
    kw_shannon$p.value,
    kw_simpson$p.value
  )
) %>%
  mutate(
    q_value_BH = p.adjust(
      p_value_raw,
      method = "BH"
    )
  )


write_csv(
  prog_kw,
  file.path(
    OUT_DIR,
    "Table_03_progression127_alpha_Kruskal_Wallis.csv"
  )
)


############################################################
## 11. Pairwise alpha diversity
##
## 10 pairwise Wilcoxon tests per metric.
##
## BH adjustment is performed separately within each metric.
############################################################

pairwise_wilcox_one_metric <- function(
    df,
    metric_name
) {
  
  pairs <- combn(
    GROUP_ORDER,
    2,
    simplify = FALSE
  )
  
  
  out <- lapply(
    pairs,
    function(pair) {
      
      group1 <- pair[1]
      
      group2 <- pair[2]
      
      
      x <- df %>%
        filter(
          Progression5 == group1
        ) %>%
        pull(
          all_of(
            metric_name
          )
        )
      
      
      y <- df %>%
        filter(
          Progression5 == group2
        ) %>%
        pull(
          all_of(
            metric_name
          )
        )
      
      
      test <- wilcox.test(
        x,
        y,
        paired = FALSE,
        exact = FALSE
      )
      
      
      tibble(
        Metric = metric_name,
        
        Group1 = group1,
        
        Group2 = group2,
        
        n_Group1 = length(x),
        
        n_Group2 = length(y),
        
        W = unname(
          test$statistic
        ),
        
        p_value_raw = test$p.value
      )
    }
  )
  
  
  bind_rows(
    out
  ) %>%
    mutate(
      q_value_BH = p.adjust(
        p_value_raw,
        method = "BH"
      )
    )
}


prog_pairwise_alpha <- bind_rows(
  
  pairwise_wilcox_one_metric(
    prog_alpha,
    "Observed_ASVs"
  ),
  
  pairwise_wilcox_one_metric(
    prog_alpha,
    "Shannon"
  ),
  
  pairwise_wilcox_one_metric(
    prog_alpha,
    "Simpson"
  )
)


write_csv(
  prog_pairwise_alpha,
  file.path(
    OUT_DIR,
    "Table_04_progression127_pairwise_Wilcoxon_BH.csv"
  )
)


############################################################
## 12. Progression relative abundance
##
## Same logic as old script:
## sample total normalization followed by Bray-Curtis.
############################################################

prog_depth <- rowSums(
  prog_sample_by_asv
)


if (any(
  prog_depth <= 0
)) {
  
  stop(
    "Zero-depth progression sample detected."
  )
}


prog_rel <- sweep(
  prog_sample_by_asv,
  1,
  prog_depth,
  "/"
)


############################################################
## 13. Progression Bray-Curtis
############################################################

prog_bray <- vegan::vegdist(
  prog_rel,
  method = "bray"
)


############################################################
## 14. Progression PCoA
##
## IMPORTANT:
## Match old script exactly:
##
## cmdscale(
##   bray_dist,
##   eig = TRUE,
##   k = 2
## )
##
## NO add = TRUE
############################################################

prog_pcoa <- cmdscale(
  prog_bray,
  eig = TRUE,
  k = 2
)


prog_eig <- prog_pcoa$eig


prog_positive_eig <- prog_eig[
  prog_eig > 0
]


prog_variance <- 100 *
  prog_eig[
    1:2
  ] /
  sum(
    prog_positive_eig
  )


prog_pcoa_df <- tibble(
  SampleID = rownames(
    prog_pcoa$points
  ),
  
  PCoA1 = prog_pcoa$points[
    ,
    1
  ],
  
  PCoA2 = prog_pcoa$points[
    ,
    2
  ]
) %>%
  left_join(
    prog_meta %>%
      select(
        SampleID,
        Progression5
      ),
    by = "SampleID"
  )


write_csv(
  prog_pcoa_df,
  file.path(
    OUT_DIR,
    "Table_05_progression127_Bray_PCoA_coordinates.csv"
  )
)


prog_pcoa_variance <- tibble(
  Axis = c(
    "PCoA1",
    "PCoA2"
  ),
  
  Variance_percent = prog_variance
)


write_csv(
  prog_pcoa_variance,
  file.path(
    OUT_DIR,
    "Table_06_progression127_Bray_PCoA_variance.csv"
  )
)


############################################################
## 15. Progression global PERMANOVA
############################################################

set.seed(
  SEED_MAIN
)


prog_permanova <- vegan::adonis2(
  prog_bray ~ Progression5,
  data = prog_meta,
  permutations = PROGRESSION_PERMUTATIONS
)


prog_permanova_summary <- extract_adonis(
  prog_permanova,
  "Progression5_global"
)


write_csv(
  prog_permanova_summary,
  file.path(
    OUT_DIR,
    "Table_07_progression127_PERMANOVA.csv"
  )
)


############################################################
## 16. Progression Betadisper
############################################################

prog_bd <- vegan::betadisper(
  prog_bray,
  group = prog_meta$Progression5
)


set.seed(
  SEED_MAIN
)


prog_bd_perm <- vegan::permutest(
  prog_bd,
  permutations = PROGRESSION_PERMUTATIONS
)


prog_bd_tab <- as.data.frame(
  prog_bd_perm$tab
)


prog_bd_summary <- tibble(
  Df = prog_bd_tab$Df[1],
  
  Sum_Sq = prog_bd_tab$`Sum Sq`[1],
  
  F = prog_bd_tab$F[1],
  
  p_value = prog_bd_tab$`Pr(>F)`[1]
)


write_csv(
  prog_bd_summary,
  file.path(
    OUT_DIR,
    "Table_08_progression127_Betadisper.csv"
  )
)


############################################################
## 17. Pairwise PERMANOVA
##
## 10 comparisons.
## 999 permutations per comparison.
## BH correction across the 10 raw P values.
############################################################

pairwise_permanova <- function(
    relative_matrix,
    metadata
) {
  
  pairs <- combn(
    GROUP_ORDER,
    2,
    simplify = FALSE
  )
  
  
  out <- lapply(
    seq_along(
      pairs
    ),
    function(i) {
      
      pair <- pairs[[i]]
      
      group1 <- pair[1]
      
      group2 <- pair[2]
      
      
      meta_sub <- metadata %>%
        filter(
          Progression5 %in%
            c(
              group1,
              group2
            )
        ) %>%
        droplevels()
      
      
      rel_sub <- relative_matrix[
        meta_sub$SampleID,
        ,
        drop = FALSE
      ]
      
      
      bray_sub <- vegan::vegdist(
        rel_sub,
        method = "bray"
      )
      
      
      set.seed(
        SEED_MAIN + i
      )
      
      
      perm_sub <- vegan::adonis2(
        bray_sub ~ Progression5,
        data = meta_sub,
        permutations = PROGRESSION_PERMUTATIONS
      )
      
      
      perm_df <- as.data.frame(
        perm_sub
      )
      
      
      tibble(
        Comparison = paste0(
          group1,
          " vs ",
          group2
        ),
        
        Group1 = group1,
        
        Group2 = group2,
        
        n_Group1 = sum(
          meta_sub$Progression5 ==
            group1
        ),
        
        n_Group2 = sum(
          meta_sub$Progression5 ==
            group2
        ),
        
        Pseudo_F = perm_df$F[1],
        
        R2 = perm_df$R2[1],
        
        p_value_raw =
          perm_df$`Pr(>F)`[1]
      )
    }
  )
  
  
  bind_rows(
    out
  ) %>%
    mutate(
      q_value_BH = p.adjust(
        p_value_raw,
        method = "BH"
      )
    )
}


prog_pairwise_perm <- pairwise_permanova(
  prog_rel,
  prog_meta
)


write_csv(
  prog_pairwise_perm,
  file.path(
    OUT_DIR,
    "Table_09_progression127_pairwise_PERMANOVA_BH.csv"
  )
)


############################################################
## 18. Read paired CA23 vs nonCA23 dataset
############################################################

paired_asv_by_sample <- read_count_matrix(
  PAIRED_COUNT_FILE
)


paired_meta <- read_tsv(
  PAIRED_META_FILE,
  show_col_types = FALSE,
  progress = FALSE,
  name_repair = "minimal"
)


required_paired_cols <- c(
  "SampleID",
  "PairID",
  "Group"
)


missing_paired_cols <- setdiff(
  required_paired_cols,
  names(paired_meta)
)


if (length(missing_paired_cols) > 0) {
  
  stop(
    paste0(
      "Missing paired metadata columns: ",
      paste(
        missing_paired_cols,
        collapse = ", "
      )
    )
  )
}


paired_meta <- paired_meta %>%
  mutate(
    SampleID = as.character(
      SampleID
    ),
    
    PairID = as.character(
      PairID
    ),
    
    Group = factor(
      Group,
      levels = c(
        "nonCAC",
        "CAC"
      )
    )
  )


############################################################
## 19. Verify paired cohort
############################################################

if (nrow(paired_meta) != 46) {
  
  stop(
    paste0(
      "Expected 46 paired samples; found ",
      nrow(paired_meta)
    )
  )
}


if (n_distinct(
  paired_meta$PairID
) != 23) {
  
  stop(
    "Expected 23 PairIDs."
  )
}


pair_check <- paired_meta %>%
  count(
    PairID
  )


if (!all(
  pair_check$n == 2
)) {
  
  stop(
    "At least one pair is incomplete."
  )
}


missing_paired_samples <- setdiff(
  paired_meta$SampleID,
  colnames(paired_asv_by_sample)
)


if (length(missing_paired_samples) > 0) {
  
  stop(
    paste0(
      "Paired samples missing from count table:\n",
      paste(
        missing_paired_samples,
        collapse = ", "
      )
    )
  )
}


paired_asv_by_sample <- paired_asv_by_sample[
  ,
  paired_meta$SampleID,
  drop = FALSE
]


stopifnot(
  identical(
    colnames(paired_asv_by_sample),
    paired_meta$SampleID
  )
)


paired_sample_by_asv <- t(
  paired_asv_by_sample
)


############################################################
## 20. Paired alpha diversity
############################################################

paired_alpha <- calculate_alpha(
  paired_sample_by_asv
) %>%
  left_join(
    paired_meta %>%
      select(
        SampleID,
        PairID,
        Group
      ),
    by = "SampleID"
  )


write_csv(
  paired_alpha,
  file.path(
    OUT_DIR,
    "Table_10_paired_CA23_nonCA23_alpha_by_sample.csv"
  )
)


############################################################
## 21. Paired alpha summary
############################################################

paired_alpha_summary <- paired_alpha %>%
  group_by(
    Group
  ) %>%
  summarise(
    n = n(),
    
    Median_depth = median(
      Sequencing_depth
    ),
    
    Q1_depth = quantile(
      Sequencing_depth,
      0.25
    ),
    
    Q3_depth = quantile(
      Sequencing_depth,
      0.75
    ),
    
    Median_Observed = median(
      Observed_ASVs
    ),
    
    Q1_Observed = quantile(
      Observed_ASVs,
      0.25
    ),
    
    Q3_Observed = quantile(
      Observed_ASVs,
      0.75
    ),
    
    Median_Shannon = median(
      Shannon
    ),
    
    Q1_Shannon = quantile(
      Shannon,
      0.25
    ),
    
    Q3_Shannon = quantile(
      Shannon,
      0.75
    ),
    
    Median_Simpson = median(
      Simpson
    ),
    
    Q1_Simpson = quantile(
      Simpson,
      0.25
    ),
    
    Q3_Simpson = quantile(
      Simpson,
      0.75
    ),
    
    .groups = "drop"
  )


write_csv(
  paired_alpha_summary,
  file.path(
    OUT_DIR,
    "Table_11_paired_CA23_nonCA23_alpha_summary.csv"
  )
)


############################################################
## 22. Paired Wilcoxon tests
############################################################

paired_wilcox_metric <- function(
    df,
    metric_name
) {
  
  wide <- df %>%
    select(
      PairID,
      Group,
      all_of(
        metric_name
      )
    ) %>%
    pivot_wider(
      names_from = Group,
      values_from = all_of(
        metric_name
      )
    ) %>%
    arrange(
      PairID
    )
  
  
  if (
    anyNA(wide$nonCAC) ||
    anyNA(wide$CAC)
  ) {
    
    stop(
      paste0(
        "Incomplete paired data for ",
        metric_name
      )
    )
  }
  
  
  test <- wilcox.test(
    wide$CAC,
    wide$nonCAC,
    paired = TRUE,
    exact = FALSE
  )
  
  
  tibble(
    Metric = metric_name,
    
    N_pairs = nrow(
      wide
    ),
    
    V = unname(
      test$statistic
    ),
    
    p_value_raw =
      test$p.value,
    
    Median_nonCAC = median(
      wide$nonCAC
    ),
    
    Median_CAC = median(
      wide$CAC
    ),
    
    Median_pair_difference_CAC_minus_nonCAC =
      median(
        wide$CAC -
          wide$nonCAC
      )
  )
}


paired_alpha_tests <- bind_rows(
  
  paired_wilcox_metric(
    paired_alpha,
    "Observed_ASVs"
  ),
  
  paired_wilcox_metric(
    paired_alpha,
    "Shannon"
  ),
  
  paired_wilcox_metric(
    paired_alpha,
    "Simpson"
  )
)


write_csv(
  paired_alpha_tests,
  file.path(
    OUT_DIR,
    "Table_12_paired_CA23_nonCA23_Wilcoxon.csv"
  )
)


############################################################
## 23. Paired relative abundance
############################################################

paired_depth <- rowSums(
  paired_sample_by_asv
)


if (any(
  paired_depth <= 0
)) {
  
  stop(
    "Zero-depth paired sample detected."
  )
}


paired_rel <- sweep(
  paired_sample_by_asv,
  1,
  paired_depth,
  "/"
)


############################################################
## 24. Paired Bray-Curtis
############################################################

paired_bray <- vegan::vegdist(
  paired_rel,
  method = "bray"
)


############################################################
## 25. Paired PCoA
##
## Match old script:
##
## cmdscale(
##   bray_dist,
##   eig = TRUE,
##   k = 2
## )
##
## NO add = TRUE
############################################################

paired_pcoa <- cmdscale(
  paired_bray,
  eig = TRUE,
  k = 2
)


paired_eig <- paired_pcoa$eig


paired_positive_eig <- paired_eig[
  paired_eig > 0
]


paired_variance <- 100 *
  paired_eig[
    1:2
  ] /
  sum(
    paired_positive_eig
  )


paired_pcoa_df <- tibble(
  SampleID = rownames(
    paired_pcoa$points
  ),
  
  PCoA1 = paired_pcoa$points[
    ,
    1
  ],
  
  PCoA2 = paired_pcoa$points[
    ,
    2
  ]
) %>%
  left_join(
    paired_meta %>%
      select(
        SampleID,
        PairID,
        Group
      ),
    by = "SampleID"
  )


write_csv(
  paired_pcoa_df,
  file.path(
    OUT_DIR,
    "Table_13_paired_CA23_nonCA23_Bray_PCoA_coordinates.csv"
  )
)


paired_pcoa_variance <- tibble(
  Axis = c(
    "PCoA1",
    "PCoA2"
  ),
  
  Variance_percent =
    paired_variance
)


write_csv(
  paired_pcoa_variance,
  file.path(
    OUT_DIR,
    "Table_14_paired_CA23_nonCA23_Bray_PCoA_variance.csv"
  )
)


############################################################
## 26. Paired PERMANOVA
##
## Match old paired script:
##   permutations = 9999
##   strata = PairID
############################################################

set.seed(
  SEED_MAIN
)


paired_permanova <- vegan::adonis2(
  paired_bray ~ Group,
  data = paired_meta,
  permutations = PAIRED_PERMUTATIONS,
  strata = paired_meta$PairID
)


paired_permanova_summary <- extract_adonis(
  paired_permanova,
  "Paired_CAC_vs_nonCAC"
)


write_csv(
  paired_permanova_summary,
  file.path(
    OUT_DIR,
    "Table_15_paired_CA23_nonCA23_PERMANOVA.csv"
  )
)


############################################################
## 27. Paired Betadisper
############################################################

paired_bd <- vegan::betadisper(
  paired_bray,
  group = paired_meta$Group
)


set.seed(
  SEED_MAIN
)


paired_bd_perm <- vegan::permutest(
  paired_bd,
  permutations = PROGRESSION_PERMUTATIONS
)


paired_bd_tab <- as.data.frame(
  paired_bd_perm$tab
)


paired_bd_summary <- tibble(
  Df = paired_bd_tab$Df[1],
  
  Sum_Sq = paired_bd_tab$`Sum Sq`[1],
  
  F = paired_bd_tab$F[1],
  
  p_value =
    paired_bd_tab$`Pr(>F)`[1]
)


write_csv(
  paired_bd_summary,
  file.path(
    OUT_DIR,
    "Table_16_paired_CA23_nonCA23_Betadisper.csv"
  )
)


############################################################
## 28. Main statistics summary
############################################################

summary_main <- tibble(
  Analysis = c(
    "Progression Observed Kruskal-Wallis",
    "Progression Shannon Kruskal-Wallis",
    "Progression Simpson Kruskal-Wallis",
    "Progression Bray PERMANOVA",
    "Progression Bray Betadisper",
    "Paired CAC vs nonCAC Observed Wilcoxon",
    "Paired CAC vs nonCAC Shannon Wilcoxon",
    "Paired CAC vs nonCAC Simpson Wilcoxon",
    "Paired CAC vs nonCAC Bray PERMANOVA",
    "Paired CAC vs nonCAC Bray Betadisper"
  ),
  
  Effect = c(
    NA_real_,
    NA_real_,
    NA_real_,
    
    prog_permanova_summary$R2,
    
    prog_bd_summary$F,
    
    NA_real_,
    NA_real_,
    NA_real_,
    
    paired_permanova_summary$R2,
    
    paired_bd_summary$F
  ),
  
  p_value_raw = c(
    
    prog_kw$p_value_raw[
      prog_kw$Metric ==
        "Observed_ASVs"
    ],
    
    prog_kw$p_value_raw[
      prog_kw$Metric ==
        "Shannon"
    ],
    
    prog_kw$p_value_raw[
      prog_kw$Metric ==
        "Simpson"
    ],
    
    prog_permanova_summary$p_value,
    
    prog_bd_summary$p_value,
    
    paired_alpha_tests$p_value_raw[
      paired_alpha_tests$Metric ==
        "Observed_ASVs"
    ],
    
    paired_alpha_tests$p_value_raw[
      paired_alpha_tests$Metric ==
        "Shannon"
    ],
    
    paired_alpha_tests$p_value_raw[
      paired_alpha_tests$Metric ==
        "Simpson"
    ],
    
    paired_permanova_summary$p_value,
    
    paired_bd_summary$p_value
  )
)


write_csv(
  summary_main,
  file.path(
    OUT_DIR,
    "SUMMARY_Figure2_and_SuppFigure1_statistics.csv"
  )
)


############################################################
## 29. Console report
############################################################

cat("\n")
cat("============================================================\n")
cat("7-KB ALPHA / BETA STATISTICS COMPLETE\n")
cat("============================================================\n\n")


cat("PROGRESSION127 GROUP COUNTS:\n")

print(
  table(
    prog_meta$Progression5
  )
)


cat("\n")
cat("============================================================\n")
cat("PROGRESSION ALPHA DIVERSITY\n")
cat("============================================================\n\n")


cat("PROGRESSION ALPHA SUMMARY:\n")

print(
  as.data.frame(
    prog_alpha_summary
  ),
  row.names = FALSE
)


cat("\nPROGRESSION GLOBAL ALPHA TESTS:\n")
cat("Raw P + BH q across the 3 predefined alpha metrics\n\n")


print(
  as.data.frame(
    prog_kw
  ),
  row.names = FALSE
)


cat("\nPROGRESSION PAIRWISE ALPHA TESTS:\n")
cat("BH adjustment separately within each metric\n\n")


print(
  as.data.frame(
    prog_pairwise_alpha
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("PROGRESSION BETA DIVERSITY\n")
cat("============================================================\n\n")


cat("PROGRESSION BRAY PERMANOVA:\n")

print(
  as.data.frame(
    prog_permanova_summary
  ),
  row.names = FALSE
)


cat("\nPROGRESSION BETADISPER:\n")

print(
  as.data.frame(
    prog_bd_summary
  ),
  row.names = FALSE
)


cat("\nPROGRESSION PCoA VARIANCE:\n")

print(
  as.data.frame(
    prog_pcoa_variance
  ),
  row.names = FALSE
)


cat("\nPAIRWISE PERMANOVA:\n")

print(
  as.data.frame(
    prog_pairwise_perm
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("PAIRED CA23 vs nonCA23\n")
cat("============================================================\n\n")


cat("PAIRED ALPHA SUMMARY:\n")

print(
  as.data.frame(
    paired_alpha_summary
  ),
  row.names = FALSE
)


cat("\nPAIRED WILCOXON TESTS:\n")

print(
  as.data.frame(
    paired_alpha_tests
  ),
  row.names = FALSE
)


cat("\nPAIRED BRAY PERMANOVA:\n")

print(
  as.data.frame(
    paired_permanova_summary
  ),
  row.names = FALSE
)


cat("\nPAIRED BETADISPER:\n")

print(
  as.data.frame(
    paired_bd_summary
  ),
  row.names = FALSE
)


cat("\nPAIRED PCoA VARIANCE:\n")

print(
  as.data.frame(
    paired_pcoa_variance
  ),
  row.names = FALSE
)


cat("\n")
cat("============================================================\n")
cat("MAIN SUMMARY\n")
cat("============================================================\n\n")


print(
  as.data.frame(
    summary_main
  ),
  row.names = FALSE
)


cat("\nOutput directory:\n")

cat(
  OUT_DIR,
  "\n"
)


cat("\nDone.\n")