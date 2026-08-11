


############################################################
## 06_01_inflammation3_alpha_beta_statistics.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## Supplementary analysis
## Polyp vs UC remission vs Active UC
## Alpha and Beta diversity
##
## Based on 7KB rerun progression127 pipeline
############################################################

rm(list = ls())

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(readr)


############################################################
# Project paths
############################################################

PROJECT_DIR <- file.path(PROJECT_ROOT, "output/analysis")

DATA_DIR <- file.path(
  PROJECT_DIR,
  "00_clean_data",
  "progression127"
)


OUT_DIR <- file.path(
  PROJECT_DIR,
  "05_Supp_UC_inflammation_alpha_beta_7KB"
)


if (!dir.exists(OUT_DIR)) {
  dir.create(
    OUT_DIR,
    recursive = TRUE
  )
}


############################################################
# Input
############################################################

asv_file <- file.path(
  DATA_DIR,
  "asv_count_7KB_progression127.tsv"
)

meta_file <- file.path(
  DATA_DIR,
  "metadata_7KB_progression127.tsv"
)


############################################################
# Read data
############################################################

asv <- read.delim(
  asv_file,
  row.names = 1,
  check.names = FALSE
)


meta <- read.delim(
  meta_file,
  check.names = FALSE
)


############################################################
# Select inflammation spectrum
############################################################

GROUP_ORDER <- c(
  "Polyp",
  "UC_remission",
  "UC_active"
)


meta <- meta %>%
  filter(
    Progression5 %in% GROUP_ORDER
  )


rownames(meta) <- meta$SampleID


asv <- asv[, rownames(meta)]


############################################################
# Build phyloseq
############################################################

OTU <- otu_table(
  as.matrix(asv),
  taxa_are_rows = TRUE
)


SAM <- sample_data(meta)


ps <- phyloseq(
  OTU,
  SAM
)



############################################################
# Alpha diversity
############################################################

otu <- as.data.frame(
  as(otu_table(ps), "matrix")
)


if (!taxa_are_rows(ps)) {
  otu <- t(otu)
}


alpha_df <- data.frame(
  SampleID = colnames(otu),
  Observed = colSums(otu > 0),
  Shannon = diversity(
    t(otu),
    index = "shannon"
  ),
  Simpson = diversity(
    t(otu),
    index = "simpson"
  )
)


alpha_df <- alpha_df %>%
  left_join(
    meta %>%
      select(
        SampleID,
        Progression5
      ),
    by = "SampleID"
  )


write.csv(
  alpha_df,
  file.path(
    OUT_DIR,
    "Supplement_UC_inflammation_alpha_values.csv"
  ),
  row.names = FALSE
)



############################################################
# Alpha statistics
############################################################

alpha_results <- data.frame()


for (m in c(
  "Observed",
  "Shannon",
  "Simpson"
)) {
  
  
  kw <- kruskal.test(
    alpha_df[[m]] ~ alpha_df$Progression5
  )
  
  
  alpha_results <- rbind(
    alpha_results,
    data.frame(
      Metric = m,
      KW_p = kw$p.value
    )
  )
  
}


write.csv(
  alpha_results,
  file.path(
    OUT_DIR,
    "Supplement_UC_inflammation_alpha_KW.csv"
  ),
  row.names = FALSE
)



############################################################
# Pairwise Wilcoxon
############################################################

pairwise_results <- data.frame()


for (m in c(
  "Observed",
  "Shannon",
  "Simpson"
)) {
  
  
  pw <- pairwise.wilcox.test(
    alpha_df[[m]],
    alpha_df$Progression5,
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  
  tmp <- as.data.frame(
    as.table(
      pw$p.value
    )
  )
  
  
  colnames(tmp) <- c(
    "Group1",
    "Group2",
    "BH_p"
  )
  
  
  tmp$Metric <- m
  
  
  pairwise_results <- rbind(
    pairwise_results,
    tmp
  )
  
}


write.csv(
  pairwise_results,
  file.path(
    OUT_DIR,
    "Supplement_UC_inflammation_pairwise_Wilcoxon.csv"
  ),
  row.names = FALSE
)



############################################################
# Beta diversity
############################################################

bray <- distance(
  ps,
  method = "bray"
)


ordination <- ordinate(
  ps,
  method = "PCoA",
  distance = bray
)


pcoa_df <- as.data.frame(
  ordination$vectors
)


pcoa_df$SampleID <- rownames(
  pcoa_df
)


pcoa_df <- pcoa_df %>%
  left_join(
    meta %>%
      select(
        SampleID,
        Progression5
      ),
    by = "SampleID"
  )


write.csv(
  pcoa_df,
  file.path(
    OUT_DIR,
    "Supplement_UC_inflammation_Bray_PCoA_coordinates.csv"
  ),
  row.names = FALSE
)

############################################################
# PCoA variance
############################################################

pcoa_var <- data.frame(
  Axis = c(
    "PCoA1",
    "PCoA2"
  ),
  Variance_percent = round(
    ordination$values$Relative_eig[1:2] * 100,
    1
  )
)


write.csv(
  pcoa_var,
  file.path(
    OUT_DIR,
    "Supplement_UC_inflammation_Bray_PCoA_variance.csv"
  ),
  row.names = FALSE
)

############################################################
# PERMANOVA
############################################################

permanova <- adonis2(
  bray ~ Progression5,
  data = meta,
  permutations = 999
)


write.csv(
  as.data.frame(permanova),
  file.path(
    OUT_DIR,
    "Supplement_UC_inflammation_PERMANOVA.csv"
  )
)



############################################################
# Betadisper
############################################################

bd <- betadisper(
  bray,
  meta$Progression5
)


bd_test <- anova(
  bd
)


write.csv(
  as.data.frame(bd_test),
  file.path(
    OUT_DIR,
    "Supplement_UC_inflammation_betadisper.csv"
  )
)



cat(
  "Finished UC inflammation alpha/beta analysis\n"
)