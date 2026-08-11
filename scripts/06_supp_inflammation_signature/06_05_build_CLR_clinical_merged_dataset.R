


############################################################
## 06_05_build_CLR_clinical_merged_dataset.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## 7KB progression127
##
## Oral-associated microbial score and UCG-005 clinical association
############################################################


rm(list = ls())

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
# 1. Packages
############################################################

library(tidyverse)
library(readxl)
library(compositions)



############################################################
# 2. Paths
############################################################


asv_file <- file.path(PROJECT_ROOT, "output/analysis/00_clean_data/progression127/asv_count_7KB_progression127.tsv")


taxonomy_file <- file.path(PROJECT_ROOT, "output/analysis/00_clean_data/progression127/taxonomy_7KB_progression127.tsv")


metadata_file <- file.path(PROJECT_ROOT, "output/analysis/00_clean_data/progression127/metadata_7KB_progression127.tsv")


clinical_file <- file.path(PROJECT_ROOT, "data/clinical/FFPE_patient_information.xlsx")



output_dir <- file.path(PROJECT_ROOT, "output/analysis/07_clinical_microbiome_association")


dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################
# 3. Read files
############################################################


asv <- read.delim(
  asv_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


taxonomy <- read.delim(
  taxonomy_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


metadata <- read.delim(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


clinical <- read_excel(
  clinical_file
)

clinical <- clinical %>%
  mutate(
    CRP = as.numeric(CRP),
    ESR = as.numeric(ESR),
    CEA = as.numeric(CEA),
    CA199 = as.numeric(CA199)
  )

clinical <- clinical %>%
  mutate(
    CRP = as.numeric(trimws(as.character(CRP))),
    ESR = as.numeric(trimws(as.character(ESR))),
    CEA = as.numeric(trimws(as.character(CEA))),
    CA199 = as.numeric(trimws(as.character(CA199)))
  )


############################################################
# 4. Extract genus from taxonomy
############################################################


taxonomy <- taxonomy %>%
  mutate(
    Genus = str_extract(
      Taxonomy,
      "g__[^;]+"
    )
  ) %>%
  mutate(
    Genus = str_remove(
      Genus,
      "g__"
    )
  )



print(head(taxonomy))



############################################################
# 5. ASV table -> genus abundance
############################################################


asv_long <- asv %>%
  pivot_longer(
    cols = -ASV,
    names_to = "SampleID",
    values_to = "Count"
  )


asv_long <- asv_long %>%
  left_join(
    taxonomy %>%
      select(
        ASV,
        Genus
      ),
    by = "ASV"
  )



asv_long <- asv_long %>%
  filter(
    !is.na(Genus),
    Genus != "",
    !grepl(
      "uncultured|unclassified",
      Genus,
      ignore.case = TRUE
    )
  )



genus_count <- asv_long %>%
  group_by(
    SampleID,
    Genus
  ) %>%
  summarise(
    Count = sum(Count),
    .groups = "drop"
  )



genus_wide <- genus_count %>%
  pivot_wider(
    names_from = Genus,
    values_from = Count,
    values_fill = 0
  )



############################################################
# 6. Relative abundance
############################################################


genus_matrix <- genus_wide %>%
  column_to_rownames(
    "SampleID"
  )


genus_relative <- genus_matrix /
  rowSums(genus_matrix)



############################################################
# 7. CLR transformation
############################################################


genus_clr <- clr(
  as.matrix(
    genus_relative + 1e-06
  )
)


genus_clr <- as.data.frame(
  genus_clr
)


genus_clr$SampleID <- rownames(genus_clr)



############################################################
# 8. Merge metadata and clinical
############################################################


metadata_clinical <- metadata %>%
  left_join(
    clinical,
    by = c(
      "SampleID" = "Sample_ID"
    )
  )

clinical <- clinical %>%
  mutate(
    CRP = as.numeric(
      gsub(
        "[^0-9\\.]+",
        "",
        CRP
      )
    ),
    ESR = as.numeric(
      gsub(
        "[^0-9\\.]+",
        "",
        ESR
      )
    ),
    CEA = as.numeric(
      gsub(
        "[^0-9\\.]+",
        "",
        CEA
      )
    ),
    CA199 = as.numeric(
      gsub(
        "[^0-9\\.]+",
        "",
        CA199
      )
    )
  )


analysis_df <- genus_clr %>%
  inner_join(
    metadata_clinical,
    by = "SampleID"
  )



write.csv(
  analysis_df,
  file.path(
    output_dir,
    "CLR_clinical_merged_dataset.csv"
  ),
  row.names = FALSE
)



############################################################
# 9. Check data
############################################################


cat(
  "Total samples:",
  nrow(analysis_df),
  "\n"
)


cat(
  "Groups:\n"
)

print(
  table(
    analysis_df$Group_main
  )
)



############################################################
# 10. Spearman function
############################################################


spearman_test <- function(
    data,
    x,
    y
){
  
  tmp <- data %>%
    select(
      all_of(
        c(x,y)
      )
    ) %>%
    drop_na()
  
  
  if(nrow(tmp)<5){
    
    return(
      tibble(
        N=nrow(tmp),
        rho=NA,
        p=NA
      )
    )
    
  }
  
  
  res <- cor.test(
    tmp[[x]],
    tmp[[y]],
    method="spearman",
    exact = FALSE
  )
  
  
  tibble(
    N=nrow(tmp),
    rho=as.numeric(res$estimate),
    p=res$p.value
  )
  
}



############################################################
# Analysis 1
# Oral-associated microbial score
############################################################


uc_df <- analysis_df %>%
  filter(
    Group_main %in%
      c(
        "UC_remission",
        "UC_active"
      )
  )



oral_taxa <- c(
  "Olsenella",
  "Granulicatella",
  "Aerococcus",
  "Peptostreptococcus",
  "Lancefieldella"
)



cat("Missing oral taxa:\n")

print(
  oral_taxa[
    !oral_taxa %in% colnames(uc_df)
  ]
)



############################################################
# Calculate score
############################################################


oral_matrix <- uc_df %>%
  select(
    all_of(oral_taxa)
  )


oral_z <- as.data.frame(
  scale(
    oral_matrix
  )
)


uc_df$oral_score <- rowMeans(
  oral_z,
  na.rm = TRUE
)



############################################################
# Score vs CRP/ESR
############################################################


oral_score_result <- bind_rows(
  
  spearman_test(
    uc_df,
    "oral_score",
    "CRP"
  ) %>%
    mutate(
      Marker="CRP"
    ),
  
  
  spearman_test(
    uc_df,
    "oral_score",
    "ESR"
  ) %>%
    mutate(
      Marker="ESR"
    )
  
)



write.csv(
  oral_score_result,
  file.path(
    output_dir,
    "oral_score_vs_CRP_ESR.csv"
  ),
  row.names = FALSE
)



############################################################
# Individual oral taxa vs CRP
############################################################


oral_taxa_result <- map_dfr(
  oral_taxa,
  function(x){
    
    spearman_test(
      uc_df,
      x,
      "CRP"
    ) %>%
      mutate(
        Genus=x
      )
    
  }
)



oral_taxa_result$q <- p.adjust(
  oral_taxa_result$p,
  method="BH"
)



write.csv(
  oral_taxa_result,
  file.path(
    output_dir,
    "oral_taxa_vs_CRP_BH.csv"
  ),
  row.names = FALSE
)



############################################################
# Analysis 2
# UCG-005
############################################################


ucg005 <- grep(
  "UCG-005",
  colnames(analysis_df),
  value = TRUE
)[1]


cat(
  "UCG005 column:",
  ucg005,
  "\n"
)



############################################################
# UCG005 vs inflammation
############################################################


uc_ucg <- analysis_df %>%
  filter(
    Group_main %in%
      c(
        "UC_remission",
        "UC_active"
      )
  )


ucg_inflammation <- bind_rows(
  
  spearman_test(
    uc_ucg,
    ucg005,
    "CRP"
  ) %>%
    mutate(
      Marker="CRP"
    ),
  
  
  spearman_test(
    uc_ucg,
    ucg005,
    "ESR"
  ) %>%
    mutate(
      Marker="ESR"
    )
  
)



write.csv(
  ucg_inflammation,
  file.path(
    output_dir,
    "UCG005_vs_inflammation_UC.csv"
  ),
  row.names = FALSE
)



############################################################
# UCG005 vs tumor markers
############################################################


ca_df <- analysis_df %>%
  filter(
    Group_main=="CA"
  )



if("CEA" %in% colnames(ca_df)){
  
  
  write.csv(
    
    spearman_test(
      ca_df,
      ucg005,
      "CEA"
    ),
    
    file.path(
      output_dir,
      "UCG005_vs_CEA_CA.csv"
    ),
    
    row.names = FALSE
    
  )
  
}



if("CA199" %in% colnames(ca_df)){
  
  
  write.csv(
    
    spearman_test(
      ca_df,
      ucg005,
      "CA199"
    ),
    
    file.path(
      output_dir,
      "UCG005_vs_CA19_9_CA.csv"
    ),
    
    row.names = FALSE
    
  )
  
}

############################################################
# Analysis 3
# Oral-associated score:
# UC active vs UC remission
############################################################


oral_group_test <- uc_df %>%
  select(
    Group_main,
    oral_score
  ) %>%
  drop_na()



print(
  table(
    oral_group_test$Group_main
  )
)



wilcox_result <- wilcox.test(
  oral_score ~ Group_main,
  data = oral_group_test,
  exact = FALSE
)



print(wilcox_result)



summary_oral_score <- oral_group_test %>%
  group_by(
    Group_main
  ) %>%
  summarise(
    n = n(),
    median = median(
      oral_score
    ),
    IQR = IQR(
      oral_score
    ),
    .groups = "drop"
  )



write.csv(
  summary_oral_score,
  file.path(
    output_dir,
    "oral_score_active_vs_remission_summary.csv"
  ),
  row.names = FALSE
)



write.csv(
  data.frame(
    statistic = wilcox_result$statistic,
    p = wilcox_result$p.value
  ),
  file.path(
    output_dir,
    "oral_score_active_vs_remission_wilcox.csv"
  ),
  row.names = FALSE
)

cat("Finished.\n")