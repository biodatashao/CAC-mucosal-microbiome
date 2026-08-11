


############################################################
## 06_06_oral_score_vs_ESR_CRP.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## Panel C/D:
## Oral-associated microbial score vs ESR/CRP
############################################################


rm(list = ls())

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


library(tidyverse)
library(ggplot2)



############################################################
# 1. Path
############################################################


input_file <- file.path(PROJECT_ROOT, "output/analysis/07_clinical_microbiome_association/CLR_clinical_merged_dataset.csv")


output_dir <- file.path(PROJECT_ROOT, "output/analysis/07_clinical_microbiome_association/SuppFig_UC_inflammation_signature")



dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################
# 2. Read
############################################################


df <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)



############################################################
# 3. UC only
############################################################


uc_df <- df %>%
  filter(
    Group_main %in%
      c(
        "UC_active",
        "UC_remission"
      )
  )



############################################################
# 4. Oral-associated score
############################################################


oral_taxa <- c(
  "Olsenella",
  "Granulicatella",
  "Aerococcus",
  "Peptostreptococcus",
  "Lancefieldella"
)



print(
  oral_taxa[
    !oral_taxa %in% colnames(uc_df)
  ]
)



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
# 5. Correlation function
############################################################


cor_result <- function(
    data,
    marker
){
  
  tmp <- data %>%
    select(
      oral_score,
      all_of(marker)
    ) %>%
    drop_na()
  
  
  test <- cor.test(
    tmp$oral_score,
    tmp[[marker]],
    method = "spearman",
    exact = FALSE
  )
  
  
  list(
    n = nrow(tmp),
    rho = round(
      as.numeric(test$estimate),
      3
    ),
    p = signif(
      test$p.value,
      3
    )
  )
  
}



esr_result <- cor_result(
  uc_df,
  "ESR"
)



crp_result <- cor_result(
  uc_df,
  "CRP"
)



print(esr_result)
print(crp_result)



############################################################
# 6. Plot function
############################################################


make_plot <- function(
    data,
    marker,
    result,
    plot_title
){
  
  plot_data <- data %>%
    select(
      oral_score,
      all_of(marker)
    ) %>%
    drop_na()
  
  
  ggplot(
    plot_data,
    aes(
      x = .data[[marker]],
      y = oral_score
    )
  ) +
    
    geom_point(
      size = 2,
      alpha = 0.75
    ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.6
    ) +
    
    theme_classic(
      base_size = 11
    ) +
    
    theme(
      plot.title = element_text(
        size = 12,
        face = "bold",
        hjust = 0.5
      )
    ) +
    
    theme(
      axis.title = element_text(
        size = 12
      ),
      axis.text = element_text(
        size = 10
      ),
      plot.title = element_blank()
    ) +
    
      labs(
        x = marker,
        y = "Oral-associated microbial score",
        title = plot_title
      )+
    
    annotate(
      "text",
      x = Inf,
      y = Inf,
      hjust = 1.1,
      vjust = 1.5,
      label =
        paste0(
          "Spearman r = ",
          result$rho,
          "\n",
          "P = ",
          result$p
        ),
      size = 3.8
    )
  
}



############################################################
# 7. Generate Panel C
############################################################


p_esr <- make_plot(
  uc_df,
  "ESR",
  esr_result,
  "Oral-associated score vs ESR"
)



############################################################
# 8. Generate Panel D
############################################################


p_crp <- make_plot(
  uc_df,
  "CRP",
  crp_result,
  "D. Oral-associated score vs CRP"
)



############################################################
# 9. Save
############################################################


ggsave(
  file.path(
    output_dir,
    "PanelC_oral_score_ESR.pdf"
  ),
  p_esr,
  width = 3.5,
  height = 3
)


ggsave(
  file.path(
    output_dir,
    "PanelD_oral_score_CRP.pdf"
  ),
  p_crp,
  width = 4.5,
  height = 4
)



ggsave(
  file.path(
    output_dir,
    "PanelC_oral_score_ESR.png"
  ),
  p_esr,
  width = 4.5,
  height = 4,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "PanelD_oral_score_CRP.png"
  ),
  p_crp,
  width = 4.5,
  height = 4,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "PanelC_oral_score_ESR.tiff"
  ),
  p_esr,
  width = 4.5,
  height = 4,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "PanelD_oral_score_CRP.tiff"
  ),
  p_crp,
  width = 4.5,
  height = 4,
  dpi = 300
)

saveRDS(
  p_esr,
  file.path(
    output_dir,
    "PanelC_oral_ESR.rds"
  )
)


saveRDS(
  p_crp,
  file.path(
    output_dir,
    "PanelD_oral_CRP.rds"
  )
)

############################################################
# END
############################################################