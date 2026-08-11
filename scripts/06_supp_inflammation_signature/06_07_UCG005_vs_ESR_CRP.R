


############################################################
## 06_07_UCG005_vs_ESR_CRP.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## UCG-005 vs CRP/ESR in UC patients
############################################################


rm(list = ls())

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


library(tidyverse)



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
# 2. Read data
############################################################


df <- read.csv(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)



############################################################
# 3. UC samples
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
# 4. UCG-005 column
############################################################


ucg005_col <- "UCG-005"


print(
  ucg005_col %in% colnames(uc_df)
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
      all_of(c(ucg005_col, marker))
    ) %>%
    drop_na()
  
  
  test <- cor.test(
    tmp[[ucg005_col]],
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



crp_result <- cor_result(
  uc_df,
  "CRP"
)


esr_result <- cor_result(
  uc_df,
  "ESR"
)



print(crp_result)

print(esr_result)



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
      all_of(c(ucg005_col, marker))
    ) %>%
    drop_na()
  
  
  ggplot(
    plot_data,
    aes(
      x = .data[[marker]],
      y = .data[[ucg005_col]]
    )
  ) +
    
    geom_point(
      size = 2.2,
      alpha = 0.8
    ) +
    
    geom_smooth(
      method = "lm",
      se = TRUE,
      linewidth = 0.6
    ) +
    
    labs(
      x = marker,
      y = "UCG-005 CLR abundance",
      title = plot_title
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
# 7. Generate plots
############################################################


p_crp <- make_plot(
  uc_df,
  "CRP",
  crp_result,
  "UCG-005 vs CRP"
)



p_esr <- make_plot(
  uc_df,
  "ESR",
  esr_result,
  "UCG-005 vs ESR"
)



############################################################
# 8. Save
############################################################


ggsave(
  file.path(
    output_dir,
    "PanelE1_UCG005_CRP.pdf"
  ),
  p_crp,
  width = 3.5,
  height = 3
)


ggsave(
  file.path(
    output_dir,
    "PanelE2_UCG005_ESR.pdf"
  ),
  p_esr,
  width = 3.5,
  height = 3
)


ggsave(
  file.path(
    output_dir,
    "PanelE1_UCG005_CRP.tiff"
  ),
  p_crp,
  width = 3.5,
  height = 3,
  dpi = 300
)


ggsave(
  file.path(
    output_dir,
    "PanelE2_UCG005_ESR.tiff"
  ),
  p_esr,
  width = 3.5,
  height = 3,
  dpi = 300
)

saveRDS(
  p_crp,
  file.path(
    output_dir,
    "PanelE_UCG005_CRP.rds"
  )
)


saveRDS(
  p_esr,
  file.path(
    output_dir,
    "PanelF_UCG005_ESR.rds"
  )
)