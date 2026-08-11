#!/usr/bin/env Rscript



############################################################
## 06_08_assemble_SuppFig_inflammation_signature.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## Assemble the supplementary inflammation-signature figure from saved panels.
############################################################




rm(list = ls())

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


suppressPackageStartupMessages({
  
  library(cowplot)
  
})


PROJECT_DIR <- file.path(PROJECT_ROOT, "output/analysis")


LEFSE_DIR <- file.path(
  PROJECT_DIR,
  "06_LEfSe_UC_inflammation3_7KB_GENUS_ONLY",
  "SupplementaryFigure_inflammation3_LEfSe_clean"
)


CLINICAL_DIR <- file.path(
  PROJECT_DIR,
  "07_clinical_microbiome_association",
  "SuppFig_UC_inflammation_signature"
)



OUT_DIR <- CLINICAL_DIR



############################################################
## Read RDS
############################################################


pA <- readRDS(
  file.path(
    LEFSE_DIR,
    "SupplementaryFig_inflammation3_LEfSe_LDA_barplot_7KB.rds"
  )
)


pB <- readRDS(
  file.path(
    LEFSE_DIR,
    "SupplementaryFig_inflammation3_LEfSe_heatmap_7KB.rds"
  )
)


pC <- readRDS(
  file.path(
    CLINICAL_DIR,
    "PanelC_oral_ESR.rds"
  )
)


pD <- readRDS(
  file.path(
    CLINICAL_DIR,
    "PanelD_oral_CRP.rds"
  )
)


pE <- readRDS(
  file.path(
    CLINICAL_DIR,
    "PanelE_UCG005_CRP.rds"
  )
)


pF <- readRDS(
  file.path(
    CLINICAL_DIR,
    "PanelF_UCG005_ESR.rds"
  )
)



############################################################
## Right columns
############################################################


right_top <- plot_grid(
  
  pC,
  pD,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1
  )
  
)



right_bottom <- plot_grid(
  
  pE,
  pF,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1
  )
  
)



############################################################
## Left column
############################################################


left_column <- plot_grid(
  
  pA,
  
  pB,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1
  )
  
)



############################################################
## Right column
############################################################


right_column <- plot_grid(
  
  right_top,
  
  right_bottom,
  
  ncol = 1,
  
  rel_heights = c(
    1,
    1
  )
  
)



############################################################
## Final
############################################################


final_fig <- plot_grid(
  
  left_column,
  
  right_column,
  
  ncol = 2,
  
  rel_widths = c(
    1.45,
    0.85
  )
  
)



############################################################
## Export
############################################################


pdf_file <- file.path(
  OUT_DIR,
  "Supplementary_Figure_UC_oral_inflammation_signature_v3.pdf"
)


tiff_file <- file.path(
  OUT_DIR,
  "Supplementary_Figure_UC_oral_inflammation_signature_v3.tiff"
)



pdf(
  pdf_file,
  width = 9,
  height = 10,
  family="Helvetica",
  useDingbats = FALSE
)

print(final_fig)

dev.off()



if(requireNamespace("ragg",quietly=TRUE)){
  
  ragg::agg_tiff(
    tiff_file,
    width=9,
    height=10,
    units="in",
    res=600,
    compression="lzw"
  )
  
  print(final_fig)
  
  dev.off()
  
}



message("Finished")

message(pdf_file)

message(tiff_file)