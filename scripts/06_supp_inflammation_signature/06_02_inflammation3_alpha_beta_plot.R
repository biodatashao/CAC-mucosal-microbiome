#!/usr/bin/env Rscript




############################################################
## 06_02_inflammation3_alpha_beta_plot.R
##
## Module 06 - Inflammation-spectrum supplementary analysis
##
## Supplementary Figure
## Polyp vs UC remission vs UC active
##
## Match Figure 2 submission style
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
  library(ggplot2)
  library(cowplot)
  library(grid)
  
})


############################################################
## 1. Paths
############################################################

PROJECT_DIR <- file.path(PROJECT_ROOT, "output/analysis")


STAT_DIR <- file.path(
  PROJECT_DIR,
  "05_Supp_UC_inflammation_alpha_beta_7KB"
)


OUT_DIR <- file.path(
  STAT_DIR,
  "figures"
)


dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################
## 2. Input files
############################################################


f_alpha <- file.path(
  STAT_DIR,
  "Supplement_UC_inflammation_alpha_values.csv"
)


f_kw <- file.path(
  STAT_DIR,
  "Supplement_UC_inflammation_alpha_KW.csv"
)


f_pcoa <- file.path(
  STAT_DIR,
  "Supplement_UC_inflammation_Bray_PCoA_coordinates.csv"
)


f_perm <- file.path(
  STAT_DIR,
  "Supplement_UC_inflammation_PERMANOVA.csv"
)



############################################################
## 3. Figure 2 style
############################################################


group_levels <- c(
  "Polyp",
  "UC remission",
  "UC active"
)



group_colors <- c(
  
  "Polyp" = "#7DA9B7",
  
  "UC remission" = "#8DBA91",
  
  "UC active" = "#E3A04F"
  
)



font_family <- "Helvetica"


base_size <- 8



fig_theme <- theme_classic(
  base_size = base_size,
  base_family = font_family
) +
  
  theme(
    
    text = element_text(
      color = "black"
    ),
    
    axis.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    axis.text = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.45
    ),
    
    axis.ticks = element_line(
      linewidth = 0.45
    ),
    
    plot.title = element_text(
      size = 8.5,
      hjust = 0.5
    ),
    
    legend.title = element_blank(),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.margin = margin(
      5,
      5,
      5,
      5
    )
    
  )



############################################################
## 4. Export function
############################################################


export_plot <- function(
    plot,
    filename,
    width,
    height,
    dpi = 600
){
  
  pdf_file <- file.path(
    OUT_DIR,
    paste0(filename, ".pdf")
  )
  
  
  png_file <- file.path(
    OUT_DIR,
    paste0(filename, ".png")
  )
  
  
  tif_file <- file.path(
    OUT_DIR,
    paste0(filename, ".tiff")
  )
  
  
  grDevices::pdf(
    pdf_file,
    width = width,
    height = height,
    family = font_family,
    useDingbats = FALSE
  )
  
  print(plot)
  
  grDevices::dev.off()
  
  
  
  if(
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ){
    
    ragg::agg_png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    
    
    ragg::agg_tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    
  } else {
    
    
    grDevices::png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    
    
    grDevices::tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
  }
  
}



############################################################
## 5. Read data
############################################################


alpha <- read_csv(
  f_alpha,
  show_col_types = FALSE
)



alpha <- alpha %>%
  
  mutate(
    
    Group_display = recode(
      Progression5,
      "Polyp" = "Polyp",
      "UC_remission" = "UC remission",
      "UC_active" = "UC active"
    ),
    
    Group_display = factor(
      Group_display,
      levels = group_levels
    )
    
  )



kw <- read_csv(
  f_kw,
  show_col_types = FALSE
)



pcoa <- read_csv(
  f_pcoa,
  show_col_types = FALSE
)



pcoa <- pcoa %>%
  
  rename(
    
    PCoA1 = Axis.1,
    
    PCoA2 = Axis.2
    
  ) %>%
  
  mutate(
    
    Group_display = recode(
      Progression5,
      "Polyp" = "Polyp",
      "UC_remission" = "UC remission",
      "UC_active" = "UC active"
    ),
    
    Group_display = factor(
      Group_display,
      levels = group_levels
    )
    
  )



perm <- read_csv(
  f_perm,
  show_col_types = FALSE
)



############################################################
## 6. Labels
############################################################


get_kw <- function(metric){
  
  kw %>%
    
    filter(
      Metric == metric
    ) %>%
    
    pull(
      KW_p
    )
  
}



format_p <- function(x){
  
  if(
    x < 0.001
  ){
    
    return(
      "P < 0.001"
    )
    
  }
  
  paste0(
    "P = ",
    formatC(
      x,
      format = "f",
      digits = 3
    )
  )
  
}



label_shannon <- format_p(
  get_kw("Shannon")
)


label_simpson <- format_p(
  get_kw("Simpson")
)


label_observed <- format_p(
  get_kw("Observed")
)



perm_r2 <- perm$R2[1]

perm_p <- perm$`Pr(>F)`[1]


label_beta <- paste0(
  
  "PERMANOVA R² = ",
  
  formatC(
    perm_r2,
    format = "f",
    digits = 3
  ),
  
  ", P = ",
  
  formatC(
    perm_p,
    format = "f",
    digits = 3
  )
  
)



############################################################
## 7. Alpha plotting function
############################################################


make_alpha <- function(
    dat,
    metric,
    ylabel,
    title
){
  
  
  p <- ggplot(
    dat,
    aes(
      x = Group_display,
      y = .data[[metric]],
      color = Group_display,
      fill = Group_display
    )
  ) +
    
    geom_boxplot(
      width = 0.58,
      linewidth = 0.55,
      outlier.shape = NA,
      alpha = 0.25,
      color = "black"
    ) +
    
    geom_jitter(
      width = 0.11,
      size = 1.8,
      alpha = 0.82
    ) +
    
    scale_color_manual(
      values = group_colors
    ) +
    
    scale_fill_manual(
      values = group_colors
    ) +
    
    labs(
      x = NULL,
      y = ylabel,
      title = title
    ) +
    
    fig_theme +
    
    theme(
      legend.position = "none",
      axis.text.x = element_text(
        angle = 40,
        hjust = 1,
        vjust = 1,
        size = 9
      )
    )
  
  
  return(p)
  
}


############################################################
## 8. Beta plot
############################################################


p_beta <- ggplot(
  pcoa,
  aes(
    x = PCoA1,
    y = PCoA2,
    color = Group_display
  )
) +
  
  stat_ellipse(
    aes(
      group = Group_display
    ),
    type = "norm",
    level = 0.95,
    linewidth = 0.55,
    show.legend = FALSE
  ) +
  
  geom_point(
    size = 1.6,
    alpha = 0.82
  ) +
  
  scale_color_manual(
    values = group_colors
  ) +
  
  labs(
    x = "Bray-Curtis PCoA 1",
    y = "Bray-Curtis PCoA 2",
    title = label_beta
  ) +
  
  fig_theme +
  
  theme(
    legend.position = "right"
  )

############################################################
## 9. Assemble Supplementary Figure
############################################################


pA <- make_alpha(
  
  alpha,
  
  "Shannon",
  
  "Shannon Index",
  
  label_shannon
  
)



pB <- make_alpha(
  
  alpha,
  
  "Simpson",
  
  "Simpson diversity",
  
  label_simpson
  
)



pC <- make_alpha(
  
  alpha,
  
  "Observed",
  
  "Observed features",
  
  label_observed
  
)



pD <- p_beta



supp_fig <- cowplot::plot_grid(
  
  pA,
  
  pB,
  
  pC,
  
  pD,
  
  ncol = 2,
  
  labels = c(
    "A",
    "B",
    "C",
    "D"
  ),
  
  label_size = 14,
  
  label_fontface = "bold",
  
  label_fontfamily = font_family
  
)



############################################################
## 10. Export
############################################################


export_plot(
  
  supp_fig,
  
  "Supplementary_Figure_UC_inflammation_alpha_beta_7KB",
  
  width = 7.2,
  
  height = 6.4
  
)



export_plot(
  
  pA,
  
  "Supp_UC_inflammation_A_Shannon",
  
  width = 3.2,
  
  height = 3.0
  
)


export_plot(
  
  pB,
  
  "Supp_UC_inflammation_B_Simpson",
  
  width = 3.2,
  
  height = 3.0
  
)


export_plot(
  
  pC,
  
  "Supp_UC_inflammation_C_Observed",
  
  width = 3.2,
  
  height = 3.0
  
)


export_plot(
  
  pD,
  
  "Supp_UC_inflammation_D_Bray_PCoA",
  
  width = 3.8,
  
  height = 3.0
  
)



message(
  "Finished Supplementary UC inflammation alpha/beta plotting"
)

message(
  OUT_DIR
)