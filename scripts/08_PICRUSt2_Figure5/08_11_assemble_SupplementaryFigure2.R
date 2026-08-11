#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_11_assemble_SupplementaryFigure2.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Assemble Supplementary Figure 2: MetaCyc pathway LEfSe markers for
## polyp versus UC remission.
############################################################



suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

base_dir <- OUTPUT_ROOT

input_file <- file.path(
  base_dir,
  "figures_progression127_7KB_retention_flagged",
  "MetaCyc_pathway_LEfSe_progression5_selected",
  "supplementary_MetaCyc_pathways_Polyp_UCremission_progression127_7KB.csv"
)

out_dir <- file.path(
  base_dir,
  "figures",
  "main_figures",
  "Figure5_MetaCyc_pathway_LEfSe"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Helvetica"

module_colors <- c(
  "Amino acid metabolism" = "#6F7FA6",
  "Carbohydrate / glycan metabolism" = "#A99E79",
  "Central carbon / fermentation" = "#7DA98B",
  "Cofactor / quinone metabolism" = "#8B79A6",
  "Ethanolamine / choline metabolism" = "#D69A4A",
  "Cell envelope / lipid metabolism" = "#C77DA8",
  "Lipid / membrane / cell envelope" = "#C77DA8",
  "Nitrogen metabolism" = "#7DA9B7",
  "Nucleotide / translation" = "#8DBA91",
  "Sulfur metabolism" = "#B8916B",
  "Other" = "#B8B8B8"
)

format_class <- function(x) {
  recode(
    as.character(x),
    "Polyp" = "Polyp",
    "UC_remission" = "UC remission",
    "UC remission" = "UC remission",
    .default = as.character(x)
  )
}

fix_pathway_module <- function(description, module) {
  desc <- str_to_lower(as.character(description))
  module <- as.character(module)
  
  case_when(
    str_detect(desc, "queuosine biosynthesis") ~ "Nucleotide / translation",
    str_detect(desc, "nad de novo biosynthesis") |
      str_detect(desc, "nad salvage") |
      str_detect(desc, "nad\\+") ~ "Cofactor / quinone metabolism",
    str_detect(desc, "galactarate") |
      str_detect(desc, "glucarate") ~ "Carbohydrate / glycan metabolism",
    str_detect(desc, "lipid iva") |
      str_detect(desc, "peptidoglycan") |
      str_detect(desc, "cdp-diacylglycerol") ~ "Cell envelope / lipid metabolism",
    TRUE ~ module
  )
}

dat0 <- read_csv(input_file, show_col_types = FALSE)

plot_dat <- dat0 %>%
  mutate(
    Enriched_group = format_class(Class),
    Enriched_group = factor(Enriched_group, levels = c("Polyp", "UC remission")),
    LDA_score = if ("abs_LDA" %in% names(.)) as.numeric(abs_LDA) else abs(as.numeric(LDA)),
    Pathway_module = ifelse(is.na(Pathway_module) | Pathway_module == "", "Other", Pathway_module),
    Description_for_module = ifelse(
      !is.na(Description_display) & Description_display != "",
      Description_display,
      Description
    ),
    Pathway_module = fix_pathway_module(Description_for_module, Pathway_module),
    Pathway_label = str_squish(str_replace_all(Description_for_module, "\\s+", " ")),
    Pathway_label = str_wrap(Pathway_label, width = 46)
  ) %>%
  filter(!is.na(Enriched_group), !is.na(LDA_score), !is.na(Pathway_label))

if ("Main_plot_eligible" %in% names(plot_dat)) {
  plot_dat <- plot_dat %>% filter(Main_plot_eligible)
}
if ("Manual_artifact_removed" %in% names(plot_dat)) {
  plot_dat <- plot_dat %>% filter(!Manual_artifact_removed)
}
if ("artifact_flag" %in% names(plot_dat)) {
  plot_dat <- plot_dat %>% filter(!artifact_flag)
}

plot_dat <- plot_dat %>%
  arrange(Enriched_group, LDA_score) %>%
  mutate(
    row_id = row_number(),
    Pathway_label_unique = paste0(Pathway_label, "___", row_id),
    Pathway_label_unique = factor(Pathway_label_unique, levels = unique(Pathway_label_unique))
  )

modules_use <- unique(as.character(plot_dat$Pathway_module))
colors_use <- module_colors[modules_use]
missing_modules <- modules_use[is.na(colors_use)]
if (length(missing_modules) > 0) {
  extra_cols <- rep("#B8B8B8", length(missing_modules))
  names(extra_cols) <- missing_modules
  colors_use <- c(colors_use[!is.na(colors_use)], extra_cols)
}

p <- ggplot(plot_dat, aes(x = LDA_score, y = Pathway_label_unique, fill = Pathway_module)) +
  geom_col(width = 0.70) +
  facet_grid(
    Enriched_group ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_y_discrete(labels = function(x) str_replace(x, "___[0-9]+$", "")) +
  scale_fill_manual(values = colors_use, drop = FALSE) +
  scale_x_continuous(breaks = seq(0, 4, 1), expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "MetaCyc pathway markers in Polyp and UC remission",
    x = "LDA score",
    y = NULL,
    fill = "Pathway module"
  ) +
  theme_classic(base_size = 8, base_family = font_family) +
  theme(
    text = element_text(color = "black"),
    axis.title = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 7.2, color = "black", lineheight = 0.9),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.title = element_text(size = 10, face = "bold", hjust = 0, margin = margin(b = 5)),
    strip.placement = "outside",
    strip.background.y = element_blank(),
    strip.text.y.left = element_text(angle = 0, size = 8.2, face = "bold", margin = margin(r = 5)),
    legend.position = "right",
    legend.title = element_text(size = 8.5, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(4.0, "mm"),
    panel.spacing.y = unit(3.2, "mm"),
    plot.margin = margin(5, 5, 5, 5)
  )

export_plot <- function(plot, filename, width, height, dpi = 600) {
  pdf_file <- file.path(out_dir, paste0(filename, ".pdf"))
  png_file <- file.path(out_dir, paste0(filename, ".png"))
  tif_file <- file.path(out_dir, paste0(filename, ".tiff"))
  
  grDevices::pdf(pdf_file, width = width, height = height, family = font_family, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
  
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(png_file, width = width, height = height, units = "in", res = dpi, background = "white")
    print(plot)
    grDevices::dev.off()
    
    ragg::agg_tiff(tif_file, width = width, height = height, units = "in", res = dpi, compression = "lzw", background = "white")
    print(plot)
    grDevices::dev.off()
  } else {
    bitmap_type <- if (Sys.info()[["sysname"]] == "Darwin") "quartz" else "cairo"
    
    grDevices::png(png_file, width = width, height = height, units = "in", res = dpi, bg = "white", type = bitmap_type)
    print(plot)
    grDevices::dev.off()
    
    grDevices::tiff(tif_file, width = width, height = height, units = "in", res = dpi, bg = "white", compression = "lzw", type = bitmap_type)
    print(plot)
    grDevices::dev.off()
  }
}

export_plot(
  p,
  "Supplementary_Figure_MetaCyc_pathway_LEfSe_Polyp_UCremission",
  width = 9.8,
  height = 5.2
)

write_csv(
  plot_dat,
  file.path(out_dir, "Supplementary_Figure_MetaCyc_pathway_LEfSe_Polyp_UCremission_plot_data.csv")
)

message("Done. Outputs written to: ", out_dir)