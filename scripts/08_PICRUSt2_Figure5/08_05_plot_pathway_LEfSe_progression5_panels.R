
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_05_plot_pathway_LEfSe_progression5_panels.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Write the panel source tables: main panel (UC active / dysplasia / CAC)
## and supplementary panel (polyp / UC remission).
############################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
  library(ggplot2)
})

PROJECT_ROOT <- PROJECT_ROOT
MAIN_OUT <- OUTPUT_ROOT  ## 7KB rerun output root

LEFSE_DIR <- file.path(
  MAIN_OUT,
  "LEfSe_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB"
)

IN_CLEAN <- file.path(
  LEFSE_DIR,
  "FINAL_LEfSe_clean_markers_MetaCyc_pathways_progression5_oldstyle_currentdata_progression127_7KB.csv"
)

OUT_DIR <- file.path(
  MAIN_OUT,
  "figures_progression127_7KB_retention_flagged/MetaCyc_pathway_LEfSe_progression5_selected"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_ALL_SELECTED <- file.path(
  OUT_DIR,
  "selected_MetaCyc_pathways_progression5_all_groups_top2_by_group_module_progression127_7KB.csv"
)
OUT_MAIN_TABLE <- file.path(
  OUT_DIR,
  "main_figure_MetaCyc_pathways_UCactive_Dysplasia_CA_progression127_7KB.csv"
)
OUT_SUPP_TABLE <- file.path(
  OUT_DIR,
  "supplementary_MetaCyc_pathways_Polyp_UCremission_progression127_7KB.csv"
)
OUT_AUDIT <- file.path(
  OUT_DIR,
  "selection_audit_MetaCyc_pathways_progression5_progression127_7KB.csv"
)

OUT_MAIN_PDF <- file.path(
  OUT_DIR,
  "Fig_MetaCyc_pathway_LEfSe_progression5_main_UCactive_Dysplasia_CA.pdf"
)
OUT_MAIN_PNG <- file.path(
  OUT_DIR,
  "Fig_MetaCyc_pathway_LEfSe_progression5_main_UCactive_Dysplasia_CA.png"
)
OUT_MAIN_TIFF <- file.path(
  OUT_DIR,
  "Fig_MetaCyc_pathway_LEfSe_progression5_main_UCactive_Dysplasia_CA.tiff"
)

OUT_SUPP_PDF <- file.path(
  OUT_DIR,
  "FigS_MetaCyc_pathway_LEfSe_progression5_supp_Polyp_UCremission.pdf"
)
OUT_SUPP_PNG <- file.path(
  OUT_DIR,
  "FigS_MetaCyc_pathway_LEfSe_progression5_supp_Polyp_UCremission.png"
)
OUT_SUPP_TIFF <- file.path(
  OUT_DIR,
  "FigS_MetaCyc_pathway_LEfSe_progression5_supp_Polyp_UCremission.tiff"
)

stopifnot(file.exists(IN_CLEAN))

GROUP_LEVELS <- c("Polyp", "UC_remission", "UC_active", "Dysplasia", "CA")
MAIN_GROUPS <- c("UC_active", "Dysplasia", "CA")
SUPP_GROUPS <- c("Polyp", "UC_remission")

MODULE_LEVELS <- c(
  "Amino acid metabolism",
  "Carbohydrate / glycan metabolism",
  "Central carbon / fermentation",
  "Cofactor / quinone metabolism",
  "Ethanolamine / choline metabolism",
  "Lipid / membrane / cell envelope",
  "Cell envelope / lipid metabolism",
  "Nitrogen metabolism",
  "Nucleotide / translation",
  "Sulfur metabolism"
)

module_colors <- c(
  "Amino acid metabolism" = "#7B3294",
  "Carbohydrate / glycan metabolism" = "#A6611A",
  "Central carbon / fermentation" = "#018571",
  "Cofactor / quinone metabolism" = "#5E3C99",
  "Ethanolamine / choline metabolism" = "#D95F02",
  "Lipid / membrane / cell envelope" = "#C51B7D",
  "Cell envelope / lipid metabolism" = "#C51B7D",
  "Nitrogen metabolism" = "#1F78B4",
  "Nucleotide / translation" = "#4D9221",
  "Sulfur metabolism" = "#E6AB02"
)

class_colors <- c(
  "Polyp" = "#8DA0CB",
  "UC_remission" = "#66C2A5",
  "UC_active" = "#FC8D62",
  "Dysplasia" = "#E78AC3",
  "CA" = "#E41A1C"
)

make_plot_label <- function(pathway_id, description_display) {
  paste0(description_display, " [", pathway_id, "]")
}

is_manual_artifact <- function(pathway_id, description_display) {
  x <- str_to_lower(paste(pathway_id, description_display))
  str_detect(x, "c4 photosynthetic carbon assimilation cycle") |
    str_detect(x, "superpathway of sulfide oxidation") |
    str_detect(x, "phototrophic sulfur bacteria")
}

wrap_label <- function(x, width = 48) {
  stringr::str_wrap(x, width = width)
}

plot_selected_pathways <- function(df, title_text, subtitle_text, base_size = 9) {
  plot_df <- df %>%
    mutate(
      Class = factor(Class, levels = GROUP_LEVELS),
      Pathway_module = factor(Pathway_module, levels = MODULE_LEVELS),
      Plot_label = wrap_label(make_plot_label(Pathway_ID, Description_display), 52),
      Plot_label_ordered = factor(
        Plot_label,
        levels = rev(unique(Plot_label[order(as.numeric(Class), LDA)]))
      )
    )
  
  ggplot(
    plot_df,
    aes(x = LDA, y = Plot_label_ordered, fill = Pathway_module)
  ) +
    geom_col(width = 0.72, color = "grey25", linewidth = 0.15) +
    facet_grid(
      Class ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_fill_manual(values = module_colors, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = "LDA score",
      y = NULL,
      fill = "Pathway module"
    ) +
    theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 3),
      plot.subtitle = element_text(size = base_size + 1),
      strip.placement = "outside",
      strip.background = element_rect(fill = "grey92", color = "grey50"),
      strip.text.y.left = element_text(angle = 0, face = "bold"),
      axis.text.y = element_text(size = base_size - 1),
      axis.text.x = element_text(size = base_size),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

clean_markers <- read_csv(IN_CLEAN, show_col_types = FALSE) %>%
  mutate(
    Class = factor(as.character(Class), levels = GROUP_LEVELS),
    Pathway_module = as.character(Pathway_module),
    Description_display = as.character(Description_display),
    Pathway_ID = as.character(Pathway_ID),
    LDA = as.numeric(LDA),
    abs_LDA = abs(as.numeric(LDA)),
    p_value = as.numeric(p_value),
    q_value_BH = if ("q_value_BH" %in% names(.)) as.numeric(q_value_BH) else p.adjust(p_value, method = "BH")
  )

selection_pool <- clean_markers %>%
  filter(
    !is.na(Class),
    !is.na(Pathway_module),
    Pathway_module != "Other"
  )

top2_by_group_module <- selection_pool %>%
  group_by(Class, Pathway_module) %>%
  arrange(desc(LDA), p_value, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

selected_all <- top2_by_group_module %>%
  mutate(
    Manual_artifact_removed = is_manual_artifact(Pathway_ID, Description_display)
  ) %>%
  filter(!Manual_artifact_removed) %>%
  mutate(
    Selection_rule = "LEfSe marker; Other module excluded; top 2 by LDA within Class x Pathway_module; predefined photosynthesis/phototrophic sulfur artifacts removed without replacement",
    Class = factor(Class, levels = GROUP_LEVELS)
  ) %>%
  arrange(Class, Pathway_module, desc(LDA), p_value)

main_selected <- selected_all %>%
  filter(as.character(Class) %in% MAIN_GROUPS) %>%
  mutate(Class = factor(as.character(Class), levels = MAIN_GROUPS)) %>%
  arrange(Class, Pathway_module, desc(LDA), p_value)

supp_selected <- selected_all %>%
  filter(as.character(Class) %in% SUPP_GROUPS) %>%
  mutate(Class = factor(as.character(Class), levels = SUPP_GROUPS)) %>%
  arrange(Class, Pathway_module, desc(LDA), p_value)

audit <- bind_rows(
  tibble(
    Step = "Initial LEfSe markers",
    n_pathways = nrow(clean_markers)
  ),
  tibble(
    Step = "After excluding Other module",
    n_pathways = nrow(selection_pool)
  ),
  tibble(
    Step = "After top2 by Class x Pathway_module",
    n_pathways = nrow(top2_by_group_module)
  ),
  tibble(
    Step = "After removing predefined artifacts",
    n_pathways = nrow(selected_all)
  ),
  tibble(
    Step = "Main figure UC_active/Dysplasia/CA",
    n_pathways = nrow(main_selected)
  ),
  tibble(
    Step = "Supplementary Polyp/UC_remission",
    n_pathways = nrow(supp_selected)
  )
)

write_csv(selected_all, OUT_ALL_SELECTED)
write_csv(main_selected, OUT_MAIN_TABLE)
write_csv(supp_selected, OUT_SUPP_TABLE)
write_csv(audit, OUT_AUDIT)

main_plot <- plot_selected_pathways(
  main_selected,
  title_text = "Selected predicted MetaCyc pathways across UC activity, dysplasia, and carcinoma",
  subtitle_text = "Top pathways by LDA within each enriched group and pathway module; Other module excluded"
)

supp_plot <- plot_selected_pathways(
  supp_selected,
  title_text = "Supplementary selected predicted MetaCyc pathways in Polyp and UC remission",
  subtitle_text = "Top pathways by LDA within each enriched group and pathway module; Other module excluded",
  base_size = 8.5
)

ggsave(OUT_MAIN_PDF, main_plot, width = 9.2, height = 8.6, device = "pdf")
ggsave(OUT_MAIN_PNG, main_plot, width = 9.2, height = 8.6, dpi = 600)
ggsave(OUT_MAIN_TIFF, main_plot, width = 9.2, height = 8.6, dpi = 600, compression = "lzw")

ggsave(OUT_SUPP_PDF, supp_plot, width = 9.2, height = 7.2, device = "pdf")
ggsave(OUT_SUPP_PNG, supp_plot, width = 9.2, height = 7.2, dpi = 600)
ggsave(OUT_SUPP_TIFF, supp_plot, width = 9.2, height = 7.2, dpi = 600, compression = "lzw")

cat("\n============================================================\n")
cat("Progression5 pathway LEfSe main/supp figure selection finished\n")
cat("============================================================\n\n")

cat("Selection audit:\n")
print(audit, n = Inf)

cat("\nMain figure selected pathways by class:\n")
print(main_selected %>% count(Class), n = Inf)

cat("\nSupplementary selected pathways by class:\n")
print(supp_selected %>% count(Class), n = Inf)

cat("\nMain figure selected pathway table:\n")
print(
  main_selected %>%
    select(Class, Pathway_ID, Description_display, Pathway_module, LDA, p_value, q_value_BH),
  n = Inf,
  width = Inf
)

cat("\nFiles written:\n")
cat("All selected table: ", OUT_ALL_SELECTED, "\n")
cat("Main table:         ", OUT_MAIN_TABLE, "\n")
cat("Supplement table:   ", OUT_SUPP_TABLE, "\n")
cat("Audit:              ", OUT_AUDIT, "\n")
cat("Main PDF:           ", OUT_MAIN_PDF, "\n")
cat("Main PNG:           ", OUT_MAIN_PNG, "\n")
cat("Main TIFF:          ", OUT_MAIN_TIFF, "\n")
cat("Supp PDF:           ", OUT_SUPP_PDF, "\n")
cat("Supp PNG:           ", OUT_SUPP_PNG, "\n")
cat("Supp TIFF:          ", OUT_SUPP_TIFF, "\n\n")
cat("Done.\n")