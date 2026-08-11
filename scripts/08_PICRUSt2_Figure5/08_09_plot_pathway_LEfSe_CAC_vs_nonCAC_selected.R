
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 08_09_plot_pathway_LEfSe_CAC_vs_nonCAC_selected.R
##
## Module 08 - PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)
##
## Select the top pathways per group and functional module, remove artefactual
## pathways, and write the panel source table.
############################################################


suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(forcats)
})

analysis_root <- OUTPUT_ROOT

input_dir <- file.path(
  analysis_root,
  "LEfSe_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB"
)

input_file <- file.path(
  input_dir,
  "FINAL_LEfSe_clean_markers_MetaCyc_pathways_CA23_vs_nonCA23_oldstyle_7KB.csv"
)

out_dir <- file.path(
  analysis_root,
  "figures_progression127_7KB_retention_flagged",
  "MetaCyc_pathway_LEfSe_CA23_vs_nonCA23_selected"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_file))

markers <- readr::read_csv(input_file, show_col_types = FALSE)

required_cols <- c(
  "Class", "Pathway_ID", "Description_display",
  "Pathway_module", "LDA", "p_value"
)

missing_cols <- setdiff(required_cols, colnames(markers))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

if (!"q_value_BH" %in% colnames(markers)) {
  markers <- markers %>%
    mutate(q_value_BH = p.adjust(p_value, method = "BH"))
}

markers_clean <- markers %>%
  mutate(
    Class = factor(Class, levels = c("nonCA", "CA")),
    Pathway_module_original = Pathway_module,
    Pathway_module = case_when(
      Pathway_module %in% c(
        "Lipid / membrane / cell envelope",
        "Cell envelope / lipid metabolism"
      ) ~ "Cell envelope / lipid metabolism",
      TRUE ~ Pathway_module
    ),
    LDA = abs(as.numeric(LDA)),
    p_value = as.numeric(p_value),
    q_value_BH = as.numeric(q_value_BH),
    Description_display = str_replace_all(Description_display, "&beta;", "beta")
  ) %>%
  filter(!is.na(Class), !is.na(Pathway_module), !is.na(LDA))

module_order <- c(
  "Amino acid metabolism",
  "Carbohydrate / glycan metabolism",
  "Central carbon / fermentation",
  "Cofactor / quinone metabolism",
  "Ethanolamine / choline metabolism",
  "Cell envelope / lipid metabolism",
  "Nitrogen metabolism",
  "Nucleotide / translation",
  "Sulfur metabolism"
)

artifact_pattern <- regex(
  paste(
    c(
      "photosynthetic",
      "photosynthesis",
      "phototrophic",
      "plastidic",
      "\\bplants?\\b",
      "plant peroxisome",
      "C4 photosynthetic",
      "sulfide oxidation.*phototrophic",
      "sucrose biosynthesis"
    ),
    collapse = "|"
  ),
  ignore_case = TRUE
)

selection_audit <- tibble(
  Step = "Initial LEfSe markers",
  n_pathways = nrow(markers_clean)
)

markers_no_other <- markers_clean %>%
  filter(Pathway_module != "Other")

selection_audit <- bind_rows(
  selection_audit,
  tibble(
    Step = "After excluding Other module",
    n_pathways = nrow(markers_no_other)
  )
)

selected_top2 <- markers_no_other %>%
  group_by(Class, Pathway_module) %>%
  arrange(desc(LDA), p_value, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

selection_audit <- bind_rows(
  selection_audit,
  tibble(
    Step = "After top2 by Class x Pathway_module",
    n_pathways = nrow(selected_top2)
  )
)

selected_final <- selected_top2 %>%
  mutate(
    artifact_flag = str_detect(Description_display, artifact_pattern)
  ) %>%
  filter(!artifact_flag) %>%
  mutate(
    Pathway_module = factor(Pathway_module, levels = module_order),
    Class = factor(Class, levels = c("nonCA", "CA"))
  ) %>%
  arrange(Class, Pathway_module, desc(LDA))

removed_artifacts <- selected_top2 %>%
  mutate(
    artifact_flag = str_detect(Description_display, artifact_pattern)
  ) %>%
  filter(artifact_flag) %>%
  arrange(Class, Pathway_module, desc(LDA))

selection_audit <- bind_rows(
  selection_audit,
  tibble(
    Step = "After removing predefined artifacts",
    n_pathways = nrow(selected_final)
  ),
  tibble(
    Step = "CA selected",
    n_pathways = sum(selected_final$Class == "CA")
  ),
  tibble(
    Step = "nonCA selected",
    n_pathways = sum(selected_final$Class == "nonCA")
  )
)

selected_ca <- selected_final %>%
  filter(Class == "CA")

selected_nonca <- selected_final %>%
  filter(Class == "nonCA")

write_csv(
  selected_final,
  file.path(out_dir, "selected_MetaCyc_pathways_CA23_vs_nonCA23_top2_by_group_module_7KB.csv")
)

write_csv(
  selected_ca,
  file.path(out_dir, "selected_CA_enriched_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv")
)

write_csv(
  selected_nonca,
  file.path(out_dir, "selected_nonCA_enriched_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv")
)

write_csv(
  selection_audit,
  file.path(out_dir, "selection_audit_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv")
)

write_csv(
  removed_artifacts,
  file.path(out_dir, "removed_artifacts_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv")
)

plot_df <- selected_final %>%
  mutate(
    signed_LDA = if_else(Class == "nonCA", -LDA, LDA),
    Plot_label = str_wrap(Description_display, width = 42)
  ) %>%
  arrange(Class, Pathway_module, LDA) %>%
  mutate(
    Plot_label_unique = make.unique(Plot_label),
    Plot_label_unique = factor(Plot_label_unique, levels = Plot_label_unique)
  )

module_colors <- c(
  "Amino acid metabolism" = "#B45F06",
  "Carbohydrate / glycan metabolism" = "#8E7CC3",
  "Central carbon / fermentation" = "#CC0000",
  "Cofactor / quinone metabolism" = "#6AA84F",
  "Ethanolamine / choline metabolism" = "#E69138",
  "Cell envelope / lipid metabolism" = "#3C78D8",
  "Nitrogen metabolism" = "#674EA7",
  "Nucleotide / translation" = "#45818E",
  "Sulfur metabolism" = "#7F6000"
)

p <- ggplot(plot_df, aes(x = signed_LDA, y = Plot_label_unique, fill = Pathway_module)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.15) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey30") +
  scale_x_continuous(
    name = "LDA score",
    labels = abs,
    expand = expansion(mult = c(0.04, 0.08))
  ) +
  scale_fill_manual(values = module_colors, drop = FALSE, name = "Functional module") +
  labs(
    y = NULL,
    title = "Predicted MetaCyc pathway differences between CA and nonCA mucosa",
    subtitle = "Selected by top 2 LDA pathways within each group × functional module; Other module excluded"
  ) +
  annotate(
    "text",
    x = min(plot_df$signed_LDA) * 0.85,
    y = nrow(plot_df) + 0.6,
    label = "nonCA-enriched",
    hjust = 0.5,
    size = 3.8,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = max(plot_df$signed_LDA) * 0.85,
    y = nrow(plot_df) + 0.6,
    label = "CA-enriched",
    hjust = 0.5,
    size = 3.8,
    fontface = "bold"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 12, 8, 8)
  )

pdf_file <- file.path(out_dir, "Fig_MetaCyc_pathway_LEfSe_CA23_vs_nonCA23_selected.pdf")
png_file <- file.path(out_dir, "Fig_MetaCyc_pathway_LEfSe_CA23_vs_nonCA23_selected.png")
tiff_file <- file.path(out_dir, "Fig_MetaCyc_pathway_LEfSe_CA23_vs_nonCA23_selected.tiff")

ggsave(pdf_file, p, width = 9.5, height = 8.5, device = "pdf")
ggsave(png_file, p, width = 9.5, height = 8.5, dpi = 600)
ggsave(tiff_file, p, width = 9.5, height = 8.5, dpi = 600, compression = "lzw")

cat("\n============================================================\n")
cat("CA vs nonCA pathway LEfSe representative selection finished\n")
cat("============================================================\n\n")

cat("Selection audit:\n")
print(selection_audit)

cat("\nSelected pathways by class:\n")
print(selected_final %>% count(Class))

cat("\nSelected pathways by class and module:\n")
print(selected_final %>% count(Class, Pathway_module))

if (nrow(removed_artifacts) > 0) {
  cat("\nRemoved predefined artifacts:\n")
  print(removed_artifacts %>% select(Class, Pathway_ID, Description_display, Pathway_module, LDA))
}

cat("\nSelected pathway table:\n")
print(selected_final %>% select(Class, Pathway_ID, Description_display, Pathway_module, LDA, p_value, q_value_BH))

cat("\nFiles written:\n")
cat("Selected table:  ", file.path(out_dir, "selected_MetaCyc_pathways_CA23_vs_nonCA23_top2_by_group_module_7KB.csv"), "\n")
cat("CA table:        ", file.path(out_dir, "selected_CA_enriched_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv"), "\n")
cat("nonCA table:     ", file.path(out_dir, "selected_nonCA_enriched_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv"), "\n")
cat("Audit:           ", file.path(out_dir, "selection_audit_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv"), "\n")
cat("Removed artifact:", file.path(out_dir, "removed_artifacts_MetaCyc_pathways_CA23_vs_nonCA23_7KB.csv"), "\n")
cat("PDF:             ", pdf_file, "\n")
cat("PNG:             ", png_file, "\n")
cat("TIFF:            ", tiff_file, "\n")
cat("\nDone.\n")