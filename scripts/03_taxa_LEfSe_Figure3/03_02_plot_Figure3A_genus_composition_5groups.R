#!/usr/bin/env Rscript

## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))



############################################################
## 03_02_plot_Figure3A_genus_composition_5groups.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
##
## Final manuscript-style 5-group genus composition panel
## for the 7KB progression127 rerun.
############################################################

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
})

############################################################
## 1. Paths
############################################################

base_dir <- file.path(PROJECT_ROOT, "output/analysis")

color_mapping_file <- file.path(
  rprojroot::find_root(rprojroot::has_file("config.R")),
  "scripts", "03_taxa_LEfSe_Figure3", "genus_color_mapping.R"
)

if (!file.exists(color_mapping_file)) {
  stop(
    "Shared genus color mapping file not found: ",
    color_mapping_file
  )
}

source(color_mapping_file)

if (!exists("GENUS_COLORS")) {
  stop(
    "GENUS_COLORS was not created after sourcing: ",
    color_mapping_file
  )
}

input_file <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "genus_composition_5groups_top20_remove_others_renormalized",
  "genus_group_mean_top20_Others_removed_renormalized100_7KB_progression127.csv"
)

out_dir <- file.path(
  base_dir,
  "02_Figure3_taxa_LEfSe",
  "Figure3_final_panels",
  "Figure3A_genus_composition_5groups"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(input_file)) {
  stop(
    "Input file not found: ",
    input_file
  )
}

message(
  "Using input file: ",
  input_file
)

############################################################
## 2. Read data
############################################################

dat0 <- read_csv(
  input_file,
  show_col_types = FALSE
)

message(
  "Input columns: ",
  paste(
    names(dat0),
    collapse = ", "
  )
)

############################################################
## 3. Standardize columns
##
## Kept the same as the original final plotting script.
############################################################

needed_cols <- c(
  "Progression5",
  "Genus_plot",
  "renormalized_percent"
)

missing_cols <- setdiff(
  needed_cols,
  names(dat0)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing expected column(s): ",
    paste(
      missing_cols,
      collapse = ", "
    ),
    "\nAvailable columns are: ",
    paste(
      names(dat0),
      collapse = ", "
    )
  )
}

plot_dat <- dat0 %>%
  transmute(
    Group_raw = as.character(
      Progression5
    ),
    Genus = as.character(
      Genus_plot
    ),
    Display_abundance = as.numeric(
      renormalized_percent
    )
  ) %>%
  mutate(
    Group = recode(
      Group_raw,
      "Polyp" = "Polyp",
      "UC_remission" = "UC remission",
      "UC remission" = "UC remission",
      "UC_active" = "UC active",
      "UC active" = "UC active",
      "Dysplasia" = "Dysplasia",
      "CA" = "CAC",
      "CAC" = "CAC",
      .default = Group_raw
    ),
    
    Group = factor(
      Group,
      levels = c(
        "Polyp",
        "UC remission",
        "UC active",
        "Dysplasia",
        "CAC"
      )
    ),
    
    ########################################################
    ## Harmonize 7KB SILVA genus names to the naming
    ## convention used in the original manuscript palette.
    ##
    ## IMPORTANT:
    ## These are naming aliases only.
    ## No taxa are removed, merged, or changed biologically.
    ########################################################
    
    Genus = recode(
      Genus,
      "UCG-005" = "UCG_005",
      "Escherichia-Shigella" = "Escherichia_Shigella",
      "Christensenellaceae_R-7_group" =
        "Christensenellaceae_R_7_group",
      .default = Genus
    ),
    
    ########################################################
    ## Original final Figure 2A display convention:
    ## underscores -> spaces
    ########################################################
    
    Genus = str_replace_all(
      Genus,
      "_",
      " "
    ),
    
    Genus = str_squish(
      Genus
    )
  ) %>%
  filter(
    !is.na(Group),
    !is.na(Genus),
    !is.na(Display_abundance)
  ) %>%
  filter(
    !str_detect(
      str_to_lower(Genus),
      "^others?$|^other$"
    )
  )

############################################################
## 4. Abundance scale check
##
## Kept exactly from the original final plotting script.
############################################################

if (
  max(
    plot_dat$Display_abundance,
    na.rm = TRUE
  ) <= 1.5
) {
  plot_dat <- plot_dat %>%
    mutate(
      Display_abundance =
        Display_abundance * 100
    )
}

############################################################
## 5. Re-normalize within each group to exactly 100%
##
## This is also present in the original final plotting script.
############################################################

plot_dat <- plot_dat %>%
  group_by(
    Group
  ) %>%
  mutate(
    Display_abundance =
      Display_abundance /
      sum(
        Display_abundance,
        na.rm = TRUE
      ) *
      100
  ) %>%
  ungroup()

############################################################
## 6. Order genera
##
## Original rule:
## overall mean abundance, ascending factor order.
## Small genera are at bottom; dominant genera are at top.
############################################################

genus_order <- plot_dat %>%
  group_by(
    Genus
  ) %>%
  summarise(
    mean_abundance = mean(
      Display_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    mean_abundance
  ) %>%
  pull(
    Genus
  )

plot_dat <- plot_dat %>%
  mutate(
    Genus = factor(
      Genus,
      levels = genus_order
    )
  )

############################################################
## 7. Shared genus palette
##
## Reuse original manuscript palette unchanged.
############################################################

genus_palette_all <- GENUS_COLORS

names(genus_palette_all) <- str_replace_all(
  names(genus_palette_all),
  "_",
  " "
)

names(genus_palette_all) <- str_squish(
  names(genus_palette_all)
)

displayed_genera <- levels(
  plot_dat$Genus
)

missing_colors <- setdiff(
  displayed_genera,
  names(genus_palette_all)
)

if (length(missing_colors) > 0) {
  stop(
    "Missing shared colors for the following genera:\n",
    paste(
      missing_colors,
      collapse = ", "
    )
  )
}

genus_palette <- genus_palette_all[
  displayed_genera
]

############################################################
## 8. Audit before plotting
############################################################

message(
  "Number of displayed genera: ",
  length(displayed_genera)
)

message(
  "Displayed genera:"
)

print(
  tibble(
    Order = seq_along(displayed_genera),
    Genus = displayed_genera
  ),
  n = Inf
)

group_sum_check <- plot_dat %>%
  group_by(
    Group
  ) %>%
  summarise(
    Total_percent = sum(
      Display_abundance
    ),
    .groups = "drop"
  )

message(
  "Group abundance sums after final re-normalization:"
)

print(
  group_sum_check,
  n = Inf
)

############################################################
## 9. Plot
##
## Kept the same as the original final manuscript script.
############################################################

font_family <- "Helvetica"

p <- ggplot(
  plot_dat,
  aes(
    x = Group,
    y = Display_abundance,
    fill = Genus
  )
) +
  geom_col(
    width = 0.68,
    color = "white",
    linewidth = 0.22
  ) +
  scale_fill_manual(
    values = genus_palette,
    drop = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(
      0,
      100,
      25
    ),
    expand = expansion(
      mult = c(
        0,
        0.01
      )
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      100
    ),
    clip = "off"
  ) +
  labs(
    x = NULL,
    y = "Relative abundance of top 20 genera (%)",
    fill = "Genus"
  ) +
  theme_classic(
    base_size = 8,
    base_family = font_family
  ) +
  theme(
    text = element_text(
      color = "black"
    ),
    
    axis.title.y = element_text(
      size = 9,
      face = "bold",
      margin = margin(
        r = 6
      )
    ),
    
    axis.text.x = element_text(
      size = 9,
      color = "black",
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    
    axis.line = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.45,
      color = "black"
    ),
    
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      size = 9,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 7.4
    ),
    
    legend.key.size = unit(
      3.3,
      "mm"
    ),
    
    legend.spacing.y = unit(
      0.8,
      "mm"
    ),
    
    plot.margin = margin(
      6,
      8,
      6,
      6
    )
  ) +
  guides(
    fill = guide_legend(
      ncol = 1,
      reverse = FALSE,
      keyheight = unit(
        3.5,
        "mm"
      ),
      keywidth = unit(
        4.5,
        "mm"
      )
    )
  )

print(p)

############################################################
## 10. Export
##
## Same final dimensions as original:
## 6.8 x 4.4 inches
##
## PDF uses base grDevices::pdf(), not Cairo.
############################################################

export_plot <- function(
    plot,
    filename,
    width,
    height,
    dpi = 600
) {
  
  pdf_file <- file.path(
    out_dir,
    paste0(
      filename,
      ".pdf"
    )
  )
  
  png_file <- file.path(
    out_dir,
    paste0(
      filename,
      ".png"
    )
  )
  
  tif_file <- file.path(
    out_dir,
    paste0(
      filename,
      ".tiff"
    )
  )
  
  ##########################################################
  ## PDF
  ##########################################################
  
  grDevices::pdf(
    pdf_file,
    width = width,
    height = height,
    family = font_family,
    useDingbats = FALSE
  )
  
  print(plot)
  
  grDevices::dev.off()
  
  ##########################################################
  ## PNG / TIFF
  ##
  ## Same preference as original final script:
  ## ragg first.
  ##
  ## On macOS fallback uses quartz, not Cairo/X11.
  ##########################################################
  
  if (
    requireNamespace(
      "ragg",
      quietly = TRUE
    )
  ) {
    
    ragg::agg_png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      background = "white"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    ragg::agg_tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      background = "white"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
  } else {
    
    if (
      Sys.info()[["sysname"]] != "Darwin"
    ) {
      stop(
        paste0(
          "Package 'ragg' is not installed. ",
          "This fallback is configured only for macOS quartz. ",
          "Please install ragg before running on a non-macOS system."
        )
      )
    }
    
    grDevices::png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white",
      type = "quartz"
    )
    
    print(plot)
    
    grDevices::dev.off()
    
    grDevices::tiff(
      tif_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white",
      compression = "lzw",
      type = "quartz"
    )
    
    print(plot)
    
    grDevices::dev.off()
  }
}

############################################################
## 11. Export final panel
############################################################

output_stem <- paste0(
  "Figure3A_genus_composition_5groups_",
  "top20_others_removed_renormalized_7KB"
)

export_plot(
  p,
  output_stem,
  width = 6.8,
  height = 4.4
)

############################################################
## 12. Save plotting data
############################################################

plot_data_file <- file.path(
  out_dir,
  "Figure3A_genus_composition_5groups_plot_data_7KB.csv"
)

write_csv(
  plot_dat,
  plot_data_file
)

############################################################
## 13. Save final ggplot object
##
## This RDS will later be used in the final multi-panel
## assembly script, exactly as in the original workflow.
############################################################

figure3A_rds_file <- file.path(
  out_dir,
  "Figure3A_genus_composition_5groups_shared_colors_7KB.rds"
)

saveRDS(
  p,
  figure3A_rds_file
)

############################################################
## 14. Finish
############################################################

message(
  "Done. Outputs written to: ",
  out_dir
)

message(
  "Saved Figure 3A plot object: ",
  figure3A_rds_file
)