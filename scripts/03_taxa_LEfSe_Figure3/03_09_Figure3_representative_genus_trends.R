
## ---------------------------------------------------------------
## Repository configuration (added during repository preparation)
## Defines PROJECT_ROOT. Set the FFPE_PROJECT_ROOT environment
## variable, or edit config.R, to point at your local data copy.
## ---------------------------------------------------------------
source(file.path(rprojroot::find_root(rprojroot::has_file("config.R")), "config.R"))


############################################################
## 03_09_Figure3_representative_genus_trends.R
##
## Module 03 - Taxonomic composition and LEfSe (Figure 3)
##
## Purpose:
## Re-run the original manuscript representative-genus panels using the
## fixed 7KB progression127 dataset.
##
## IMPORTANT:
## This script follows the ACTUAL plotting code of the original final script.
##
## Display:
## - individual samples
## - median
## - IQR
## - overall Kruskal-Wallis P
############################################################


# ==============================================================================
# 0. Required packages
# ==============================================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tidyr",
  "tibble",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "The following R packages are missing:\n",
      paste(
        missing_packages,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# ==============================================================================
# 1. Load the original unified manuscript figure style
# ==============================================================================

style_script <- file.path(
  PROJECT_ROOT,
  "script",
  "analysis",
  "00_figure_style_master.R"
)

if (!file.exists(style_script)) {
  stop(
    paste0(
      "The shared manuscript figure-style script was not found:\n",
      style_script
    ),
    call. = FALSE
  )
}

source(style_script)


# ==============================================================================
# 2. Input and output paths
# ==============================================================================

input_file <- file.path(
  PROJECT_ROOT,
  "output",
  "analysis",
  "02_Figure3_taxa_LEfSe",
  "genus_composition_5groups_top20_remove_others_renormalized",
  "sample_level_genus_relative_abundance_7KB_progression127.csv"
)

output_dir <- file.path(
  PROJECT_ROOT,
  "output",
  "analysis",
  "02_Figure3_taxa_LEfSe",
  "Figure3_final_panels",
  "Figure3D_representative_genus_trends"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==============================================================================
# 3. Confirm input file exists
# ==============================================================================

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file was not found:\n",
      input_file
    ),
    call. = FALSE
  )
}

message("Input file found:")
message(input_file)


# ==============================================================================
# 4. Read sample-level genus abundance table
# ==============================================================================

sample_level_raw <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)


# ==============================================================================
# 5. Validate required columns
# ==============================================================================

required_columns <- c(
  "Genus",
  "SampleID",
  "Relative_abundance",
  "Progression5",
  "Relative_abundance_percent"
)

missing_columns <- setdiff(
  required_columns,
  colnames(sample_level_raw)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "The input table is missing the following required columns:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 6. Define progression groups
# ==============================================================================

progression_raw_levels <- c(
  "Polyp",
  "UC_remission",
  "UC_active",
  "Dysplasia",
  "CA"
)

progression_display_levels <- c(
  "Polyp",
  "UC remission",
  "UC active",
  "Dysplasia",
  "CAC"
)

progression_display_labels <- c(
  "Polyp" = "Polyp",
  "UC_remission" = "UC remission",
  "UC_active" = "UC active",
  "Dysplasia" = "Dysplasia",
  "CA" = "CAC"
)


# ==============================================================================
# 7. Define representative genera
#
# 7KB taxonomy uses UCG-005.
# Display remains UCG-005.
# ==============================================================================

selected_genera_raw <- c(
  "UCG-005",
  "Desulfovibrio",
  "Lactococcus",
  "Mediterraneibacter"
)

selected_genera_display <- c(
  "UCG-005" = "UCG-005",
  "Desulfovibrio" = "Desulfovibrio",
  "Lactococcus" = "Lactococcus",
  "Mediterraneibacter" = "Mediterraneibacter"
)

selected_genera_display_order <- c(
  "UCG-005",
  "Desulfovibrio",
  "Lactococcus",
  "Mediterraneibacter"
)


# ==============================================================================
# 8. Define fixed group colors
# ==============================================================================

progression_colors <- c(
  "Polyp" = "#7DA9B7",
  "UC remission" = "#8DBA91",
  "UC active" = "#E3A04F",
  "Dysplasia" = "#C77DA8",
  "CAC" = "#A32635"
)


# ==============================================================================
# 9. Clean sample-level data
# ==============================================================================

sample_level <- sample_level_raw %>%
  filter(
    Progression5 %in% progression_raw_levels
  ) %>%
  transmute(
    Genus = as.character(Genus),
    SampleID = as.character(SampleID),
    Progression5 = as.character(Progression5),
    Relative_abundance = as.numeric(Relative_abundance),
    Relative_abundance_percent = as.numeric(Relative_abundance_percent)
  )

if (anyNA(sample_level$Relative_abundance_percent)) {
  stop(
    "Missing relative-abundance percentages were detected.",
    call. = FALSE
  )
}

if (any(sample_level$Relative_abundance_percent < 0)) {
  stop(
    "Negative relative-abundance percentages were detected.",
    call. = FALSE
  )
}


# ==============================================================================
# 10. Strict progression127 audit
# ==============================================================================

sample_count_audit <- sample_level %>%
  distinct(
    SampleID,
    Progression5
  ) %>%
  count(
    Progression5,
    name = "n"
  )

expected_group_counts <- tibble(
  Progression5 = progression_raw_levels,
  expected_n = c(
    26,
    36,
    25,
    17,
    23
  )
)

sample_count_audit <- expected_group_counts %>%
  left_join(
    sample_count_audit,
    by = "Progression5"
  )

if (
  any(is.na(sample_count_audit$n)) ||
  any(sample_count_audit$n != sample_count_audit$expected_n)
) {
  print(
    sample_count_audit,
    n = Inf
  )
  
  stop(
    "The progression127 group counts do not match 26/36/25/17/23.",
    call. = FALSE
  )
}

if (
  dplyr::n_distinct(
    sample_level$SampleID
  ) != 127
) {
  stop(
    "The sample-level table does not contain exactly 127 samples.",
    call. = FALSE
  )
}


# ==============================================================================
# 11. Validate selected genera
# ==============================================================================

available_genera <- unique(
  sample_level$Genus
)

missing_selected_genera <- setdiff(
  selected_genera_raw,
  available_genera
)

if (length(missing_selected_genera) > 0) {
  stop(
    paste0(
      "The following selected genera were not found in the input table:\n",
      paste(
        missing_selected_genera,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 12. Confirm one progression group per sample
# ==============================================================================

sample_metadata <- sample_level %>%
  distinct(
    SampleID,
    Progression5
  )

duplicated_sample_groups <- sample_metadata %>%
  count(
    SampleID,
    name = "Number_of_groups"
  ) %>%
  filter(
    Number_of_groups > 1
  )

if (nrow(duplicated_sample_groups) > 0) {
  stop(
    "One or more samples were assigned to multiple progression groups.",
    call. = FALSE
  )
}


# ==============================================================================
# 13. Confirm one row per selected sample-genus combination
# ==============================================================================

selected_sample_level <- sample_level %>%
  filter(
    Genus %in% selected_genera_raw
  )

duplicate_sample_genus <- selected_sample_level %>%
  count(
    SampleID,
    Genus,
    name = "Number_of_rows"
  ) %>%
  filter(
    Number_of_rows > 1
  )

if (nrow(duplicate_sample_genus) > 0) {
  stop(
    paste0(
      "Duplicated SampleID-Genus combinations were detected.\n",
      "The input table must contain one row per sample and genus."
    ),
    call. = FALSE
  )
}


# ==============================================================================
# 14. Complete missing sample-genus combinations with zero
# ==============================================================================

plot_data <- tidyr::expand_grid(
  sample_metadata,
  Genus = selected_genera_raw
) %>%
  left_join(
    selected_sample_level %>%
      select(
        SampleID,
        Progression5,
        Genus,
        Relative_abundance,
        Relative_abundance_percent
      ),
    by = c(
      "SampleID",
      "Progression5",
      "Genus"
    )
  ) %>%
  mutate(
    Relative_abundance = tidyr::replace_na(
      Relative_abundance,
      0
    ),
    
    Relative_abundance_percent = tidyr::replace_na(
      Relative_abundance_percent,
      0
    ),
    
    Progression_display = factor(
      unname(
        progression_display_labels[
          Progression5
        ]
      ),
      levels = progression_display_levels,
      ordered = TRUE
    ),
    
    Genus_display = factor(
      unname(
        selected_genera_display[
          Genus
        ]
      ),
      levels = selected_genera_display_order,
      ordered = TRUE
    )
  )

if (anyNA(plot_data$Progression_display)) {
  stop(
    "One or more progression-group labels could not be converted.",
    call. = FALSE
  )
}

if (anyNA(plot_data$Genus_display)) {
  stop(
    "One or more genus labels could not be converted.",
    call. = FALSE
  )
}


# ==============================================================================
# 15. Descriptive statistics
#
# IMPORTANT:
# The original final plotting code displays MEDIAN + IQR.
# Mean and SEM are calculated here because the original script exported them,
# but they are NOT the central estimate/error bars drawn in the panel.
# ==============================================================================

descriptive_statistics <- plot_data %>%
  group_by(
    Genus,
    Genus_display,
    Progression5,
    Progression_display
  ) %>%
  summarise(
    Number_of_samples = n(),
    
    Number_detected = sum(
      Relative_abundance_percent > 0
    ),
    
    Detection_rate_percent = 100 * mean(
      Relative_abundance_percent > 0
    ),
    
    Mean_percent = mean(
      Relative_abundance_percent,
      na.rm = TRUE
    ),
    
    SD_percent = stats::sd(
      Relative_abundance_percent,
      na.rm = TRUE
    ),
    
    SEM_percent = SD_percent /
      sqrt(
        Number_of_samples
      ),
    
    Median_percent = stats::median(
      Relative_abundance_percent,
      na.rm = TRUE
    ),
    
    Q1_percent = as.numeric(
      stats::quantile(
        Relative_abundance_percent,
        probs = 0.25,
        na.rm = TRUE,
        names = FALSE
      )
    ),
    
    Q3_percent = as.numeric(
      stats::quantile(
        Relative_abundance_percent,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE
      )
    ),
    
    Minimum_percent = min(
      Relative_abundance_percent,
      na.rm = TRUE
    ),
    
    Maximum_percent = max(
      Relative_abundance_percent,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    IQR_lower = Q1_percent,
    IQR_upper = Q3_percent
  )


# ==============================================================================
# 16. Overall Kruskal-Wallis tests
# ==============================================================================

kruskal_results <- lapply(
  selected_genera_raw,
  function(current_genus) {
    
    current_data <- plot_data %>%
      filter(
        Genus == current_genus
      )
    
    test_result <- stats::kruskal.test(
      Relative_abundance_percent ~ Progression_display,
      data = current_data
    )
    
    tibble(
      Genus = current_genus,
      
      Genus_display = unname(
        selected_genera_display[
          current_genus
        ]
      ),
      
      Number_of_samples = nrow(
        current_data
      ),
      
      Kruskal_Wallis_chi_squared = unname(
        test_result$statistic
      ),
      
      Degrees_of_freedom = unname(
        test_result$parameter
      ),
      
      Kruskal_Wallis_P = test_result$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    BH_adjusted_Q = stats::p.adjust(
      Kruskal_Wallis_P,
      method = "BH"
    )
  )


# ==============================================================================
# 17. P-value formatter
# ==============================================================================

format_panel_p <- function(
    p_value
) {
  
  if (is.na(p_value)) {
    return(
      "Kruskal-Wallis P = NA"
    )
  }
  
  if (p_value < 0.001) {
    return(
      "Kruskal-Wallis P < 0.001"
    )
  }
  
  paste0(
    "Kruskal-Wallis P = ",
    formatC(
      p_value,
      format = "f",
      digits = 3
    )
  )
}


# ==============================================================================
# 18. Function for one representative-genus panel
#
# Kept from the original final plotting script.
# ==============================================================================

create_genus_panel <- function(
    genus_name,
    fixed_y_max = NULL,
    fixed_y_breaks = NULL,
    mark_above_limit = FALSE
) {
  
  current_data <- plot_data %>%
    filter(
      as.character(
        Genus_display
      ) == genus_name
    )
  
  current_summary <- descriptive_statistics %>%
    filter(
      as.character(
        Genus_display
      ) == genus_name
    )
  
  current_test <- kruskal_results %>%
    filter(
      Genus_display == genus_name
    )
  
  if (nrow(current_test) != 1) {
    stop(
      paste0(
        "Could not uniquely identify the Kruskal-Wallis result for ",
        genus_name,
        "."
      ),
      call. = FALSE
    )
  }
  
  current_p_label <- format_panel_p(
    current_test$Kruskal_Wallis_P
  )
  
  if (
    !is.null(fixed_y_max) &&
    isTRUE(mark_above_limit)
  ) {
    
    current_data <- current_data %>%
      mutate(
        Above_display_limit =
          Relative_abundance_percent >
          fixed_y_max,
        
        Plot_abundance = pmin(
          Relative_abundance_percent,
          fixed_y_max
        )
      )
    
    current_summary <- current_summary %>%
      mutate(
        Plot_median = pmin(
          Median_percent,
          fixed_y_max
        ),
        
        Plot_IQR_lower = pmin(
          IQR_lower,
          fixed_y_max
        ),
        
        Plot_IQR_upper = pmin(
          IQR_upper,
          fixed_y_max
        )
      )
    
  } else {
    
    current_data <- current_data %>%
      mutate(
        Above_display_limit = FALSE,
        Plot_abundance =
          Relative_abundance_percent
      )
    
    current_summary <- current_summary %>%
      mutate(
        Plot_median = Median_percent,
        Plot_IQR_lower = IQR_lower,
        Plot_IQR_upper = IQR_upper
      )
  }
  
  ordinary_point_data <- current_data %>%
    filter(
      !Above_display_limit
    )
  
  above_limit_data <- current_data %>%
    filter(
      Above_display_limit
    )
  
  set.seed(
    20260716 +
      match(
        genus_name,
        selected_genera_display_order
      )
  )
  
  current_plot <- ggplot(
    current_data,
    aes(
      x = Progression_display
    )
  ) +
    geom_jitter(
      data = ordinary_point_data,
      aes(
        y = Plot_abundance,
        color = Progression_display
      ),
      width = FIG_JITTER_WIDTH,
      height = 0,
      size = FIG_SMALL_POINT_SIZE,
      alpha = 0.58,
      stroke = 0,
      show.legend = FALSE
    ) +
    geom_errorbar(
      data = current_summary,
      aes(
        x = Progression_display,
        ymin = Plot_IQR_lower,
        ymax = Plot_IQR_upper
      ),
      inherit.aes = FALSE,
      width = 0.18,
      linewidth = FIG_BOX_LINE_WIDTH,
      color = "black"
    ) +
    geom_point(
      data = current_summary,
      aes(
        x = Progression_display,
        y = Plot_median
      ),
      inherit.aes = FALSE,
      shape = 95,
      size = 7,
      stroke = 1,
      color = "black"
    ) +
    scale_color_manual(
      values = progression_colors,
      limits = progression_display_levels,
      drop = FALSE
    ) +
    labs(
      title = genus_name,
      x = NULL,
      y = "Relative abundance (%)"
    ) +
    theme_manuscript() +
    theme(
      plot.title = element_text(
        family = FIG_FONT_FAMILY,
        size = FIG_PLOT_TITLE_SIZE,
        face = "plain",
        hjust = 0.5,
        color = "black",
        margin = margin(
          b = 3
        )
      ),
      
      axis.title = element_text(
        family = FIG_FONT_FAMILY,
        size = FIG_AXIS_TITLE_SIZE,
        face = "bold",
        color = "black"
      ),
      
      axis.text.x = element_text(
        family = FIG_FONT_FAMILY,
        size = FIG_AXIS_TEXT_SIZE,
        angle = 35,
        hjust = 1,
        vjust = 1,
        color = "black"
      ),
      
      axis.text.y = element_text(
        family = FIG_FONT_FAMILY,
        size = FIG_AXIS_TEXT_SIZE,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = FIG_AXIS_LINE_WIDTH,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = FIG_TICK_LINE_WIDTH,
        color = "black"
      ),
      
      legend.position = "none",
      
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 5,
        l = 5
      )
    )
  
  if (is.null(fixed_y_max)) {
    
    observed_maximum <- max(
      c(
        current_data$Relative_abundance_percent,
        current_summary$IQR_upper
      ),
      na.rm = TRUE
    )
    
    annotation_y <- if (
      observed_maximum > 0
    ) {
      observed_maximum * 1.08
    } else {
      0.1
    }
    
    current_plot <- current_plot +
      annotate(
        geom = "text",
        x = 3,
        y = annotation_y,
        label = current_p_label,
        hjust = 0.5,
        vjust = 1,
        size = FIG_SMALL_ANNOTATION_SIZE,
        family = FIG_FONT_FAMILY,
        fontface = "plain",
        color = "black"
      ) +
      scale_y_continuous(
        expand = expansion(
          mult = c(
            0.02,
            0.16
          )
        )
      )
    
  } else {
    
    if (is.null(fixed_y_breaks)) {
      fixed_y_breaks <- pretty(
        c(
          0,
          fixed_y_max
        ),
        n = 5
      )
    }
    
    current_plot <- current_plot +
      annotate(
        geom = "text",
        x = 3,
        y = fixed_y_max * 0.87,
        label = current_p_label,
        hjust = 0.5,
        vjust = 1,
        size = FIG_SMALL_ANNOTATION_SIZE,
        family = FIG_FONT_FAMILY,
        fontface = "plain",
        color = "black"
      ) +
      scale_y_continuous(
        breaks = fixed_y_breaks,
        expand = expansion(
          mult = c(
            0,
            0.03
          )
        )
      ) +
      coord_cartesian(
        ylim = c(
          0,
          fixed_y_max
        ),
        clip = "off"
      )
  }
  
  if (
    !is.null(fixed_y_max) &&
    isTRUE(mark_above_limit) &&
    nrow(above_limit_data) > 0
  ) {
    
    current_plot <- current_plot +
      geom_point(
        data = above_limit_data,
        aes(
          x = Progression_display,
          y = fixed_y_max,
          color = Progression_display
        ),
        inherit.aes = FALSE,
        shape = 17,
        size = FIG_POINT_SIZE,
        alpha = 0.85,
        show.legend = FALSE
      ) +
      geom_text(
        data = above_limit_data,
        aes(
          x = Progression_display,
          y = fixed_y_max,
          label = formatC(
            Relative_abundance_percent,
            format = "f",
            digits = 1
          )
        ),
        inherit.aes = FALSE,
        vjust = -0.50,
        size = 2.2,
        family = FIG_FONT_FAMILY,
        color = "black",
        show.legend = FALSE
      )
  }
  
  current_plot
}


# ==============================================================================
# 19. Create the four original panels
# ==============================================================================

plot_ucg005 <- create_genus_panel(
  genus_name = "UCG-005"
)

plot_desulfovibrio <- create_genus_panel(
  genus_name = "Desulfovibrio",
  fixed_y_max = 2.5,
  fixed_y_breaks = seq(
    0,
    2.5,
    by = 0.5
  ),
  mark_above_limit = TRUE
)

plot_lactococcus <- create_genus_panel(
  genus_name = "Lactococcus"
)

plot_mediterraneibacter <- create_genus_panel(
  genus_name = "Mediterraneibacter"
)


# ==============================================================================
# 20. Combine the original four panels
# ==============================================================================

trajectory_plot <- (
  plot_ucg005 |
    plot_desulfovibrio
) /
  (
    plot_lactococcus |
      plot_mediterraneibacter
  ) +
  patchwork::plot_layout(
    ncol = 2,
    widths = c(
      1,
      1
    ),
    heights = c(
      1,
      1
    )
  )


# ==============================================================================
# 21. Save combined figure
# ==============================================================================

save_manuscript_figure(
  plot_object = trajectory_plot,
  file_stem =
    "Figure3D_representative_genus_4panel_7KB",
  output_dir = output_dir,
  width = 7.2,
  height = 6.2,
  dpi = 600
)


# ==============================================================================
# 22. Save each panel separately
# ==============================================================================

save_manuscript_panel(
  plot_object = plot_ucg005,
  file_stem = "Figure3D_panel_UCG_005_7KB",
  output_dir = output_dir,
  width = 3.6,
  height = 3.1,
  dpi = 600
)

save_manuscript_panel(
  plot_object = plot_desulfovibrio,
  file_stem = "Figure3D_panel_Desulfovibrio_7KB",
  output_dir = output_dir,
  width = 3.6,
  height = 3.1,
  dpi = 600
)

save_manuscript_panel(
  plot_object = plot_lactococcus,
  file_stem = "Figure3D_panel_Lactococcus_7KB",
  output_dir = output_dir,
  width = 3.6,
  height = 3.1,
  dpi = 600
)

save_manuscript_panel(
  plot_object = plot_mediterraneibacter,
  file_stem = "Figure3D_panel_Mediterraneibacter_7KB",
  output_dir = output_dir,
  width = 3.6,
  height = 3.1,
  dpi = 600
)


# ==============================================================================
# 23. Save ggplot objects for final multi-panel assembly
# ==============================================================================

saveRDS(
  plot_ucg005,
  file.path(
    output_dir,
    "Figure3D_panel_UCG_005_7KB.rds"
  )
)

saveRDS(
  plot_desulfovibrio,
  file.path(
    output_dir,
    "Figure3D_panel_Desulfovibrio_7KB.rds"
  )
)

saveRDS(
  plot_lactococcus,
  file.path(
    output_dir,
    "Figure3D_panel_Lactococcus_7KB.rds"
  )
)

saveRDS(
  plot_mediterraneibacter,
  file.path(
    output_dir,
    "Figure3D_panel_Mediterraneibacter_7KB.rds"
  )
)

saveRDS(
  trajectory_plot,
  file.path(
    output_dir,
    "Figure3D_representative_genus_4panel_combined_7KB.rds"
  )
)

message("Saved Figure 3D 7KB ggplot objects as RDS files.")


# ==============================================================================
# 24. Export sample-level plotting data
# ==============================================================================

readr::write_tsv(
  plot_data %>%
    select(
      SampleID,
      Progression5,
      Progression_display,
      Genus,
      Genus_display,
      Relative_abundance,
      Relative_abundance_percent
    ) %>%
    arrange(
      Genus_display,
      Progression_display,
      SampleID
    ),
  file.path(
    output_dir,
    "Figure3D_representative_genus_sample_level_data_7KB.tsv"
  )
)


# ==============================================================================
# 25. Export descriptive statistics
# ==============================================================================

readr::write_tsv(
  descriptive_statistics %>%
    arrange(
      Genus_display,
      Progression_display
    ),
  file.path(
    output_dir,
    "Figure3D_representative_genus_descriptive_statistics_7KB.tsv"
  )
)


# ==============================================================================
# 26. Export Kruskal-Wallis results
# ==============================================================================

readr::write_tsv(
  kruskal_results,
  file.path(
    output_dir,
    "Figure3D_representative_genus_Kruskal_Wallis_results_7KB.tsv"
  )
)


# ==============================================================================
# 27. Export Desulfovibrio display-limit audit
#
# The actual display ceiling in the original plotting code is 2.5%.
# ==============================================================================

desulfovibrio_display_audit <- plot_data %>%
  filter(
    as.character(
      Genus_display
    ) == "Desulfovibrio"
  ) %>%
  mutate(
    Above_2_5_percent_display_limit =
      Relative_abundance_percent > 2.5
  ) %>%
  arrange(
    desc(
      Relative_abundance_percent
    )
  )

readr::write_tsv(
  desulfovibrio_display_audit,
  file.path(
    output_dir,
    "Figure3D_Desulfovibrio_display_limit_audit_7KB.tsv"
  )
)


# ==============================================================================
# 28. Save session information
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "sessionInfo_Figure3D_representative_genus_trends_7KB.txt"
  )
)


# ==============================================================================
# 29. Console summary
# ==============================================================================

message("")
message("==============================================================")
message("7KB representative genus abundance panels completed.")
message("==============================================================")
message("")
message("Progression127 group counts:")
print(
  sample_count_audit,
  n = Inf
)

message("")
message("Kruskal-Wallis results:")
print(
  kruskal_results,
  n = Inf,
  width = Inf
)

message("")
message("Input file:")
message(input_file)

message("")
message("Output directory:")
message(output_dir)

message("")
message("Actual display used by original plotting code:")
message("  Points: individual samples")
message("  Central estimate: median")
message("  Error bars: IQR")
message("  Test: overall Kruskal-Wallis")
message("")
message("Final manuscript assembly:")
message("  Keep Mediterraneibacter")
message("  Keep UCG-005")
message("  Keep Lactococcus")
message("  Desulfovibrio remains generated but is not used in the final combined figure.")