# Script → figure / table map

Run modules in numerical order. Within a module, run scripts in numerical order.

Legend: **source** = computes and writes a source table; **plot** = renders a
panel; **assemble** = combines saved panels into the final figure file.

---

## Module 01 — Data preparation

| Script | Role | Produces |
|---|---|---|
| `01_01_build_master_clean_tables.R` | source | Master decontaminated ASV / taxonomy / metadata tables; reconstructed 127-sample progression cohort; candidate CAC / nonCAC pool |
| `01_02_build_paired_CAC23_nonCAC23_cohort.R` | source | Locked paired 23 CAC vs 23 nonCAC cohort + pair manifest |

## Module 02 — Diversity (Figure 2, Supplementary Figure 1)

| Script | Role | Produces |
|---|---|---|
| `02_01_alpha_beta_statistics.R` | source | Alpha/beta diversity statistics for progression127 and the paired cohort; Supplementary Table S4 |
| `02_02_plot_Figure2_and_SuppFigure1.R` | plot + assemble | **Figure 2**, **Supplementary Figure 1** |

## Module 03 — Composition and LEfSe (Figure 3)

| Script | Role | Produces |
|---|---|---|
| `03_01_genus_composition_5groups_source.R` | source | Top-20 genus composition, five progression groups |
| `genus_color_mapping.R` | helper | Shared genus-to-colour mapping, sourced by `03_02` and `03_04` |
| `03_02_plot_Figure3A_genus_composition_5groups.R` | plot | **Figure 3A** |
| `03_03_genus_composition_CAC_vs_nonCAC_source.R` | source | Top-20 genus composition, paired CAC vs nonCAC |
| `03_04_plot_Figure3B_genus_composition_CAC_vs_nonCAC.R` | plot | **Figure 3B** |
| `03_05_LEfSe_progression5_genus.R` | source | Genus-level LEfSe, five groups |
| `03_06_plot_Figure3_LEfSe_progression5_markers.R` | plot | LDA barplot + row-scaled heatmap panels |
| `03_07_LEfSe_CAC_vs_nonCAC_genus.R` | source | Genus-level LEfSe, CAC vs nonCAC (unpaired) |
| `03_08_plot_Figure3C_LEfSe_CAC_vs_nonCAC_markers.R` | plot | **Figure 3C** |
| `03_09_Figure3_representative_genus_trends.R` | plot | Representative genus trend panels |
| `03_10_assemble_Figure3_8panels.R` | assemble | **Figure 3** (8 panels) |

## Module 04 — DMM community states (Figure 4)

| Script | Role | Produces |
|---|---|---|
| `04_01_DMM_fit_K1toK7_progression127.R` | source | DMM fit, K = 1–7 × 20 repeats; best-converged model objects, assignments, posteriors |
| `04_02_DMM_characterize_K3_states.R` | source | C1/C2/C3 characterisation, dominant genera, association with progression stage |
| `04_03_DMM_cluster_mean_heatmap_source.R` | source | State-defining genus selection, row-normalised cluster-mean matrix |
| `04_04_DMM_project_paired_cohort_to_K3.R` | source | Projection of paired nonCAC/CAC samples onto the fitted K = 3 model (no refit) |
| `04_05_plot_Figure4A_Laplace_model_selection.R` | plot | **Figure 4A** |
| `04_06_plot_Figure4B_state_distribution.R` | plot | **Figure 4B** |
| `04_07_plot_Figure4C_state_heatmap.R` | plot | **Figure 4C** |
| `04_08_plot_Figure4D_ecological_metrics.R` | plot | **Figure 4D** |
| `04_09_plot_Figure4E_paired_state_correspondence.R` | plot | **Figure 4E** |
| `04_10_assemble_Figure4.R` | assemble | **Figure 4** (A–E) |

> **Note.** This script reruns the panel code inline rather than reading the
> saved panel objects, so `04_05`–`04_09` are not strictly required to
> regenerate Figure 4. They are retained because they also write the panel
> source tables.

## Module 05 — Clinical association and recurrence (Figure 7, Supplementary Figure 3)

| Script | Role | Produces |
|---|---|---|
| `05_01_lock_clinical_metadata.R` | source | Locked CAC clinical metadata: AJCC stage, recurrence, DFS status and time |
| `05_02_core_analysis_Figure7_SuppFig3.R` | source | All core statistics and source tables for Figure 7 and Supplementary Figure 3 |
| `05_03_Figure7DHI_RMI_original_method.R` | source | Figure 7D/H/I recomputed with the original methodology: Bray-Curtis PCoA, PERMANOVA (9999 permutations), betadisper, stage-adjusted PERMANOVA, RMI |
| `05_04_assemble_Figure7.R` | assemble | **Figure 7** (A–I) |
| `05_05_assemble_SupplementaryFigure3.R` | assemble | **Supplementary Figure 3** |

## Module 06 — Supplementary inflammation signature

| Script | Role | Produces |
|---|---|---|
| `06_01_inflammation3_alpha_beta_statistics.R` | source | Alpha/beta diversity: polyp vs UC remission vs active UC |
| `06_02_inflammation3_alpha_beta_plot.R` | plot | Diversity panels |
| `06_03_LEfSe_inflammation3_genus.R` | source | Genus-level LEfSe across the inflammation spectrum |
| `06_04_plot_inflammation3_LEfSe_panels.R` | plot | Panels A (LDA barplot) and B (heatmap) |
| `06_05_build_CLR_clinical_merged_dataset.R` | source | CLR-transformed abundance merged with clinical variables |
| `06_06_oral_score_vs_ESR_CRP.R` | plot | Panels C and D |
| `06_07_UCG005_vs_ESR_CRP.R` | plot | Panels E and F |
| `06_08_assemble_SuppFig_inflammation_signature.R` | assemble | **Supplementary Figure [x]** |

## Module 07 — Sensitivity analyses and negative controls

Read-only. These scripts test the robustness of the main findings and
characterise the negative controls; they do not modify any analysis output.

| Script | Produces |
|---|---|
| `07_01_JonckheereTerpstra_trend_top20_genera.R` | Monotonic trend test across progression groups, BH-corrected over the top 20 genera |
| `07_02_DMM_rarefaction_sensitivity_ARI.R` | DMM cluster stability across rarefaction depths (adjusted Rand index) |
| `07_03_DMM_rarefaction_K3_multirepeat_ARI.R` | K = 3 stability across repeated rarefactions |
| `07_04_sequencing_depth_rarefaction_sensitivity.R` | Sequencing depth summary and post-rarefaction alpha diversity → **rarefaction supplementary table** |
| `07_05_Figure7H_Cliffs_delta.R` | Cliff's delta for the Figure 7H recurrence comparison |
| `07_06_Figure7H_RMI_Cliffs_delta_CI.R` | Cliff's delta with confidence interval for the RMI |
| `07_07_negative_control_read_depth.R` | Read depth and source of each negative control |
| `07_08_negative_control_raw_composition.R` | Negative-control composition, raw reads |
| `07_09_negative_control_predecontam_composition.R` | Negative-control composition before decontamination |
| `07_10_final_negative_control_table.R` | **Negative-control supplementary table** |

## Module 08 — PICRUSt2 functional pathways (Figure 5, Supplementary Figure 2)

MetaCyc pathway abundances were predicted with PICRUSt2 and compared with LEfSe.
PICRUSt2 and LEfSe themselves are run at the command line; the R scripts below
prepare their inputs and process their outputs. Each preparation script prints
the exact command to run at the end of its execution.

| Script | Role | Produces |
|---|---|---|
| `08_01_locate_PICRUSt2_pathway_files.R` | utility | Locates PICRUSt2 pathway abundance files and checks sample-name concordance with the metadata |
| `08_02_prepare_PICRUSt2_input_progression127.R` | source | Feature table, representative-sequence FASTA and metadata for the 127-sample progression cohort → **run PICRUSt2** |
| `08_03_prepare_LEfSe_input_pathways_progression5.R` | source | CPM-normalised MetaCyc pathway table, prevalence filtering, LEfSe input → **run LEfSe** |
| `08_04_parse_and_select_pathways_progression5.R` | source | Parses the LEfSe `.res`, maps pathway descriptions and modules, selects representative pathways per group |
| `08_05_plot_pathway_LEfSe_progression5_panels.R` | plot | Main panel source (UC active / dysplasia / CAC) and supplementary panel source (polyp / UC remission) |
| `08_06_prepare_PICRUSt2_input_CAC_vs_nonCAC.R` | source | PICRUSt2 input for the paired CAC vs nonCAC cohort → **run PICRUSt2** |
| `08_07_prepare_LEfSe_input_pathways_CAC_vs_nonCAC.R` | source | LEfSe input for CAC vs nonCAC → **run LEfSe** |
| `08_08_parse_LEfSe_pathways_CAC_vs_nonCAC.R` | source | Parses the CAC vs nonCAC `.res` into clean marker tables |
| `08_09_plot_pathway_LEfSe_CAC_vs_nonCAC_selected.R` | plot | Selects top pathways per group/module, removes artefactual pathways, writes the panel source |
| `08_10_assemble_Figure5.R` | assemble | **Figure 5A** (progression groups) + **Figure 5B** (CAC vs nonCAC) |
| `08_11_assemble_SupplementaryFigure2.R` | assemble | **Supplementary Figure 2** (polyp vs UC remission) |

`08_PICRUSt2_Figure5/superseded/` contains three earlier scripts whose outputs
are not read by any script in the final pipeline. See the README in that folder.

> **Dataset version.** Module 08 was rerun on the same 7-negative-control
> decontaminated tables used by modules 01–07 (`output/final_7KB_rerun_20260726`).
> The pathway-level results were unchanged from the earlier decontamination
> version. All module-08 inputs, outputs and filenames now carry the `7KB` /
> `progression127` naming used elsewhere in the repository.

---

## Not covered by this repository

The following manuscript items have no script in this repository. Add them, or
state in the manuscript how they were produced:

- **Figure 1** and **Figure 6** — not produced in R; state in the manuscript or
  figure legends how they were generated
- **Upstream QIIME 2 / DADA2 / decontam processing** — no scripts present;
  document the software versions and parameters in `environment/`
- **PICRUSt2 and LEfSe command-line invocations** — run outside R; the exact
  commands are printed by the corresponding preparation scripts and should be
  recorded in `environment/`
