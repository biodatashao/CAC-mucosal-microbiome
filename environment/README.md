# Environment

## R

Analyses were performed in **R version 4.5.1** with **Bioconductor 3.22**.

Package versions recorded from the analysis environment:

| Package | Version | Used for |
|---|---|---|
| phyloseq | 1.54.2 | community data container, alpha diversity |
| vegan | 2.7.5 | Bray-Curtis dissimilarity, PCoA, PERMANOVA, betadisper, rarefaction |
| DirichletMultinomial | 1.52.0 | Dirichlet multinomial mixture community-state modelling |
| compositions | 2.0.9 | centred log-ratio transformation |
| clinfun | 1.1.5 | Jonckheere-Terpstra trend test |
| survival | 3.8.6 | Kaplan-Meier estimates, log-rank test |
| pROC | 1.19.0.1 | ROC analysis, AUC with bootstrap confidence intervals |
| tidyverse | 2.0.0 | data handling (meta-package) |
| dplyr | 1.2.1 | data handling |
| tidyr | 1.3.2 | data handling |
| readr | 2.2.0 | delimited file I/O |
| tibble | 3.3.1 | data handling |
| stringr | 1.6.0 | string handling |
| forcats | 1.0.1 | factor handling |
| data.table | 1.18.4 | large table I/O and aggregation |
| readxl | 1.5.0 | clinical metadata import |
| ggplot2 | 4.0.3 | figures |
| patchwork | 1.3.2 | figure assembly |
| cowplot | 1.2.0 | figure assembly |
| pheatmap | 1.0.13 | heatmaps |
| ggalluvial | 0.12.6 | alluvial state-transition plots |
| RColorBrewer | 1.1.3 | colour palettes |
| scales | 1.4.0 | axis scaling and formatting |
| ragg | 1.5.2 | high-resolution TIFF export |
| rprojroot | 2.1.1 | path resolution in `config.R` |

To record the full session state, run in R and save the output to
`sessionInfo.txt`:

```r
sessionInfo()
```

### Statistics implemented directly in the scripts

The following were computed from first principles rather than through a
package, so no version applies:

- **Cliff's delta** — all pairwise comparisons between groups
  (`07_09`, `07_10`, `07_11`)
- **Adjusted Rand index** — for DMM cluster stability across rarefaction
  depths (`07_07`, `07_08`)
- **Recurrence-associated Microbial Index (RMI)** — standardised
  log-abundance difference with leave-one-out cross-validation (`05_03`)

## Upstream processing

Steps from raw reads to the ASV table were performed by the sequencing provider.
Versions below are as reported by the provider; the corresponding software is not
part of the analysis environment used for this repository.

| Step | Software | Version |
|---|---|---|
| Read merging | FLASH | v1.2.11 |
| Reverse-primer trimming | Cutadapt | v3.3 |
| Quality filtering | fastp | v0.23.1 |
| Chimera removal | vsearch | v2.16.0 |
| ASV inference | DADA2 | as implemented in QIIME 2 v2025.4 |
| Pipeline | QIIME 2 | v2025.4 |
| Reference database | SILVA | 138.2 |


Software run locally for the analyses in this repository:

| Step | Software | Version |
|---|---|---|
| Contaminant identification | decontam (R/Bioconductor) | 1.30.0 |
| Functional prediction | PICRUSt2 | 2.6.3 |
| Biomarker discovery | LEfSe | 1.1.2 (Bioconda) |
| Additional plotting | GraphPad Prism | 10.5.0 for macOS |

Key parameters:

- **decontam**: prevalence-based method, threshold = 0.5, using seven negative
  controls (tissue-free FFPE blocks and reagent-only blanks) as reference
- **LEfSe**: Kruskal-Wallis alpha = 0.05, LDA score threshold > 2.0
- **Rarefaction sensitivity analysis**: minimum library size 25,033 reads

To capture the local conda environments in full:

```bash
conda env export -n picrust2 > environment/conda_picrust2.yml
conda env export -n lefse    > environment/conda_lefse.yml
```

## Random seeds

`config.R` sets a global seed (`GLOBAL_SEED <- 20260726`). Scripts performing
permutation, bootstrap, rarefaction or DMM fitting set their own seed where the
analysis requires it.
