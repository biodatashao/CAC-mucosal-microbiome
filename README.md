# Mucosal microbiome in ulcerative colitis-associated colorectal carcinogenesis

Analysis code for the manuscript:

> **[MANUSCRIPT TITLE]**
> [Author list]
> *Gut Microbes*, [year]. DOI: [to be added on acceptance]

This repository contains the complete R analysis pipeline used to generate every
figure, supplementary figure and reported statistic in the manuscript, based on
16S rRNA gene sequencing of 150 FFPE colorectal mucosal specimens across five
clinical groups.

---

## Data availability

| Resource | Location |
|---|---|
| Raw 16S sequencing reads (300 paired-end FASTQ files) | NCBI SRA, BioProject **[PRJNAxxxxxx]** |
| Sample-level sequencing metadata | NCBI BioSample, linked to the BioProject above |
| Clinical metadata (patient-level) | Not publicly deposited — see below |
| Analysis code | This repository |

Patient-level clinical metadata (disease duration, Mayo score, AJCC stage,
recurrence status, disease-free survival, ESR/CRP) are not deposited publicly
because of patient-privacy restrictions and the terms of the ethics approval
granted by the Institutional Review Board of Peking Union Medical College
Hospital ([approval number]). De-identified clinical variables required to
reproduce specific analyses are available from the corresponding author on
reasonable request.

---

## Repository layout

```
.
├── config.R                     # single place where all paths are defined
├── environment/                 # session info, package and QIIME 2 versions
├── data/example/                # small de-identified example input tables
├── docs/
│   ├── script_figure_map.md     # which script produces which figure/table
│   └── script_rename_map.csv    # original filename -> current filename
└── scripts/
    ├── 01_data_preparation/
    ├── 02_diversity_Figure2/
    ├── 03_taxa_LEfSe_Figure3/
    ├── 04_DMM_Figure4/
    ├── 05_clinical_recurrence_Figure7/
    ├── 06_supp_inflammation_signature/
    ├── 07_sensitivity_analyses/
    └── 08_PICRUSt2_Figure5/
```

Scripts are numbered `MM_NN_description.R`, where `MM` is the module and `NN` is
the execution order within that module. Modules must be run in order: module 01
builds the clean tables that every later module reads.

---

## Reproducing the analysis

### 1. Obtain the sequencing data

```bash
# Requires the SRA Toolkit: https://github.com/ncbi/sra-tools
prefetch PRJNAxxxxxx
fasterq-dump --split-files SRRxxxxxxx
```

### 2. Upstream processing

Raw reads were processed into an ASV table using QIIME 2 v2025.4 with DADA2 and
classified against SILVA 138.2. Reagent and laboratory contaminants were removed
with `decontam` (prevalence method, threshold 0.5) using seven negative controls.
Functional potential was predicted with PICRUSt2, and pathway- and genus-level
biomarker discovery used LEfSe (Kruskal-Wallis alpha 0.05, LDA > 2.0); both are
run at the command line, and the exact invocations are printed by the
corresponding preparation scripts in modules 03, 06 and 08.

Full software and package versions are listed in
[`environment/README.md`](environment/README.md). Reagent and laboratory contaminants were removed using `decontam` against
seven negative controls. The resulting ASV count table, taxonomy table and
sample metadata are the inputs to this repository.

Scripts in `07_sensitivity_analyses/` characterise the negative controls and
test the robustness of the main findings.

### 3. Configure paths

All scripts resolve paths through `config.R`. Set the project root once:

```bash
export FFPE_PROJECT_ROOT=/path/to/your/project
```

Expected directory layout under that root:

```
$FFPE_PROJECT_ROOT/
├── raw/                                     # raw sequences, taxonomy assignments
├── data/clinical/                           # clinical metadata (not distributed)
└── output/
    ├── decontamination_8controls/           # first decontamination run
    ├── decontamination_7controls_final/     # decontamination used in the manuscript
    └── analysis/                            # all analysis outputs (modules 01-08)
```

Two decontamination directories are listed because the study used seven negative
controls in its final decontamination. An earlier run had included an eighth
control (`kb2`) that failed quality checks and was excluded. The earlier run is
retained only because module 01 reads its metadata to reconstruct the exact
127-sample progression cohort; no result in the manuscript derives from it.

### 4. Run the pipeline

```bash
Rscript scripts/01_data_preparation/01_01_build_master_clean_tables.R
Rscript scripts/01_data_preparation/01_02_build_paired_CAC23_nonCAC23_cohort.R
# ... then modules 02 through 07 in order
```

Each script is self-contained: it clears the workspace, sources `config.R`,
reads its inputs from disk and writes its outputs to disk.

---

## Analysis overview

| Module | Content |
|---|---|
| **01** | Build the master decontaminated ASV/taxonomy/metadata tables; reconstruct the 127-sample progression cohort and the paired 23 CAC vs 23 nonCAC cohort |
| **02** | Alpha diversity (Shannon, observed genera) and beta diversity (Bray-Curtis, PCoA, PERMANOVA) → **Figure 2**, **Supplementary Figure 1** |
| **03** | Genus-level composition and LEfSe biomarker discovery across the five progression groups and in paired CAC vs nonCAC → **Figure 3** |
| **04** | Dirichlet multinomial mixture (DMM) probabilistic community-state modelling: model selection over K = 1–7, characterisation of the K = 3 states, and projection of the paired cohort onto the fitted model → **Figure 4** |
| **05** | Clinical association analyses in CAC: diversity by AJCC stage, dysplasia grade and recurrence; Bray-Curtis PCoA by recurrence; disease-free survival; and the two-genus Recurrence-associated Microbial Index (RMI) → **Figure 7**, **Supplementary Figure 3** |
| **06** | Supplementary inflammation-spectrum analysis (polyp vs UC remission vs active UC), oral-associated microbial score and UCG-005 versus ESR/CRP → **Supplementary Figure [x]** |
| **07** | Sensitivity analyses (monotonic trend testing, rarefaction-depth robustness of the DMM solution, effect sizes) and negative-control characterisation → **Supplementary Tables** |
| **08** | PICRUSt2 MetaCyc pathway prediction and pathway-level LEfSe, across the five progression groups and in paired CAC vs nonCAC → **Figure 5**, **Supplementary Figure 2** |

A per-script table is in [`docs/script_figure_map.md`](docs/script_figure_map.md).

---

## Reproducibility notes

- Analyses were run in R 4.5.1; package versions are in `environment/README.md`.
- Random seeds are set in `config.R` and, where the analysis requires it, inside
  the individual scripts (DMM fitting, PERMANOVA, rarefaction, bootstrap).
- The DMM model is fitted **once** in `04_01`; all downstream DMM scripts read
  the saved model object rather than refitting it.
- Scripts in module 07 are read-only: they test robustness and characterise the
  negative controls from the frozen source tables, and do not modify any
  analysis output.
- All modules, including the PICRUSt2 functional prediction in module 08, read
  from the same 7-negative-control decontaminated tables.

---

## Citation

If you use this code, please cite the manuscript above.

## Contact

[Name], [email] — Department of Gastroenterology, Peking Union Medical College
Hospital, Chinese Academy of Medical Sciences and Peking Union Medical College,
Beijing, China.
