# 🧬 DIABIMMUNE Microbiome Causal Analysis Workflow

## Overview

This repository contains the full reproducible workflow for the **Diabimmune cohort analysis**, examining **early-life microbiome development** and its **causal relationships with autoimmune risk (e.g., Type 1 Diabetes)**.  
The analysis integrates metagenomic and 16S rRNA sequencing data with participant metadata, dietary exposures, and clinical outcomes using **multi-task learning and instrumental variable (IV) modeling**.

All analyses are performed in **R (≥ 4.3)** and rely on open-source packages (full list below).

---

## Folder Structure

```
Diabimmune/
│
├── data/
│   └── data_diabimmune/
│       ├── DiabImmune_Sample.txt
│       ├── DiabImmune_Participant.txt
│       ├── DiabImmune_participant_repeated_measure.txt
│       ├── DiabImmune_Metagenomic_sequencing_assay.txt
│       ├── DiabImmune_16S_rRNA_(V4)_assay.txt
│       ├── DiabImmune_OntologyMetadata.txt
│       └── processed_data_diabiimune.rds
│
├── output/
│   ├── microb_diabimmune_mtg_gofar*.RData
│   ├── microb_diabimmune_16s_gofar*.RData
│   └── causal_effect_*_gofar*.RData
│
├── plots/
│   ├── causal_heatmap_metagenomics_diabimmune.pdf
│   ├── causal_barplot_metagenomics_diabimmune.pdf
│   └── causal_boxplot_metagenomics_diabimmune.pdf
│
├── dataPrep_diabimmune.R
├── microb_diabimmune.R
└── compile_figure_diabimmune.R
```
[For data in the output folder please use the Zenodo link of the manuscript. ](https://zenodo.org/uploads/17968387)


---

## Data Source and Context

The **DIABIMMUNE project** investigates how infant gut microbiome composition and function influence the development of autoimmune diseases, including **Type 1 Diabetes (T1D)**.  It includes multi-omic measurements (metagenomics, 16S rRNA), repeated longitudinal sampling, and comprehensive participant metadata (e.g., birth mode, diet, antibiotics, growth metrics).

Each sample includes:
- **Microbiome data:** metagenomic species-level abundances and 16S operational taxonomic units (OTUs).  
- **Participant metadata:** perinatal, anthropometric, and dietary features.  
- **Clinical outcomes:** autoantibody status, HLA risk, and diabetes diagnosis timing.

---

## Analysis Workflow

The pipeline is divided into **three main R scripts**, run sequentially.

### 1️⃣ Data Processing — `dataPrep_diabimmune.R`

**Purpose:**  
To clean, harmonize, and integrate all raw DIABIMMUNE data files into a structured R object for downstream modeling.

**Main steps:**
- Load and sanitize sample, participant, and ontology metadata.
- Handle missing values (mean/mode imputation).
- Filter redundant or invariant variables.
- Match metagenomic and 16S samples with participant data.
- Generate per-participant microbiome matrices across age windows:
  - 0–12 months
  - 12–24 months
  - Greater than 24 months
- Save processed objects: `data/data_diabimmune/processed_data_diabiimune.rds`

**Output:** Harmonized microbiome matrices, participant covariates, and outcomes for each time window.

---

### 2️⃣ Model Fitting — `microb_diabimmune.R`

**Purpose:**  
Perform multi-task learning and causal modeling of microbiome–phenotype associations using **GOFAR** and **IV-Probit regression**.

**Main steps:**
1. Centered log-ratio (CLR) normalization of microbiome abundances.
2. Standardization of participant covariates and instruments.
3. Fit low-rank multi-task regression (GOFAR) linking microbiome composition to autoimmune outcomes.
4. Use latent factors as valid instruments in IV-Probit causal inference framework.

**Outputs:**
- GOFAR model fits: `output/microb_diabimmune_mtg_gofar*.RData`, `output/microb_diabimmune_16s_gofar*.RData`
- Causal effect estimates: `output/causal_effect_mtg_gofar*.RData`, `output/causal_effect_16s_gofar*.RData`

---

### 3️⃣ Figure Compilation — `compile_figure_diabimmune.R`

**Purpose:**  
Compile and visualize model outputs for manuscript figures.

**Main Figures:**
1. Heatmap of normalized causal effects of microbial species on autoimmune phenotypes.  
2. Bar plot of top early-life instruments contributing to latent microbial factors.  
3. Box plot of latent factor scores stratified by autoimmune outcomes.  
4. Summary table linking microbial factors with autoimmune outcomes across age windows.

**All figures are saved under:** `Diabimmune/output/plots/`

---

## Reproducibility Notes

1. Analyses are deterministic (`set.seed(123)`).
2. Run from the **repository root** (paths are hard-coded with `Diabimmune/...` prefixes). Running from inside `Diabimmune/` will break file lookups.
3. Raw data must be stored under `Diabimmune/data/data_diabimmune/`.
4. Generate results for each age window by looping `index_select` (1: 0–12 mo, 2: 12–24 mo, 3: >24 mo). Example:


5. `microb_diabimmune.R` saves causal-effect objects to `Diabimmune/output/`. Before running `compile_figure_diabimmune.R`, ensure these files are accessible at `Diabimmune/data/data_diabimmune/` (the script loads from there). 

6. Figures are written to `Diabimmune/plots/`.

---

## Key Dependencies

| Package | Purpose |
|----------|----------|
| phyloseq, microbiome | Microbial data manipulation |
| compositions | CLR transformation |
| gofar | Low-rank multitask regression |
| ivprobit | Instrumental variable probit regression |
| ggplot2, ComplexHeatmap, pheatmap, ggpubr | Visualization |
| tidyverse, dplyr, tidyr, stringr, magrittr | Data manipulation and pipes |
| DESeq2 | Differential abundance utilities |
| qdapRegex, biomformat, readxl | Data cleaning/IO helpers |
