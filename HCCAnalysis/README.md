# 🧬 HCCAnalysis: Factor IV Microbiome–Metabolome Causal Pipeline

Workflow for the MD Anderson HCC mouse cohort, integrating gut microbiome and metabolomics profiles to infer latent factors and causal links to liver phenotypes using the **Factor IV (GOFAR + IV-probit)** framework. Run from the **repository root** so the hard-coded `HCCAnalysis/...` paths resolve correctly.

---

## Folder Layout

```
HCCAnalysis/
├── data/MDAnderson/
│   ├── 25970_otu_table.biom                # microbiome OTUs (used)
│   ├── sample_information_from_prep_2458.tsv
│   ├── 2914e76aee36462c8fcfe11a2b244917.csv # metabolomics intensities
│   ├── gnps_metabolite_annotations.tsv      # GNPS IDs → names/classes
│   └── processed_data_clref.rds             # output of dataPrep_md_otu.R
│       processed_data_metabolites_otu.rds   # output of dataPrep_metabolitesGNP.R
│       microbiome_metabolites_data.xlsx     # normalized matrices export
│       ... (other intermediate RDS files)
│
│── dataPrep_md_otu.R            # OTU/metadata cleaning, instruments/controls/phenotypes
│── dataPrep_metabolitesGNP.R    # metabolite processing + microbiome alignment, CLR, export
│── microb_metabolites_otu.R     # GOFAR fit (prefit load), IV-probit, figures, Excel outputs
│
├── output/ 
│   ├── microb_metabolites_full_otu_wsd1.RData   # prefit GOFAR object loaded by the analysis script
│   ├── microb_metabolites_analysis.xlsx         # U/V/C matrices, causal effects, annotations
│   ├── causal_effect_mtb.xlsx                   # metabolite-only causal effects
│   ├── metabolites_annotation.xlsx              # GNPS annotations with compound types
│   ├── microbe_metabolites_analysis_otu.txt     # IV-probit console log
│   └── deseq_microbiome.csv                     # DESeq2 results (created when DESeq section runs)
│
├── plots/
│   ├── balance_v1.pdf, balance_v2.pdf           # latent balance boxplots
│   ├── lf_effect.pdf, lf_effect2.pdf            # effect-size and log-p heatmaps (IV-probit)
│   ├── rel_imp_intervention.pdf                 # instrument loadings (U)
│   └── causal_effect_mtb_mb.pdf                 # combined metabolite/taxa causal map
└── README_HCC.md
```
[For data in the output folder please use the Zenodo link of the manuscript. ](https://zenodo.org/uploads/17968387)

---

## Workflow (run from repo root)

1) **Prep microbiome OTUs and metadata**  
   `Rscript HCCAnalysis/scripts/dataPrep_md_otu.R`  
   - Imports `25970_otu_table.biom` and `sample_information_from_prep_2458.tsv`  
   - Filters low-depth samples (>10k reads) and cleans taxonomy  
   - Builds **instruments (Z)** from treatment/diet/antibiotic indicators and **controls (Q)** from age, weight, strain  
   - Encodes phenotypes (HCC, NASH, immune perturbations)  
   - Saves `processed_data_clref.rds` in `data/MDAnderson/`

2) **Integrate metabolites + microbiome**  
   `Rscript HCCAnalysis/scripts/dataPrep_metabolitesGNP.R`  
   - Reads metabolomics table `2914e76aee36462c8fcfe11a2b244917.csv`, aligns to host IDs  
   - Filters low-variance metabolites, applies CLR to metabolites and OTUs  
   - Harmonizes samples across Z/Q/phenotypes/microbiome/metabolites  
   - Exports `processed_data_metabolites_otu.rds` and `microbiome_metabolites_data.xlsx`

3) **Factor IV modeling, IV-probit, and figures**  
   `Rscript HCCAnalysis/scripts/microb_metabolites_otu.R`  
   - Loads `processed_data_metabolites_otu.rds` and the **prefit GOFAR object** `output/microb_metabolites_full_otu_wsd1.RData`  
   - If you need to refit, uncomment the `gofar_s` block and comment out the `load(ofile)` line  
   - Constructs latent scores (U_hat, IV_hat), runs IV-probit per phenotype, and computes causal effects  
   - Generates plots (`plots/`) and Excel summaries (`output/`)

4) **Differential abundance (optional, end of script)**  
   The DESeq2 section in `microb_metabolites_otu.R` runs after the causal analysis and writes `output/deseq_microbiome.csv`.

---

## Key Outputs

| File | Description |
|------|-------------|
| `output/microb_metabolites_full_otu_wsd1.RData` | Prefit GOFAR model (`fit.seq`) used by the pipeline |
| `output/microb_metabolites_analysis.xlsx` | U (instruments), V (features), C (covariates), causal effects, annotations |
| `output/causal_effect_mtb.xlsx` | Metabolite-only causal effect matrix with GNPS annotations |
| `output/metabolites_annotation.xlsx` | GNPS hits with compound class labels |
| `plots/causal_effect_mtb_mb.pdf` | Joint metabolite–taxa causal heatmap (selected bile acids, antibiotics, *R. gnavus*, Enterobacteriaceae) |
| `plots/lf_effect*.pdf`, `plots/rel_imp_intervention.pdf` | IV-probit effect summaries and instrument loadings |
| `plots/balance_v*.pdf` | Latent balance boxplots across liver phenotypes |
| `output/deseq_microbiome.csv` | DESeq2 differential abundance results (if run) |

---

## Factor IV Notes

- **Instruments (Z):** treatment/diet/antibiotic indicators plus parsed group labels.  
- **Controls (Q):** age, weight, strain.  
- **Responses (Y):** CLR-transformed metabolites concatenated with CLR-transformed OTUs.  
- **Model:** GOFAR estimates low-rank loadings (U, V, C), then latent scores feed an **IV-probit** to recover causal effects on phenotypes.  
- **Significance:** coefficients with p < 0.05 are retained when forming causal effects; heatmaps show normalized effects.

---

## Dependencies (R ≥ 4.2)

`phyloseq`, `microbiome`, `gofar`, `nbfar`, `ivprobit`, `DESeq2`,  
`tidyverse`, `tidyr`, `tidyselect`, `magrittr`, `ggplot2`, `ggpubr`, `rstatix`,  
`ComplexHeatmap`, `pheatmap`, `circlize`, `grid`, `RColorBrewer`,  
`readxl`, `openxlsx`, `biomformat`, `qdapRegex`

Install example:

```r
install.packages(c(
  "tidyverse","phyloseq","microbiome","ggplot2","ggpubr","rstatix",
  "ComplexHeatmap","pheatmap","circlize","grid","RColorBrewer",
  "readxl","openxlsx","biomformat","qdapRegex","DESeq2","ivprobit"
))
devtools::install_github("amishra-stats/gofar")
devtools::install_github("amishra-stats/nbfar")
```

