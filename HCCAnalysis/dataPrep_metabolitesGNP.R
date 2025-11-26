# ================================================================
# File: dataPrep_metabolitesGNP.R
# Project: FactorIV - HCCAnalysis
# Purpose: Integrate metabolomics and microbiome data (OTU) for 
#          FactorIV causal inference in HCC mouse study.
#
# This script:
#   1. Loads processed microbiome data from dataPrep_md_otu.R
#   2. Reads metabolite abundance tables
#   3. Harmonizes sample IDs across data types
#   4. Filters low-variance metabolites
#   5. Applies log-centered normalization (CLR)
#   6. Constructs joint metabolite–microbiome matrices
#   7. Saves combined data as .RDS and .XLSX for further modeling
#
# Dependencies: phyloseq, microbiome, magrittr, tidyverse, ComplexHeatmap, openxlsx
# ================================================================

rm(list = ls())

# ---------------------------
# Load Required Packages
# ---------------------------
library(phyloseq)
library(tidyr)
library(tidyselect)
library(tidyverse)
library(ggplot2)
library(readxl)
library(ComplexHeatmap)
library(pheatmap)
library(microbiome)
library(qdapRegex)
library(magrittr)
library(biomformat)
library(openxlsx)

# ---------------------------
# Define Data Folder
# ---------------------------
dfol <- "HCCAnalysis/data/MDAnderson/"

# ---------------------------
# Load Preprocessed Microbiome Data
# ---------------------------
load(file.path(dfol, "processed_data_clref.rds"))

# ---------------------------
# Load and Process Metabolomics Data
# ---------------------------
# Prepare file for analysis
# metabolites_raw <- read.delim(file.path(dfol,"metabolites.tsv"),skip = 1)
# metabolites_raw %<>% rename(sampleid=1)

metabolites_raw <- read.csv(file.path(dfol, "2914e76aee36462c8fcfe11a2b244917.csv"))

# Extract column annotations (first 10 columns typically contain metadata)
ColAnot <- metabolites_raw %>%
  select_at(1:10)

# Extract metabolite intensity values
metabolites_raw %<>% select_at(-(1:10)) %>% t()

# Clean sample IDs
rownames(metabolites_raw) %<>% gsub("_R.*area", "", .)

# ---------------------------------------------------------------
# (Preserve commented original code for future reference)
# ---------------------------------------------------------------
# temp <- read.delim(file.path(dfol,"HFD/metabolite_feature_metadata.txt"))
# temp2 <- temp
# all(temp$sampleid %in% metabolites_raw$sampleid)
# temp <- temp %>%
#   select(sampleid, Compound_Name) %>%
#   left_join(metabolites_raw) %>%
#   filter(nchar(Compound_Name) >= 0) %>%
#   select(-sampleid)
# rownames(temp) <- NULL
# all(temp2$Compound_Name == temp$Compound_Name)
# ColAnot <- temp %>%
#   select('Compound_Name')
# length(unique(ColAnot$Compound_Name))
# metabolites_raw <- temp %>% select(-Compound_Name) %>% t() %>%
#   data.frame()
# ---------------------------------------------------------------

# ---------------------------
# Align Metabolomics with Microbiome Data
# ---------------------------

# Identify common samples
r_sel <- intersect(var2$host_subject_id, rownames(metabolites_raw))
# r_sel <- setdiff(var2$host_subject_id, rownames(metabolites_raw))

# Align rownames for all datasets
rownames(instruments)   <- var2$host_subject_id
rownames(phenotypes)    <- var2$host_subject_id
rownames(control_var)   <- var2$host_subject_id
sample_names(microbiome_raw) <- var2$host_subject_id

# Subset datasets by matched samples
metabolites_raw <- metabolites_raw[r_sel, ]
instruments     <- instruments[r_sel, ]
phenotypes      <- phenotypes[r_sel, ]
control_var     <- control_var[r_sel, ]
microbiome_raw_meta <- prune_samples(r_sel, microbiome_raw)

# Consistency checks across matched data
stopifnot(
  all(rownames(metabolites_raw) == rownames(instruments)),
  all(rownames(metabolites_raw) == rownames(phenotypes)),
  all(rownames(metabolites_raw) == rownames(control_var)),
  all(rownames(metabolites_raw) == sample_names(microbiome_raw_meta))
)

# ---------------------------
# Filter and Normalize Metabolites
# ---------------------------

# Filter low-variance metabolites
sd_metabolites <- (apply(metabolites_raw, 2, sd) %>% log()) > 0
ColAnot <- ColAnot[sd_metabolites, , drop = FALSE]
metabolites_raw <- metabolites_raw[, sd_metabolites]

# Centered log-ratio (CLR) normalization
metabolites_raw_clr <- apply(metabolites_raw, 1, function(x) log(x + (x == 0)) - mean(log(x + (x == 0))))
metabolites_raw_clr <- t(metabolites_raw_clr)

# ---------------------------
# Process Microbiome Data (Subset + CLR)
# ---------------------------

microbiome_use <- microbiome_raw_meta  # could aggregate at Family level if needed
# microbiome_use <- aggregate_taxa(microbiome_raw_meta, 'Family')

# Filter taxa based on abundance variability
x_filt <- subset_taxa(microbiome_use, apply(otu_table(microbiome_use), 1, sd) > 3)

# Rarefaction-like normalization (rescale counts)
mlib <- min(sample_sums(microbiome_raw_meta))
x_filt <- transform_sample_counts(x_filt, function(x) round(mlib * x / sum(x)))

# Extract taxonomy annotation
tax_anot <- tax_table(x_filt) %>% data.frame()

# Extract OTU table and apply CLR
x_filt <- t(x_filt@otu_table@.Data + 0)
Y <- x_filt  # raw abundance matrix

microbiome_use_clr <- apply(x_filt, 1, function(x) log(x + (x == 0)) - mean(log(x + (x == 0))))
microbiome_use_clr <- t(microbiome_use_clr)

# meta_metabolites <- temp2[sd_metabolites,]

# ---------------------------
# Save Preprocessed Data
# ---------------------------

save(instruments, phenotypes, ColAnot, metabolites_raw, microbiome_raw_meta,
     control_var, microbiome_use_clr, metabolites_raw_clr, # meta_metabolites,
     tax_anot,
     file = file.path(dfol, "processed_data_metabolites_otu.rds"))

message("✅ Saved combined microbiome-metabolite data to processed_data_metabolites_otu.rds")

# ================================================================
# Export Normalized Data to Excel (for reproducibility / inspection)
# ================================================================

fmicrobiome <- file.path(dfol , "microbiome_metabolites_data.xlsx")
if (file.exists(fmicrobiome)) file.remove(fmicrobiome)

wb <- createWorkbook()

# Add worksheets for each data component
addWorksheet(wb, "microbiome_clr")
addWorksheet(wb, "metabolites_normalized")
addWorksheet(wb, "Instruments")
addWorksheet(wb, "Confounder")
addWorksheet(wb, "Phenotypes")
addWorksheet(wb, "Metabolites_meta")
addWorksheet(wb, "microbiome_meta_taxa_table")

# Write each dataset to workbook
writeData(wb, sheet = "microbiome_clr", microbiome_use_clr %>% as.matrix() %>% scale(), rowNames = TRUE)
writeData(wb, sheet = "metabolites_normalized", metabolites_raw_clr %>% as.matrix() %>% scale(), rowNames = TRUE)
writeData(wb, sheet = "Instruments", instruments, rowNames = TRUE)
writeData(wb, sheet = "Confounder", control_var, rowNames = TRUE)
writeData(wb, sheet = "Phenotypes", phenotypes, rowNames = TRUE)
writeData(wb, sheet = "Metabolites_meta", ColAnot, rowNames = TRUE)
writeData(wb, sheet = "microbiome_meta_taxa_table", tax_anot, rowNames = TRUE)

# Final save
saveWorkbook(wb, fmicrobiome, overwrite = TRUE)

message("✅ Exported normalized datasets to: ", fmicrobiome)
