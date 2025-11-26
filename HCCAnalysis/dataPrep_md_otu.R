# ================================================================
# File: dataPrep_md_otu.R
# Project: FactorIV - HCCAnalysis
# Purpose: Prepare microbiome (OTU) data and experimental metadata
#          for FactorIV causal inference analysis in HCC mouse study.
#
# This script:
#   1. Imports OTU tables and metadata.
#   2. Cleans, filters, and harmonizes sample IDs.
#   3. Constructs instrument variables (treatment, antibiotics, etc.).
#   4. Defines control variables (age, weight, strain).
#   5. Builds phenotype matrix from descriptive metadata.
#   6. Saves preprocessed data as `processed_data_clref.rds`.
#
# Dependencies: phyloseq, microbiome, magrittr, tidyverse, ComplexHeatmap
# ================================================================

rm(list = ls())

# ---------------------------
# Load Required Libraries
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

# ---------------------------
# Define Folder Paths
# ---------------------------
# The script assumes it is run from the HCCAnalysis folder.
# Adjust if your folder layout differs.

base_dir <- getwd()
data_dir <- file.path(base_dir, "HCCAnalysis","data", "MDAnderson")

# ---------------------------
# Load and Filter Microbiome Data
# ---------------------------

# Import OTU table as a phyloseq object
# microbiome_raw = import_biom(file.path(data_dir, "46116_otu_table.biom"))
microbiome_raw <- import_biom(file.path(data_dir, "25970_otu_table.biom"))

# Keep samples with sufficient sequencing depth (>10,000 reads)
microbiome_raw <- prune_samples(sample_sums(microbiome_raw) > 10000, microbiome_raw)

# ---------------------------
# Load Metadata and Match Samples
# ---------------------------

# temp <- read.delim(file.path(data_dir, "10856_prep_2458_20201201-081028.txt"))
temp <- read.delim(file.path(data_dir, "sample_information_from_prep_2458.tsv"))

# Remove rows with invalid age
temp <- temp %>% filter(age != "Not applicable")

# Drop columns with no variation
temp <- temp[, apply(temp, 2, function(x) length(unique(x))) > 1]

# Prune microbiome samples to those present in metadata
microbiome_raw <- prune_samples(temp$sample_id, microbiome_raw)

# Match sample order between metadata and phyloseq object
stopifnot(all(sample_names(microbiome_raw) %in% temp$sample_id))
temp <- temp[match(sample_names(microbiome_raw), temp$sample_id), ]
stopifnot(all(sample_names(microbiome_raw) == temp$sample_id))

# Rename taxonomy table columns for readability
colnames(tax_table(microbiome_raw)) <- c(
  "Kingdom", "Phylum", "Class", "Order",
  "Family", "Genus", "Species"
)

# ---------------------------
# Construct Instrument Variables
# ---------------------------

# Derive treatment information from "mice_group"
dict <- temp$mice_group %>%
  rm_between("(", ")") %>%
  gsub(" Col compare to SI", "", .) %>%
  str_split("-")

uniq_wd <- dict %>%
  unlist() %>%
  unique() %>%
  .[-grep("Group", .)]

# Remove uninformative terms
uniq_wd <- setdiff(uniq_wd, c("NC", "HFD", "Antibiotic", "In", "L1"))

# Create binary matrix for treatments
temp_out <- matrix(0, nrow(temp), length(uniq_wd)) %>% data.frame()
colnames(temp_out) <- uniq_wd

for (i in 1:nrow(temp)) {
  temp_out[i, uniq_wd %in% dict[[i]]] <- 1
}

# Merge overlapping or equivalent treatment labels
temp_out$STAM <- temp_out$STAM + temp_out$` STAM`
temp_out$` STAM` <- NULL
temp_out$STAM[temp_out$STAM == 2] <- 1

temp_out$Oxali <- temp_out$Oxali + temp_out$`Oxali+aPD`
temp_out$Oxali[temp_out$Oxali == 2] <- 1

temp_out$aPD <- temp_out$aPD + temp_out$`Oxali+aPD`
temp_out$aPD[temp_out$aPD == 2] <- 1
temp_out$`Oxali+aPD` <- NULL

# Summarize instrument columns
colSums(temp_out)

# ---------------------------
# Construct Additional Covariates
# ---------------------------

# Select relevant variables for modeling
var2 <- temp %>%
  select(age, diet, weight, antibiotics, strain, host_subject_id) %>%
  mutate(antibiotics = ifelse(antibiotics == "No", "No", "Yes"))

# Create design matrix for diet and antibiotics (instrument variables)
instruments <- select(var2, diet, antibiotics) %>%
  model.matrix(~ diet + antibiotics - 1, data = .) %>%
  cbind(temp_out, .)

colnames(instruments)

# ---------------------------
# Define Control Variables
# ---------------------------

# Convert weight from ranges (e.g., "20-25") to numeric mean
var2$weight <- var2$weight %>%
  str_split("-") %>%
  lapply(function(x) mean(as.numeric(x))) %>%
  unlist()

# Scale numeric controls and encode strain
control_var <- select(var2, age, weight) %>%
  apply(2, as.numeric) %>%
  scale() %>%
  cbind(model.matrix(~ strain, var2)[, 2]) %>%
  data.frame() %>%
  dplyr::rename("Strain" = 3)

# ---------------------------
# Construct Phenotype Matrix
# ---------------------------

dict <- temp$phenotype %>%
  # rm_between("(", ")") %>%   # intentionally commented, retain as-is
  gsub("HCC,CD8KO", "HCC, CD8KO\"", .) %>%
  gsub("NASH,lean", "NASH, lean", .) %>%
  gsub("CD8KO\"", "CD8KO", .) %>%
  gsub("HCC,no IgA", "HCC, no IgA", .) %>%
  gsub("FVB lean", "FVB, lean", .) %>%
  gsub("obese only", "obese", .) %>%
  gsub("CD8ko", "CD8KO", .) %>%
  gsub("NASH driven HCC", "NASH, HCC", .) %>%
  gsub("no b and t cells", "no B cell, no T cell", .) %>%
  gsub("obese no IgA", "obese, no IgA", .) %>%
  gsub("PigRK", "PIGRK", .) %>%
  gsub("FVB ", "FVB", .) %>%
  gsub("mon B ", "mon B cells", .) %>%
  gsub("obesity", "obese", .) %>%
  gsub("mon B cellscells", "mon B cells", .) %>%
  str_split(", ")

uniq_wd <- dict %>%
  unlist() %>%
  unique()

# Exclude general strain or tissue labels
uniq_wd <- setdiff(uniq_wd, c("CD8KO", "ABX", "PIGRK", "Colon",
                              "TGFbrko", "TGFbRflox/flox", "B6", "FVB", "SI"))

# Create phenotype matrix
temp_out <- matrix(0, nrow(temp), length(uniq_wd)) %>% data.frame()
colnames(temp_out) <- uniq_wd

for (i in 1:nrow(temp)) {
  temp_out[i, uniq_wd %in% dict[[i]]] <- 1
}

phenotypes <- temp_out

# Check correlations between control vars and phenotypes
cor(control_var, phenotypes)

# ---------------------------
# Save Preprocessed Data
# ---------------------------

save(instruments, phenotypes, control_var,
     microbiome_raw, var2,
     file = file.path(data_dir, "processed_data_clref.rds"))

message("✅ Data preparation complete: processed_data_clref.rds saved in ", data_dir)