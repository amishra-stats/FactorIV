# ==========================================================
# Prepare phyloseq object for the engraftment data analysis
# Source: Diabimmune dataset (150 pick close referenced OTU)
# ==========================================================

# ----------------------------------------------------------
# Load required packages
# ----------------------------------------------------------
rm(list = ls())
library(phyloseq)
library(tidyr)
library(tidyselect)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(readxl)
library(ComplexHeatmap)
library(pheatmap)
library(microbiome)
library(qdapRegex)
library(magrittr)
library(biomformat)

# ----------------------------------------------------------
# Define data folder
# ----------------------------------------------------------
dfol <- 'Diabimmune/data/data_diabimmune/'

# ----------------------------------------------------------
# Read metadata files
# ----------------------------------------------------------

# Sample-level metadata
sample_data <- read.delim(file.path(dfol, "DiabImmune_Sample.txt")) %>%
  .[, apply(., 2, function(x) length(unique(x)) > 1)]
dim(sample_data)

# Participant-level metadata
participant_data <- read.delim(file.path(dfol, "DiabImmune_Participant.txt")) %>%
  .[, apply(., 2, function(x) length(unique(x)) > 1)]

# Clean column names for participant data
colnames(participant_data) %<>%
  gsub("\\.EUPATH.*", "", .) %>%     # Remove .EUPATH and everything after
  gsub("\\.\\.\\..*", "", .) %>%     # Remove ...EUPATH and everything after
  gsub("\\.", "_", .)                # Replace . with _

# ----------------------------------------------------------
# Select autoimmune-related outcomes
# ----------------------------------------------------------
autoimmune_outcomes <- c(
  "Age_at_two_autoantibodies_positive__days__",
  "Age_at_type_1_diabetes_diagnosis__days__",
  "Type_1_diabetes_diagnosed_",
  "HLA_risk__by_HLA_haplotyping_",
  "Glutamic_acid_decarboxylase_antibodies_",
  "Insulin_autoantibodies_",
  "Insulinoma_associated_protein_2_autoantibodies_",
  "Islet_cell_autoantibodies_",
  "Zinc_transporter_8_autoantibodies_"
)

# Extract autoimmune outcome table
outcome_table <- participant_data %>%
  select(Participant_Id, all_of(autoimmune_outcomes)) %>%
  column_to_rownames('Participant_Id')

# Keep participant metadata excluding outcome variables
participant_data <- participant_data %>%
  select(-all_of(autoimmune_outcomes)) %>%
  mutate(Participant_Id = as.character(Participant_Id))

# ----------------------------------------------------------
# Handle missing data in participant metadata
# ----------------------------------------------------------
missing_cols <- sapply(participant_data, function(x) any(is.na(x)))
sum(missing_cols)
colnames(participant_data[, missing_cols])

# Replace missing numeric values with mean and categorical with mode
participant_data %<>%
  mutate(across(where(is.numeric), ~ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
  mutate(across(where(is.character),
                ~ifelse(is.na(.), names(sort(table(.), decreasing = TRUE)[1]), .))) %>%
  column_to_rownames('Participant_Id')

# Sanity check for matching participant IDs
all(rownames(participant_data) == rownames(outcome_table))

# ----------------------------------------------------------
# Process participant repeated measure data
# ----------------------------------------------------------
participant_repeated_measure <- read.delim(file.path(dfol, "DiabImmune_Participant_Repeated_Measure.txt")) %>%
  .[, apply(., 2, function(x) length(unique(x)) > 1)]
colnames(participant_repeated_measure)[3] <- "Age_months"

# ----------------------------------------------------------
# Process ontology metadata
# ----------------------------------------------------------
ontology_metadata <- read.delim(file.path(dfol, "DiabImmune_OntologyMetadata.txt")) %>%
  .[, apply(., 2, function(x) length(unique(x)) > 1)]

# Participant-level metadata description
participants_data_description <- ontology_metadata %>%
  filter(category == "Participant")

# Species annotation for metagenomic and 16S assays
species_annotation_mtg <- ontology_metadata %>%
  filter(category == "Metagenomic sequencing assay") %>%
  filter(parentlabel == "Species")

species_annotation_16s <- ontology_metadata %>%
  filter(category == "16S rRNA (V4) assay") %>%
  filter(parentlabel == "Species")

# ----------------------------------------------------------
# Process metagenomic sequencing assay data
# ----------------------------------------------------------
metagenomic_sequencing_assay <- read.delim(file.path(dfol, "DiabImmune_Metagenomic_sequencing_assay.txt")) %>%
  .[, apply(., 2, function(x) length(unique(x)) > 1)]

colAnnotation <- colnames(metagenomic_sequencing_assay) %>%
  str_extract("EUPATH.*") %>%
  gsub("\\.$", "", .)

colnames(metagenomic_sequencing_assay)[-(1:4)] <- colAnnotation[-(1:4)]
select_col <- colAnnotation %in% species_annotation_mtg$iri
select_col[1:4] <- TRUE
metagenomic_sequencing_assay %<>% .[, select_col]

# ----------------------------------------------------------
# Process 16S rRNA assay data
# ----------------------------------------------------------
microbiome16s_assay <- read.delim(file.path(dfol, "DiabImmune_16S_rRNA_(V4)_assay.txt")) %>%
  .[, apply(., 2, function(x) length(unique(x)) > 1)]

colAnnotation <- colnames(microbiome16s_assay) %>%
  str_extract("EUPATH.*") %>%
  gsub("\\.$", "", .)

colnames(microbiome16s_assay)[-(1:4)] <- colAnnotation[-(1:4)]
select_col <- colAnnotation %in% species_annotation_16s$iri
select_col[1:4] <- TRUE
microbiome16s_assay %<>% .[, select_col]

# ----------------------------------------------------------
# Match repeated measures with microbiome data
# ----------------------------------------------------------

# Create empty lists for storing matched data
matched_data_mtg <- vector("list", length = 3)
matched_data_16s <- vector("list", length = 3)
sel_sample_ <- vector("list", length = 3)

# ----- Metagenomic data -----
temp_mtg <- participant_repeated_measure %>%
  left_join(metagenomic_sequencing_assay, by = c("Participant_Id", "participant_repeated_measure_Id")) %>%
  select(-Sample_Id, -Metagenomic_sequencing_assay_Id) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), 0, .))) %>%
  group_by(Participant_Id, participant_repeated_measure_Id, Age_months) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), 0, .)))

temp_mtg <- temp_mtg[rowSums(temp_mtg[, -(1:3)]) > 0.8, ]
temp_mtg <- cbind(temp_mtg[, (1:3)], temp_mtg[, -(1:3)][, colSums(temp_mtg[, -(1:3)]) > 0])

# ----- 16S rRNA data -----
temp_16s <- participant_repeated_measure %>%
  left_join(microbiome16s_assay, by = c("Participant_Id", "participant_repeated_measure_Id")) %>%
  select(-Sample_Id, -X16S_rRNA_.V4._assay_Id) %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), 0, .))) %>%
  group_by(Participant_Id, participant_repeated_measure_Id, Age_months) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ ifelse(is.nan(.), 0, .)))

temp_16s <- temp_16s[rowSums(temp_16s[, -(1:3)]) > 2000, ]
temp_16s <- cbind(temp_16s[, (1:3)], temp_16s[, -(1:3)][, colSums(temp_16s[, -(1:3)]) > 0])

# ----------------------------------------------------------
# Define age intervals (0–12, 12–24, >24 months)
# ----------------------------------------------------------
sel_sample_[[1]] <- list(
  mtg = temp_mtg %>% filter(Age_months <= 12),
  mb16s = temp_16s %>% filter(Age_months <= 12)
)
sel_sample_[[2]] <- list(
  mtg = temp_mtg %>% filter(Age_months <= 24) %>% filter(Age_months > 12),
  mb16s = temp_16s %>% filter(Age_months <= 24) %>% filter(Age_months > 12)
)
sel_sample_[[3]] <- list(
  mtg = temp_mtg %>% filter(Age_months > 24),
  mb16s = temp_16s %>% filter(Age_months > 24)
)

# ----------------------------------------------------------
# Select the last available sample per participant
# ----------------------------------------------------------
for (i in 1:length(matched_data_mtg)) {
  print(paste("Processing slice number:", i))
  sel_sample <- sel_sample_[[i]]
  
  # --- Metagenomic data ---
  sel_sample$mtg %<>%
    arrange(desc(-Age_months)) %>%
    group_by(Participant_Id) %>%
    dplyr::slice(1) %>%
    ungroup()
  Age_month_mtg <- sel_sample$mtg$Age_months
  sel_sample$mtg %<>%
    select(-Age_months, -participant_repeated_measure_Id) %>%
    column_to_rownames('Participant_Id')
  sel_sample$mtg %<>% .[rowSums(.) > 0.8, ] %>% .[, colSums(.) > 0]
  print(dim(sel_sample$mtg))
  
  # --- 16S data ---
  sel_sample$mb16s %<>%
    arrange(desc(-Age_months)) %>%
    group_by(Participant_Id) %>%
    dplyr::slice(1) %>%
    ungroup()
  Age_month_16s <- sel_sample$mb16s$Age_months
  sel_sample$mb16s %<>%
    select(-Age_months, -participant_repeated_measure_Id) %>%
    column_to_rownames('Participant_Id')
  sel_sample$mb16s %<>% .[rowSums(.) > 2000, ] %>% .[, colSums(.) > 0]
  print(dim(sel_sample$mb16s))
  
  # --- Match microbiome data with metadata ---
  matched_data_mtg[[i]] <- list(
    microbiome = sel_sample$mtg,
    participant_data = cbind(participant_data[rownames(sel_sample$mtg), ],
                             Age_months = Age_month_mtg),
    outcome_table = outcome_table[rownames(sel_sample$mtg), ],
    annotation = species_annotation_mtg %>%
      select(iri, label) %>%
      column_to_rownames('iri') %>%
      .[colnames(sel_sample$mtg), , drop = FALSE]
  )
  
  matched_data_16s[[i]] <- list(
    microbiome = sel_sample$mb16s,
    participant_data = cbind(participant_data[rownames(sel_sample$mb16s), ],
                             Age_months = Age_month_16s),
    outcome_table = outcome_table[rownames(sel_sample$mb16s), ],
    annotation = species_annotation_16s %>%
      select(iri, label) %>%
      column_to_rownames('iri') %>%
      .[colnames(sel_sample$mb16s), , drop = FALSE]
  )
}

# [1] "Processing slice number: 1"
# [1] 220 401
# [1]  256 1018
# [1] "Processing slice number: 2"
# [1] 224 417
# [1]  255 1056
# [1] "Processing slice number: 3"
# [1] 114 320
# [1] 165 948
# ----------------------------------------------------------
# Save processed data objects
# ----------------------------------------------------------
save(matched_data_16s, matched_data_mtg,
     file = file.path(dfol, 'processed_data_diabiimune.rds'))