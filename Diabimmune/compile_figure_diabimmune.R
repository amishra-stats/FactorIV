# ==========================================================
# Script: compile_figure_diabimmune.R
# Project: Diabimmune (Factor IV Microbiome Causal Analysis)
# Purpose: Compile, visualize, and summarize causal effect results
# Author: Aditya Mishra
# ==========================================================

rm(list = ls())

# ----------------------------------------------------------
# Load dependencies
# ----------------------------------------------------------
library(phyloseq)
library(gofar)
library(microbiome)
library(magrittr)
library(DESeq2)
library(compositions)
library(tidyverse)
library(ivprobit)
library(dplyr)
library(stringr)
library(kableExtra)

# ----------------------------------------------------------
# Set folder paths
# ----------------------------------------------------------
dfol <- 'Diabimmune/data/data_diabimmune/'
outdir <- 'Diabimmune/output'
plotdir <- 'Diabimmune/plots' 
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

# ----------------------------------------------------------
# Load processed data
# ----------------------------------------------------------
load(file.path(dfol, 'processed_data_diabiimune.rds'))

# Initialize
index_select <- 1
causal_effect_16s_ <- data.frame()
causal_effect_mtg_ <- data.frame()

# ==========================================================
# 1. Aggregate Causal Effect Results Across Timepoints
# ==========================================================
for (index_select in 1:3) {
  mtg_file <- sprintf('%s/microb_diabimmune_mtg_gofar%d.RData', outdir, index_select)
  rna_file <- sprintf('%s/microb_diabimmune_16s_gofar%d.RData', outdir, index_select)
  ce_rna_file <- file.path(dfol, sprintf('causal_effect_16s_gofar%d.RData', index_select))
  ce_mtg_file <- file.path(dfol, sprintf('causal_effect_mtg_gofar%d.RData', index_select))
  
  load(mtg_file)
  load(rna_file)
  load(ce_mtg_file)
  load(ce_rna_file)
  
  mtgs <- matched_data_mtg[[index_select]]
  ctrl_var <- c('Country__ENVO_00000009_', 'Height_for_age_z_score_',
                'Sex__PATO_0000047_', 'BMI_for_age_z_score_')
  
  # --- Append 16S causal effects
  causal_effect_16s_ %<>% rbind(
    causal_effect_16s %>%
      as.data.frame() %>%
      rownames_to_column('microbiome_id') %>%
      cbind(Annotation = matched_data_16s[[index_select]]$annotation[.$microbiome_id, ], .) %>%
      mutate(index_select = index_select)
  )
  
  # --- Append metagenomic causal effects
  causal_effect_mtg_ %<>% rbind(
    causal_effect_mtg %>%
      as.data.frame() %>%
      rownames_to_column('microbiome_id') %>%
      cbind(Annotation = matched_data_mtg[[index_select]]$annotation[.$microbiome_id, ], .) %>%
      mutate(index_select = index_select)
  )
}

# ----------------------------------------------------------
# Filter causal effects for plotting
# ----------------------------------------------------------
temp <- colSums(abs(causal_effect_mtg_[, -(1:2)])) != 0
phenotypes_select <- c(colnames(causal_effect_mtg_)[1:2], names(temp)[temp])
plot_causaleffect <- causal_effect_mtg_[, phenotypes_select]

# ==========================================================
# 2. Metagenomic Causal Effects Heatmap
# ==========================================================

# Step 1: Long-format transformation
long_df <- plot_causaleffect[, -2] %>%
  pivot_longer(
    cols = -c(Annotation, index_select),
    names_to = "Phenotype",
    values_to = "CausalEffect"
  )

# Step 2: Filter and normalize
filtered_df <- long_df %>%
  group_by(index_select, Phenotype) %>%
  filter(any(CausalEffect != 0)) %>%
  ungroup()

normalized_df <- filtered_df %>%
  group_by(index_select, Phenotype) %>%
  mutate(NormEffect = CausalEffect / max(abs(CausalEffect))) %>%
  ungroup() %>%
  mutate(NormEffect = ifelse(is.na(NormEffect), 0, NormEffect))

normalized_df <- normalized_df %>%
  mutate(Annotation = fct_reorder(Annotation, NormEffect, .fun = max))

# Step 3: Select top species per timepoint
top_annotations <- normalized_df %>%
  group_by(index_select, Annotation) %>%
  summarise(TopAbs = max(abs(NormEffect)), .groups = "drop") %>%
  group_by(index_select) %>%
  slice_max(TopAbs, n = 25) %>%
  ungroup() %>%
  distinct(Annotation)

filtered_top_df <- normalized_df %>%
  filter(Annotation %in% top_annotations$Annotation)

# Step 4: Label timepoints and simplify phenotype names
filtered_top_df <- filtered_top_df %>%
  mutate(
    index_select = recode_factor(
      as.factor(index_select),
      `1` = "0–12\nMonth",
      `2` = "12–24 Month",
      `3` = ">24 Month"
    ),
    Phenotype = str_replace_all(Phenotype, c(
      "Glutamic_acid_decarboxylase_antibodies_" = "GADA",
      "Islet_cell_autoantibodies_" = "ICA",
      "Zinc_transporter_8_autoantibodies_" = "ZnT8A",
      "HLA_risk_binary" = "HLA-risk",
      "Insulin_autoantibodies_" = "IAA",
      "Insulinoma_associated_protein_2_autoantibodies_" = "IA-2A"
    ))
  )

# Step 5: Plot causal heatmap
ggplot(filtered_top_df, aes(x = Phenotype, y = fct_reorder(Annotation, NormEffect, max), fill = NormEffect)) +
  geom_tile(color = "white") +
  facet_grid(~ index_select, labeller = label_value, scales = 'free', space = 'free') +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
  labs(
    title = "Top Metagenomic Species with Causal Effects on Autoimmune Phenotypes",
    x = "Phenotype", y = "Species (Annotation)", fill = "Normalized Effect"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "top",
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0)
  )

ggsave(
  filename = file.path(plotdir, "causal_heatmap_metagenomics_diabimmune.pdf"),
  plot = last_plot(), width = 8, height = 10, dpi = 500, units = "in"
)

# ==========================================================
# 3. Instrument Contribution and Outcome Association Plots
# ==========================================================

# (original logic preserved, only re-indented and with output paths fixed)
name_map <- c(
  "Age_1st_animal_milk_or_solids_given__days__" = "Age solids (days)",
  "Age_1st_given_apple__months__" = "Age apple",
  "Age_1st_given_banana__months__" = "Age banana",
  "Age_1st_given_barley__months__" = "Age barley",
  "Age_1st_given_beef__months__" = "Age beef",
  "Age_1st_given_berries__months__" = "Age berries",
  "Age_1st_given_cabbs__months__" = "Age cabbage",
  "Age_1st_given_carrot__months__" = "Age carrot",
  "Age_1st_given_corn__months__" = "Age corn",
  "Age_1st_given_cowsmilk__months__" = "Age cow milk",
  "Age_1st_given_egg__months__" = "Age egg",
  "Age_1st_given_fish__months__" = "Age fish",
  "Age_1st_given_milk_products__months__" = "Age milk prod",
  "Age_1st_given_oat__months__" = "Age oat",
  "Age_1st_given_pear__months__" = "Age pear",
  "Age_1st_given_peas__months__" = "Age peas",
  "Age_1st_given_plum__months__" = "Age plum",
  "Age_1st_given_pork__months__" = "Age pork",
  "Age_1st_given_potato__months__" = "Age potato",
  "Age_1st_given_poultry__months__" = "Age poultry",
  "Age_1st_given_rice__months__" = "Age rice",
  "Age_1st_given_rye__months__" = "Age rye",
  "Age_1st_given_sweet_potato__months__" = "Age sweet potato",
  "Age_1st_given_tomato__months__" = "Age tomato",
  "Age_1st_given_wheat__months__" = "Age wheat",
  "Age_at_anthropometry__days__" = "Anthropometry age",
  "Antibiotics_before_delivery__by_maternal_report_No" = "AB before delivery (No)",
  "Antibiotics_before_delivery__by_maternal_report_Yes" = "AB before delivery (Yes)",
  "Breastfed_duration_" = "Breastfed duration",
  "Delivery_mode_Vaginal" = "Vaginal delivery",
  "Diet_in_first_three_days_Mother's breast milk" = "Diet 1st 3d: Breast milk",
  "Diet_in_first_three_days_Multiple types or not reported" = "Diet 1st 3d: Mixed",
  "Exclusive_breastfed_duration_" = "Exclusive BF duration",
  "Gestational_diabetes__by_maternal_report_Yes" = "Gest. diabetes (Yes)",
  "Linear_growth_during_1st_year__cm__" = "Growth 1st year",
  "Maternal_age_at_birth__year__" = "Maternal age",
  "Mean_linear_growth_during_1st_3_years__cm_year__" = "Growth 1–3 yrs",
  "Mean_weight_gain_during_1st_3_years__kg_year__" = "Wt gain 1–3 yrs",
  "Study_group_Three country cohort (Karelia)" = "Cohort: Karelia",
  "Study_group_Type I Diabetes (T1D) cohort" = "Cohort: T1D",
  "Urban_or_rural_site_Rural" = "Site: Rural",
  "Urban_or_rural_site_Urban" = "Site: Urban",
  "Weight_gain_during_1st_year__kg__" = "Wt gain 1st year",
  "Weight_for_age_z_score_" = "Wt-for-age Z",
  "Age_months" = "Age (months)",
  'Country__ENVO_00000009_Finland' = "Country Finland",
  'Country__ENVO_00000009_Russia' = "Country Russia",
  'Height_for_age_z_score_' = "Height",
  'Sex__PATO_0000047_Male' = "Sex",
  'BMI_for_age_z_score_' = "BMI"
)

# ------------------------------------------------------------------
# Relative contribution of instruments across timepoints
# (metagenomics + 16S; formatting-only cleanup, logic unchanged)
# ------------------------------------------------------------------

var_place <- c(
  "Type_1_diabetes_diagnosed_"                    = "T1D",
  "Glutamic_acid_decarboxylase_antibodies_"       = "GADA",
  "Islet_cell_autoantibodies_"                    = "ICA",
  "Zinc_transporter_8_autoantibodies_"            = "ZnT8A",
  "HLA_risk_binary"                               = "HLA-risk",
  "Insulin_autoantibodies_"                       = "IAA",
  "Insulinoma_associated_protein_2_autoantibodies_" = "IA-2A"
)

plot_instruments <- data.frame()
plot_outcome     <- data.frame()

for (index_select in 1:3) {
  # model fit + causal-effect files
  mtg_file   <- sprintf("%s/microb_diabimmune_mtg_gofar%d.RData", outdir, index_select)
  # rna_file <- sprintf("%s/microb_diabimmune_16s_gofar%d.RData", outdir, index_select)
  
  load(mtg_file)
  # load(rna_file)
  
  ce_mtg_file <- file.path(dfol, sprintf("causal_effect_mtg_gofar%d.RData", index_select))
  load(ce_mtg_file)
  
  # Latent factors with significant outcomes (labels kept compact)
  sel_latent <- model_summary_mtg %>%
    filter(grepl("U", Variable)) %>%                                   # optionally: filter(p.val < 0.1)
    dplyr::select(Variable, Outcome, p.val) %>%
    mutate(Variable = gsub("U", "V", Variable)) %>%
    mutate(Outcome  = str_replace_all(Outcome, var_place)) %>%
    mutate(Outcome  = ifelse(p.val < 0.1, paste0(Outcome, " "), "")) %>%
    group_by(Variable) %>%
    summarise(Outcome = paste(Outcome, collapse = ""), .groups = "drop") %>%
    dplyr::rename(LatentFactor = 1)
  
  # Prepare microbiome scores (metagenomics)
  mtgs <- matched_data_mtg[[index_select]]
  ctrl_var <- c("Country__ENVO_00000009_", "Height_for_age_z_score_",
                "Sex__PATO_0000047_", "BMI_for_age_z_score_")
  
  Y_mtg <- mtgs$microbiome %>%
    .[, apply(., 2, function(x) sum(x > 0) > 0.1 * nrow(.))] %>%
    sweep(1, rowSums(., na.rm = TRUE), FUN = "/") %>%
    add(1e-6) %>% t() %>% clr() %>% t() %>%
    as.matrix() %>% scale()
  
  microbiome_score <- (Y_mtg %*% fit_seq_mtg$V) %>%
    data.frame() %>%
    cbind(mtgs$outcome_table[, -(1:2)]) %>%
    mutate(
      HLA_risk_binary = case_when(
        HLA_risk__by_HLA_haplotyping_ == 2       ~ "Low",
        HLA_risk__by_HLA_haplotyping_ %in% c(3,4) ~ "High",
        TRUE                                      ~ NA_character_
      ),
      HLA_risk_binary = factor(HLA_risk_binary, levels = c("Low", "High"))
    ) %>%
    dplyr::select(-HLA_risk__by_HLA_haplotyping_) %>%
    mutate(index_select = index_select) %>%
    dplyr::rename(V1 = 1, V2 = 2) %>%
    gather("LatentFactor", "Score", -(Type_1_diabetes_diagnosed_:index_select))
  
  # Shorten phenotype column names
  matched_names <- intersect(names(var_place), colnames(microbiome_score))
  colnames(microbiome_score)[match(matched_names, colnames(microbiome_score))] <- var_place[matched_names]
  
  # Collect outcome-annotated scores for the selected latent factors
  temp <- sel_latent %>% filter(nchar(Outcome) > 0)
  
  for (k in 1:nrow(temp)) {
    lf <- temp$LatentFactor[k]
    colselect <- str_trim(temp$Outcome[k]) %>% strsplit(" ") %>% unlist() %>%
      c("index_select", "LatentFactor", "Score")
    
    plot_outcome %<>%
      rbind(
        microbiome_score %>%
          filter(LatentFactor %in% lf) %>%
          .[, colselect] %>%
          gather("Type", "Outcome", -index_select, -LatentFactor, -Score)
      )
  }
  
  # Instruments and controls (renaming for plotting only)
  instruments_mtg <- mtgs$participant_data %>%
    dplyr::select(-one_of(ctrl_var)) %>%
    mutate(across(where(is.numeric), scale)) %>%
    model.matrix(~ ., data = .) %>% .[, -1]
  
  control_var_mtg <- mtgs$participant_data %>%
    dplyr::select(one_of(ctrl_var)) %>%
    mutate(across(where(is.numeric), scale)) %>%
    model.matrix(~ ., data = .) %>% .[, -1]
  
  matched_names <- intersect(names(name_map), colnames(instruments_mtg))
  colnames(instruments_mtg)[match(matched_names, colnames(instruments_mtg))] <- name_map[matched_names]
  
  matched_names <- intersect(names(name_map), colnames(control_var_mtg))
  colnames(control_var_mtg)[match(matched_names, colnames(control_var_mtg))] <- name_map[matched_names]
  
  plot_instruments %<>%
    rbind(
      fit_seq_mtg$U %>%
        as.data.frame() %>%
        mutate(Instruments = colnames(instruments_mtg),
               index_select = index_select) %>%
        gather("LatentFactor", "Contribution", -Instruments, -index_select) %>%
        filter(LatentFactor %in% sel_latent$LatentFactor) %>%
        left_join(sel_latent, by = "LatentFactor")
    )
}

# Quick sanity outputs
plot_instruments$LatentFactor %>% unique()
head(plot_instruments)
head(plot_outcome)

# Label timepoints (for plots)
plot_df <- plot_outcome %>%
  mutate(
    Timepoint = recode_factor(
      as.factor(index_select),
      `1` = "0–12 Month",
      `2` = "12–24 Month",
      `3` = ">24 Month"
    )
  ) %>%
  filter(nchar(Outcome) > 0)

# ---------------- Boxplot of LF scores by outcome ----------------
ggplot(plot_df, aes(x = Type, y = Score, fill = Outcome)) +
  geom_boxplot(outlier.size = 0.8, width = 0.7, alpha = 0.8, notch = TRUE) +
  facet_grid(Timepoint + LatentFactor ~ ., scales = "free", space = "free") +
  ggpubr::stat_compare_means(
    aes(group = Outcome),
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE,
    label.y.npc = "top"
  ) +
  labs(
    title = "Latent Microbiome Factor Scores by Outcome",
    x = "Phenotype Type", y = "Latent Factor Score", fill = "Outcome"
  ) +
  coord_flip() +
  theme_minimal(base_size = 13) +
  theme(
    strip.text         = element_text(face = "bold"),
    axis.text.x        = element_text(angle = 0, hjust = 0.5),
    plot.title.position= "plot",
    plot.title         = element_text(hjust = 0),
    legend.position    = "top"
  )

ggsave(
  filename = file.path(plotdir, "causal_boxplot_metagenomics_diabimmune.pdf"),
  plot     = last_plot(),
  width    = 5, height = 10, dpi = 500, units = "in"
)

# ---------------- Bar plot: top instrument contributions ----------------
plot_df <- plot_instruments %>%
  mutate(
    AbsContribution = abs(Contribution),
    Timepoint = recode_factor(
      as.factor(index_select),
      `1` = "0–12 Month",
      `2` = "12–24 Month",
      `3` = ">24 Month"
    )
  )

top_df <- plot_df %>%
  group_by(Timepoint, LatentFactor) %>%
  slice_max(AbsContribution, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Instruments = fct_reorder(Instruments, AbsContribution))

ggplot(top_df, aes(x = AbsContribution, y = Instruments, fill = Contribution > 0)) +
  geom_col(show.legend = FALSE) +
  facet_grid(Timepoint ~ LatentFactor, scales = "free", space = "free") +
  scale_fill_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2")) +
  labs(
    title = "Top Instrument Contributions to Latent Microbiome Factors",
    x = "Absolute Contribution", y = "Instrument"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text         = element_text(face = "bold"),
    plot.title.position= "plot",
    plot.title         = element_text(hjust = 0)
  )

ggsave(
  filename = file.path(plotdir, "causal_barplot_metagenomics_diabimmune.pdf"),
  plot     = last_plot(),
  width    = 6, height = 10, dpi = 500, units = "in"
)


