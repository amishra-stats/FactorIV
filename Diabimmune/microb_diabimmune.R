# ==========================================================
# Script: microb_diabimmune.R
# Project: Diabimmune (Factor IV Microbiome Causal Analysis)
# Author: Aditya Mishra
# Purpose: Perform GOFAR-based multi-task modeling and
#          causal inference using microbiome and metagenomic data.
# ==========================================================

rm(list = ls())

# ----------------------------------------------------------
# HPC usage examples (kept for documentation)
# ----------------------------------------------------------
# module load GCC/13.3.0
# ml R/4.4.2-gfbf-2024a
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' -o out/out_3 -e err/err_3 "Rscript model_fit_gofar_joint_otu_wsd1.R Family 123"

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

# ----------------------------------------------------------
# Set folder paths (ensure output structure exists)
# ----------------------------------------------------------
# dfol <- 'data/data_diabimmune/'
dfol <- 'Diabimmune/data/data_diabimmune/'
outdir <- 'Diabimmune/output'
plotdir <- 'Diabimmune/plots' 
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

# ----------------------------------------------------------
# Load preprocessed data
# ----------------------------------------------------------
load(file.path(dfol, 'processed_data_diabiimune.rds'))

# ----------------------------------------------------------
# Analysis configuration
# ----------------------------------------------------------
index_select <- 1
mtg_file <- sprintf('%s/microb_diabimmune_mtg_gofar%d.RData', outdir, index_select)
rna_file <- sprintf('%s/microb_diabimmune_16s_gofar%d.RData', outdir, index_select)

mb16s <- matched_data_16s[[index_select]]
mtgs <- matched_data_mtg[[index_select]]

ctrl_var <- c('Country__ENVO_00000009_', 'Height_for_age_z_score_',
              'Sex__PATO_0000047_', 'BMI_for_age_z_score_')

# ==========================================================
# 1. Metagenomic Data Processing
# ==========================================================

# --- Prepare clr-transformed data for modeling
Y_mtg <- mtgs$microbiome %>%
  .[, apply(., 2, function(x) sum(x > 0) > 0.1 * nrow(.))] %>%
  sweep(1, rowSums(., na.rm = TRUE), FUN = "/") %>%
  add(1e-6) %>% t() %>% clr() %>% t() %>%
  as.matrix() %>% scale()

# --- Prepare instruments and controls
instruments_mtg <- mtgs$participant_data %>%
  dplyr::select(-one_of(ctrl_var)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  model.matrix(~ ., data = .) %>% .[, -1]

control_var_mtg <- mtgs$participant_data %>%
  dplyr::select(one_of(ctrl_var)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  model.matrix(~ ., data = .) %>% .[, -1]

phenotypes_mtg <- mtgs$outcome_table[, -(1:2)]
X_mtg <- cbind(control_var_mtg, instruments_mtg) %>% as.matrix()

# --- GOFAR configuration
family <- list(gaussian(), binomial(), poisson())
SD <- 123; nlam <- 30
if (index_select == 1) {
  control <- gofar_control(epsilon = 1e-10, spU = 0.99, spV = 0.98,
                           maxit = 5000, initmaxit = 10000,
                           initepsilon = 1e-12, gamma0 = 0.1,
                           equalphi = 0, lamMaxFac = 1,
                           lamMinFac = 1e-6, se = 0)
} else {
  control <- gofar_control(epsilon = 1e-10, spU = 0.99, spV = 0.98,
                           maxit = 5000, initmaxit = 10000,
                           initepsilon = 1e-20, gamma0 = 0.1,
                           equalphi = 0, lamMaxFac = 1,
                           lamMinFac = 1e-6, se = 0)
}

familygroup <- rep(1, ncol(Y_mtg))
set.seed(SD)
rank.est <- 5

# --- Load prefit model (or run GOFAR if needed)
# fit_seq_mtg <- gofar_p(Y_mtg, X_mtg, cIndex = 1:ncol(control_var_mtg),
#                    nrank = rank.est, family = family,
#                    nlambda = nlam, familygroup = familygroup,
#                    control = control, nfold = 5)
# fit_seq_mtg$fit[[1]]$dev
# save(fit_seq_mtg,file = mtg_file)
load(mtg_file)
# plot(fit_seq_mtg$cv_details[[1]]$mean_dev)
# plot(fit_seq_mtg$cv_details[[2]]$mean_dev)

# ==========================================================
# 2. 16S rRNA Data Processing
# ==========================================================
Y_16s <- mb16s$microbiome %>% round() %>%
  .[, apply(., 2, function(x) sum(x > 0) > 0.1 * nrow(.))] %>%
  sweep(1, rowSums(., na.rm = TRUE), FUN = "/") %>%
  add(1e-6) %>% t() %>% clr() %>% t() %>%
  as.matrix() %>% scale()

instruments_16s <- mb16s$participant_data %>%
  dplyr::select(-one_of(ctrl_var)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  model.matrix(~ ., data = .) %>% .[, -1]

control_var_16s <- mb16s$participant_data %>%
  dplyr::select(one_of(ctrl_var)) %>%
  mutate(across(where(is.numeric), scale)) %>%
  model.matrix(~ ., data = .) %>% .[, -1]

phenotypes_16s <- mb16s$outcome_table[, -(1:2)]
X_16s <- cbind(control_var_16s, instruments_16s) %>% as.matrix()

family <- list(gaussian(), binomial(), poisson())
SD <- 123; nlam <- 30

# --- GOFAR control setup (same logic retained)
if (index_select == 1) {
  control <- gofar_control(epsilon = 1e-10, spU = 0.99, spV = 0.98,
                           maxit = 5000, initmaxit = 10000, 
                           initepsilon = 1e-12, gamma0 = 0.1, equalphi = 0,
                           lamMaxFac = 1, lamMinFac = 1e-6, se1 = 0)
  familygroup <- rep(1,ncol(Y_16s)); set.seed(SD); rank.est <- 5
  # fit_seq_16s <- gofar_p(Y_16s, X_16s, cIndex = 1:ncol(control_var_16s),
  #                        nrank = rank.est, family = family,
  #                        nlambda = nlam, familygroup = familygroup,
  #                        control = control, nfold = 5)
} else if (index_select == 2) {
  control <- gofar_control(epsilon = 1e-10, spU = 0.99, spV = 0.98,
                           maxit = 5000, initmaxit = 10000, 
                           initepsilon = 1e-12, gamma0 = 0.1, equalphi = 0,
                           lamMaxFac = 1, lamMinFac = 1e-6, se1 = 0)
  familygroup <- rep(1,ncol(Y_16s)); set.seed(SD); rank.est <- 5
  # fit_seq_16s <- gofar_p(Y_16s, X_16s, cIndex = 1:ncol(control_var_16s),
  #                        nrank = rank.est, family = family,
  #                        nlambda = nlam, familygroup = familygroup,
  #                        control = control, nfold = 5)
} else {
  control <- gofar_control(epsilon = 1e-10, spU = 0.99, spV = 0.98,
                           maxit = 5000, initmaxit = 10000, 
                           initepsilon = 1e-20, gamma0 = 0.1, equalphi = 0,
                           lamMaxFac = 1, lamMinFac = 1e-6, se1 = 0)
  familygroup <- rep(1,ncol(Y_16s)); set.seed(SD); rank.est <- 5
  # fit_seq_16s <- gofar_rrr(Y_16s, X_16s, cIndex = 1:ncol(control_var_16s),
  #                        maxrank = rank.est, family = family,
  #                        familygroup = familygroup,
  #                        control = control, nfold = 5)
}
# fit_seq_16s$U
# save(fit_seq_16s,file = rna_file)
load(rna_file)
# plot(colMeans(fit_seq_16s$fit[[1]]$dev))
# plot(colMeans(fit_seq_16s$fit[[2]]$dev))
# plot(colMeans(fit_seq_16s$fit[[3]]$dev))


# ==========================================================
# 3. IV-Probit Model for Causal Discovery (16S)
# ==========================================================
sink(file.path(outdir, sprintf('causal_effect_16s_gofar%d.txt', index_select)))

Z <- X_16s[, -(1:ncol(control_var_16s))]  # instruments
Q <- X_16s[, 1:ncol(control_var_16s)]     # controls

phenotypes <- phenotypes_16s %>%
  mutate(
    HLA_risk_binary = case_when(
      HLA_risk__by_HLA_haplotyping_ == 2 ~ "Low",
      HLA_risk__by_HLA_haplotyping_ %in% c(3, 4) ~ "High",
      TRUE ~ NA_character_
    ),
    HLA_risk_binary = factor(HLA_risk_binary, levels = c("Low", "High"))
  ) %>%
  dplyr::select(-HLA_risk__by_HLA_haplotyping_) %>%
  mutate(across(where(is.character), ~ factor(., levels = c("No", "Yes")))) %>%
  data.frame()

r_est <- length(fit_seq_16s$D)
U_hat <- Y_16s %*% fit_seq_16s$V
IV_hat <- as.matrix(Z) %*% fit_seq_16s$U %*% diag(fit_seq_16s$D, r_est, r_est)

Q <- data.frame(Q); colnames(Q) <- paste0('Q', seq_len(ncol(Q)))
C <- data.frame(Q); colnames(C) <- paste0('C', seq_len(ncol(Q)))
IV_hat <- data.frame(IV_hat); colnames(IV_hat) <- paste0('IV', seq_len(r_est))
U_hat <- data.frame(U_hat); colnames(U_hat) <- paste0('U', seq_len(r_est))

causal_effect <- matrix(0, ncol(Y_16s), ncol(phenotypes))
model_summary <- NULL

for (i in 1:ncol(phenotypes)) {
  na_index <- !is.na(phenotypes[, i])
  df <- cbind(y = as.numeric(phenotypes[, i]) - 1, U_hat, Q, IV_hat, C)
  mdlX <- as.formula(paste(
    'y ~', paste(colnames(U_hat), collapse = '+'), '|',
    paste(colnames(Q), collapse = '+'), '|',
    paste(c(colnames(IV_hat), colnames(C)), collapse = '+')
  ))
  fit_gmm <- ivprobit(mdlX, df[na_index, ])
  wt <- fit_gmm$coefficients * (fit_gmm$pval < 0.1)
  print(colnames(phenotypes)[i])
  summary(fit_gmm)
  cat('\n--------------------------------------\n\n')
  causal_effect[, i] <- fit_seq_16s$V %*% wt[1 + (1:r_est)]
  model_summary %<>%
    rbind(data.frame(summary(fit_gmm)) %>%
            rownames_to_column('Variable') %>%
            mutate(Outcome = colnames(phenotypes)[i]))
}

colnames(causal_effect) <- colnames(phenotypes)
rownames(causal_effect) <- colnames(Y_16s)
causal_effect_16s <- as.data.frame(causal_effect)
model_summary_16s <- model_summary
sink()

save(causal_effect_16s, model_summary_16s,
     file = file.path(outdir, sprintf('causal_effect_16s_gofar%d.RData', index_select)))

# ==========================================================
# 4. IV-Probit Model for Causal Discovery (Metagenomics)
# ==========================================================
sink(file.path(outdir, sprintf('causal_effect_mtg_gofar%d.txt', index_select)))

Z <- X_mtg[, -(1:ncol(control_var_mtg))]
Q <- X_mtg[, 1:ncol(control_var_mtg)]

phenotypes <- phenotypes_mtg %>%
  mutate(
    HLA_risk_binary = case_when(
      HLA_risk__by_HLA_haplotyping_ == 2 ~ "Low",
      HLA_risk__by_HLA_haplotyping_ %in% c(3, 4) ~ "High",
      TRUE ~ NA_character_
    ),
    HLA_risk_binary = factor(HLA_risk_binary, levels = c("Low", "High"))
  ) %>%
  dplyr::select(-HLA_risk__by_HLA_haplotyping_) %>%
  mutate(across(where(is.character), ~ factor(., levels = c("No", "Yes")))) %>%
  data.frame()

r_est <- length(fit_seq_mtg$D)
U_hat <- Y_mtg %*% fit_seq_mtg$V
IV_hat <- as.matrix(Z) %*% fit_seq_mtg$U %*% diag(fit_seq_mtg$D, r_est, r_est)

Q <- data.frame(Q); colnames(Q) <- paste0('Q', seq_len(ncol(Q)))
C <- data.frame(Q); colnames(C) <- paste0('C', seq_len(ncol(Q)))
IV_hat <- data.frame(IV_hat); colnames(IV_hat) <- paste0('IV', seq_len(r_est))
U_hat <- data.frame(U_hat); colnames(U_hat) <- paste0('U', seq_len(r_est))

causal_effect <- matrix(0, ncol(Y_mtg), ncol(phenotypes))
model_summary <- NULL

for (i in 1:ncol(phenotypes)) {
  na_index <- !is.na(phenotypes[, i])
  df <- cbind(y = as.numeric(phenotypes[, i]) - 1, U_hat, Q, IV_hat, C)
  mdlX <- as.formula(paste(
    'y ~', paste(colnames(U_hat), collapse = '+'), '|',
    paste(colnames(Q), collapse = '+'), '|',
    paste(c(colnames(IV_hat), colnames(C)), collapse = '+')
  ))
  fit_gmm <- ivprobit(mdlX, df[na_index, ])
  wt <- fit_gmm$coefficients * (fit_gmm$pval < 0.1)
  print(colnames(phenotypes)[i])
  summary(fit_gmm)
  cat('\n--------------------------------------\n\n')
  causal_effect[, i] <- fit_seq_mtg$V %*% wt[1 + (1:r_est)]
  model_summary %<>%
    rbind(data.frame(summary(fit_gmm)) %>%
            rownames_to_column('Variable') %>%
            mutate(Outcome = colnames(phenotypes)[i]))
}

colnames(causal_effect) <- colnames(phenotypes)
rownames(causal_effect) <- colnames(Y_mtg)
causal_effect_mtg <- as.data.frame(causal_effect)
model_summary_mtg <- model_summary
sink()

save(causal_effect_mtg, model_summary_mtg,
     file = file.path(outdir, sprintf('causal_effect_mtg_gofar%d.RData', index_select)))


