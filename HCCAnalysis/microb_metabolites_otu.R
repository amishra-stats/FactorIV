# ==========================================================
# File: microb_metabolites_otu.R  (HCCAnalysis)
# Purpose: Joint metabolite + microbiome modeling with GOFAR,
#          followed by downstream causal analysis & visualization
# Note: This section includes setup, data loading, and GOFAR config.
# ==========================================================

rm(list = ls())

# ----------------------------------------------------------
# HPC job examples (kept verbatim for reference)
# ----------------------------------------------------------
# Latest
# bsub -W 24:00 -q medium -n 1 -M 32 -R 'rusage[mem=32]' -o out/out_2 -e err/err_2 "Rscript model_fit_gofar_joint_otu.R Family 123 1>Rout/out_2 2>Rerr/err_2"
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' -o out/out_2 -e err/err_2 "Rscript model_fit_gofar_joint_otu.R Family 123 1>Rout/out_2 2>Rerr/err_2"
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' -o out/out_1 -e err/err_1 "Rscript model_fit_gofar_mbolites.R Family 123 1>Rout/out_1 2>Rerr/err_1"
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' -o out/out_3 -e err/err_3 "Rscript model_fit_gofar_mbiom.R Family 123 1>Rout/out_3 2>Rerr/err_3"

## Latest run 
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' -o out/out_2 -e err/err_2 "Rscript model_fit_gofar_joint_otu_wsd.R Family 123 1>Rout/out_2 2>Rerr/err_2"
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' -o out/out_3 -e err/err_3 "Rscript model_fit_gofar_joint_otu_wsd1.R Family 123 1>Rout/out_3 2>Rerr/err_3"

# ----------------------------------------------------------
# Libraries
# ----------------------------------------------------------
library(phyloseq)
library(gofar)
library(nbfar)
library(microbiome)
library(magrittr)
library(DESeq2)

# ----------------------------------------------------------
# Paths & data loading
# ----------------------------------------------------------
dfol <- 'HCCAnalysis/data/MDAnderson/'             # data folder (from HCCAnalysis root)
if (!dir.exists('output')) dir.create('output', recursive = TRUE)  # ensure output/ exists

load(file.path(dfol, 'processed_data_metabolites_otu.rds'))  # provides: instruments, control_var, phenotypes, metabolites_raw_clr, microbiome_use_clr, etc.


# ----------------------------------------------------------
# Construct design matrices for GOFAR
# ----------------------------------------------------------
Z <- instruments         # instruments (e.g., diet, antibiotics + constructed indicators)
Q <- control_var         # controls    (e.g., age, weight, strain)

# Output file (two lines as in original; second assignment will be used)
# metabolites_analysis1.rds gamma0 <- 2
# metabolites_analysis.rds  gamma0 <- 1
# metabolites_analysis_full.rds gamma0 <- 1
# ofile <- 'output/microb_metabolites_full_otu.RData'
ofile <- file.path('HCCAnalysis', 'output/microb_metabolites_full_otu_wsd1.RData')

# Response (Y): joint metabolite + microbiome CLR matrices
Y  <- cbind(metabolites_raw_clr, microbiome_use_clr) %>% as.matrix()
Y2 <- scale(Y, scale = TRUE, center = TRUE)   # standardized copy (exploration)
# hist(apply(Y2, 2, sd),100)
# range(apply(metabolites_raw_clr, 2, sd))
# range(apply(microbiome_use_clr, 2, sd))

# Predictors (X): controls + instruments
X <- cbind(Q, Z) %>% as.matrix()

# ----------------------------------------------------------
# GOFAR configuration (kept identical to original)
# ----------------------------------------------------------
family  <- list(gaussian(), binomial(), poisson())
control <- gofar_control()
nlam    <- 20     # number of tuning parameters
SD      <- 123
q       <- ncol(Y)
p       <- ncol(X)

# Fine-tuning control parameters
control$epsilon   <- 1e-7
control$spU       <- 0.99
control$spV       <- 0.98
control$maxit     <- 1000
control$gamma0    <- 1
control$equalphi  <- 0
control$lamMaxFac <- 1
control$lamMinFac <- 1e-10

familygroup <- rep(1, q)

# ----------------------------------------------------------
# Model fitting: GOFAR(S) on full data
# (loading a prefit model, but original fit code retained as comments)
# ----------------------------------------------------------
set.seed(SD)
rank.est <- 5

# fit.seq <- gofar_s(Y, X, cIndex = 1:3,
#                    nrank = rank.est, family = family,
#                    nlambda = nlam, familygroup = familygroup,
#                    control = control, nfold = 5
# )
# Sys.getenv()
# save(list=ls(), file = ofile)

load(ofile)   # load precomputed fit.seq (or equivalent object) from `ofile`
# Diagnostic of the cross validation error 
fit.seq$fit[[1]]$dev
plot(colMeans(fit.seq$fit[[1]]$dev))
plot(colMeans(fit.seq$fit[[2]]$dev))
plot(colMeans(fit.seq$fit[[3]]$dev))





# ==========================================================
# Section: Extract GOFAR loadings and construct causal IV model
# ==========================================================

# Extract estimated loadings from GOFAR model
Uest <- fit.seq$U
rownames(Uest) <- colnames(Z)

Vest <- fit.seq$V
rownames(Vest) <- c(ColAnot$row.ID, colnames(microbiome_use_clr))

Cest <- fit.seq$C
rownames(Cest) <- colnames(Z)
colnames(Cest) <- rownames(Vest)

# ----------------------------------------------------------
# Summarize and filter nonzero loadings
# ----------------------------------------------------------
library(tidyverse)
rname <- rowSums(Vest != 0) > 0

# ==========================================================
# Section: Instrumental Variable (IV) causal discovery
# ==========================================================

file.create(file.path('HCCAnalysis','output/microbe_metabolites_analysis_otu.txt'))
sink(file.path('HCCAnalysis','output/microbe_metabolites_analysis_otu.txt'))

library(ivprobit)

# Prepare IV components from GOFAR outputs
r_est <- length(fit.seq$D)
U_hat <- Y %*% fit.seq$V                                     # latent response scores
IV_hat <- as.matrix(Z) %*% fit.seq$U %*% diag(fit.seq$D, r_est, r_est)  # latent instruments

# Label matrices for clarity
Q <- data.frame(Q); colnames(Q) <- paste('Q', 1:ncol(Q), sep = '')
C <- data.frame(Q); colnames(C) <- paste('C', 1:ncol(C), sep = '')
IV_hat <- data.frame(IV_hat); colnames(IV_hat) <- paste('IV', 1:r_est, sep = '')
U_hat  <- data.frame(U_hat);  colnames(U_hat)  <- paste('U', 1:r_est, sep = '')

# (Original commented formulas preserved)
# mdlX <- as.formula(paste('y ~', paste(c(colnames(U_hat),colnames(Q)), collapse = '+')))
# mdlIV <- as.formula(paste('~', paste(c(colnames(IV_hat),colnames(C)), collapse = '+')))

# ----------------------------------------------------------
# Estimate causal effects for each phenotype
# ----------------------------------------------------------
causal_effect <- matrix(0, ncol(Y), ncol(phenotypes))
model_summary <- NULL

for (i in 1:ncol(phenotypes)) {
  
  df <- cbind(y = phenotypes[, i], U_hat, Q, IV_hat, C)
  
  # IV probit model formula
  mdlX <- as.formula(paste(
    'y ~', paste(colnames(U_hat), collapse = '+'), '|',
    paste(colnames(Q), collapse = '+'), '|',
    paste(c(colnames(IV_hat), colnames(C)), collapse = '+')
  ))
  
  # Fit IV Probit model
  fit_gmm <- ivprobit(mdlX, df)
  
  # Extract significant coefficient weights
  wt <- fit_gmm$coefficients * (fit_gmm$pval < 0.05)
  
  # Print model summary
  print(colnames(phenotypes)[i])
  summary(fit_gmm)
  print('--------------------------------------')
  cat('\n\n')
  
  # Compute causal effect: linear combination of V loadings and weights
  causal_effect[, i] <- fit.seq$V %*% wt[1 + (1:r_est)]
  
  # Append model summary table
  model_summary %<>%
    rbind(
      data.frame(summary(fit_gmm)) %>%
        rownames_to_column('Variable') %>%
        mutate(Outcome = colnames(phenotypes)[i])
    )
}

# Label columns and store
colnames(causal_effect) <- colnames(phenotypes)
causal_effect <- data.frame(causal_effect)
sink()

# ==========================================================
# Section: Compute microbiome–metabolite “balance” scores
# ==========================================================

# Project Y onto latent spaces separately for metabolites and microbiome
mtb_balance <- Y[, 1:ncol(metabolites_raw_clr)] %*% 
  fit.seq$V[1:ncol(metabolites_raw_clr), ]

mb_balance <- Y[, -(1:ncol(metabolites_raw_clr))] %*% 
  fit.seq$V[-(1:ncol(metabolites_raw_clr)), ]

# Combine into long-form dataset for visualization
plot_df <- data.frame(mtb_balance) %>%
  mutate(Type = 'Metabolite balance') %>%
  cbind(phenotypes) %>%
  rbind(
    data.frame(mb_balance) %>%
      mutate(Type = 'Microbiome balance') %>%
      cbind(phenotypes)
  ) %>%
  gather('Inflamation', 'Group', -(1:3)) %>%
  mutate(Group = ifelse(Group == 0, 'No', 'Yes'))

# ==========================================================
# Section: Statistical testing and boxplot visualization
# ==========================================================

library(ggplot2)
library(ggpubr)
library(dplyr)
library(rstatix)

# Perform pairwise t-tests for each facet
stat_test <- plot_df %>%
  group_by(Type, Inflamation) %>%
  t_test(X1 ~ Group) %>%  # test for first latent dimension
  adjust_pvalue(method = "bonferroni") %>%
  rstatix::add_xy_position(x = 'Group')

# Annotate significance levels
stat_test <- stat_test %>%
  mutate(p.adj.signif = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01  ~ "**",
    p.adj <= 0.05  ~ "*",
    TRUE           ~ "ns"
  ))

stat_test2 <- stat_test %>%
  select(group1, group2, p.adj.signif, y.position, Type, Inflamation)

# ==========================================================
# Section: Visualization - Boxplots for first two latent factors
# ==========================================================

# ---------- First latent balance (X1) ----------
p <- ggplot(plot_df, aes(x = Group, y = X1)) +
  geom_boxplot(aes(fill = Group), notch = TRUE) +
  facet_grid(Type ~ Inflamation, scales = 'free', space = 'fixed') +
  stat_pvalue_manual(stat_test, size = 3, label = "p.adj.signif") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "lightblue", color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Boxplot with Statistical Significance",
    y = "First Microbiome/Metabolite Balance",
    x = "Liver Phenotype"
  )

print(p)
pdf('HCCAnalysis/plots/balance_v1.pdf', width = 10, height = 5)
print(p)
dev.off()

# ---------- Second latent balance (X2) ----------
p <- ggplot(plot_df, aes(x = Group, y = X2)) +
  geom_boxplot(aes(fill = Group), notch = TRUE) +
  facet_grid(Type ~ Inflamation, scales = 'free', space = 'fixed') +
  stat_pvalue_manual(stat_test, size = 3, label = "p.adj.signif") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "lightblue", color = "black"),
    legend.position = "none"
  ) +
  labs(
    title = "Boxplot with Statistical Significance",
    y = "Second Microbiome/Metabolite Balance",
    x = "Liver Phenotype"
  )

print(p)
pdf('HCCAnalysis/plots/balance_v2.pdf', width = 10, height = 5)
print(p)
dev.off()







# ==========================================================
# Section: Visualization — Latent Factor Heatmaps
# ==========================================================

library(ComplexHeatmap)
library(circlize)
library(grid)

# ----------------------------------------------------------
# 1️⃣ Heatmap of estimated effect sizes (coefficients)
# ----------------------------------------------------------
small_mat <- model_summary %>%
  filter(Variable != 'Intercep') %>%                 # exclude intercept term
  select(Variable, Outcome, Coef) %>%
  spread(Outcome, Coef) %>%                          # reshape to wide format
  column_to_rownames('Variable') %>%
  t() %>%                                            # transpose for visualization
  round(2)

# Define color scale for effect sizes
col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("green", "white", "red"))

# Generate heatmap object
ht <- Heatmap(
  small_mat,
  name = "Effect Size",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)

# Draw to screen
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top")

# Save to PDF
pdf('HCCAnalysis/plots/lf_effect.pdf', width = 5, height = 6)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top")
dev.off()


# ----------------------------------------------------------
# 2️⃣ Heatmap of log p-values
# ----------------------------------------------------------
small_mat <- model_summary %>%
  filter(Variable != 'Intercep') %>%
  mutate(logPvalue = log(p.val)) %>%
  select(Variable, Outcome, logPvalue) %>%
  spread(Outcome, logPvalue) %>%
  column_to_rownames('Variable') %>%
  t() %>%
  round(2)

# Define color scale for log p-values
col_fun <- circlize::colorRamp2(c(-5, 0), c("blue", "white"))

# Generate heatmap
ht <- Heatmap(
  small_mat,
  name = "log\nP-value",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)

# Draw and save
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top")
pdf('HCCAnalysis/plots/lf_effect2.pdf', width = 5, height = 6)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top")
dev.off()


# ----------------------------------------------------------
# 3️⃣ Relative importance of intervention (U loadings)
# ----------------------------------------------------------
col_fun <- circlize::colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

# Select first two latent factors for visualization
small_mat <- Uest %>%
  data.frame() %>%
  rename('U1' = 1, 'U2' = 2) %>%
  as.matrix()

# Create heatmap for intervention importance
ht <- Heatmap(
  small_mat,
  name = "Effect",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", small_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)

# Save to file
pdf('HCCAnalysis/plots/rel_imp_intervention.pdf', width = 5, height = 6)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top")
dev.off()



# ----------------------------------------------------------
# ==========================================================
# Section: Integrating GNPS Metabolite Annotations
# ==========================================================

# ----------------------------------------------------------
# Load GNPS annotation file
# ----------------------------------------------------------
# fname <- file.path(
#   '/Users/akm93206/Documents/ProjectMD/IV-manuscript/IV/data/',
#   'ProteoSAFe-FEATURE-BASED-MOLECULAR-NETWORKING-c0a717dd-view_plotting_data/',
#   'DB_result/',
#   '8c98c7d7d6b04922b80b638de83dcdee.tsv'
# )
# file.copy(
#   fname,
#   'HCCAnalysis/data/MDAnderson/gnps_metabolite_annotations.tsv',
#   overwrite = TRUE
# )
fname <- 'HCCAnalysis/data/MDAnderson/gnps_metabolite_annotations.tsv'
mtb_gnps <- read.csv(fname, sep = '\t')
colnames(mtb_gnps) %<>% gsub('[.]', '', .)

# ----------------------------------------------------------
# Define metabolite type classification function
# ----------------------------------------------------------
mtbType <- function(x) {
  y <- rep('Other', length(x))
  x <- tolower(x)
  
  # Classify metabolites by text patterns
  y[grepl('cholic', x)]        <- "Bile acid"
  y[grepl('ampicillin', x)]    <- "Antibiotic"
  y[grepl('amoxacillin', x)]   <- "Antibiotic"
  
  aminos <- c('Trp-Ile', 'Trp-Phe', 'Ile-Gly-Ile', 'Trp-Ile', 'Trp-Asp', 
              'Trp-Val', 'Phe-Trp', 'Asp-Trp') %>% tolower()
  for (i in aminos) y[grepl(i, x)] <- "Amino acid"
  
  return(y)
}

# ----------------------------------------------------------
# Annotate metabolite table with GNPS compound types
# ----------------------------------------------------------
ColAnotN <- left_join(ColAnot, mtb_gnps, by = c('row.ID' = 'XScan')) %>%
  mutate(Compound_type = mtbType(Compound_Name))

# Display joint view of causal effects and annotations
cbind(data.frame(causal_effect[1:ncol(metabolites_raw_clr), ]), ColAnotN)

# Save annotated tables
openxlsx::write.xlsx(
  causal_effect[1:ncol(metabolites_raw_clr), ],
  file = 'HCCAnalysis/output/causal_effect_mtb.xlsx'
)
openxlsx::write.xlsx(
  ColAnotN,
  file = 'HCCAnalysis/output/metabolites_annotation.xlsx'
)


# ==========================================================
# Section: Select metabolites and taxa for visualization
# ==========================================================

# Select metabolite types of interest
sel_metb <- which(ColAnotN$Compound_type %in% 
                    c("Bile acid", "Antibiotic", "Amino acid"))

# Select microbiome taxa of interest
sel_mb_gnv <- which(tax_anot$Species == "s__gnavus")
sel_mb_ent <- which(tax_anot$Family == "f__Enterobacteriaceae")
sel_mb <- c(sel_mb_gnv, sel_mb_ent)
sel_mb

# Combined annotation labels
temp_ano <- c(
  ColAnotN$Compound_type[sel_metb],
  rep("R.gnavus", length(sel_mb_gnv)),
  rep("Enterobacteriaceae", length(sel_mb_ent))
)

# ----------------------------------------------------------
# Build matrix of causal effects for selected features
# ----------------------------------------------------------
plt_mat <- causal_effect[c(sel_metb, ncol(metabolites_raw_clr) + sel_mb), ] %>%
  .[, colSums(.) != 0] %>%  # remove zero columns
  t()

rownames(plt_mat) <- gsub('[.]', '', rownames(plt_mat))
plt_mat <- apply(plt_mat, 1, function(x) x / max(abs(x))) %>% t()  # normalize rowwise


# ==========================================================
# Section: Combined Causal Effect Heatmap
# ==========================================================

set.seed(123)

# Define annotation colors
type_anot <- data.frame(Type = gsub(' ', ' ', temp_ano))
col_temp  <- RColorBrewer::brewer.pal(8, 'Dark2')
group_col <- col_temp[1:5]
names(group_col) <- unique(type_anot$Type)

# Column annotations (metabolite/microbiome type)
col_anot <- columnAnnotation(
  df = type_anot,
  col = list(Type = group_col),
  annotation_name_gp = gpar(fontsize = 10),
  simple_anno_size = unit(0.3, "cm")
)

# ----------------------------------------------------------
# Define color mapping and clustering for heatmap
# ----------------------------------------------------------
library(circlize)
colfun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Hierarchical clustering on transposed causal effect matrix
a <- hclust(dist(t(plt_mat)), method = "ward.D2")$order

# Create heatmap
ht <- Heatmap(
  plt_mat,
  name = "Causal Effect",
  col = colfun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  column_title_rot = 30,
  column_split = type_anot$Type,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_gap = unit(2, "mm"),
  jitter = TRUE,
  border = TRUE
)

# Draw and save to PDF
pdf('HCCAnalysis/plots/causal_effect_mtb_mb.pdf', width = 8, height = 5)
draw(ht, heatmap_legend_side = "right", annotation_legend_side = "top")
dev.off()





# ==========================================================
# Section: Exporting Analysis Outputs to Excel Workbook
# ==========================================================

library(openxlsx)

# ----------------------------------------------------------
# Define output Excel file path
# ----------------------------------------------------------
fmicrobiome <- sprintf("HCCAnalysis/output/microb_metabolites_analysis.xlsx")

# Remove existing file if present (to avoid appending)
if (file.exists(fmicrobiome)) {
  file.remove(fmicrobiome)
}

# ----------------------------------------------------------
# Create and structure workbook
# ----------------------------------------------------------
wb <- createWorkbook()

# Add sheets corresponding to figures and model components
addWorksheet(wb, 'Figure V')                 # Full V loadings (metabolites + microbiome)
addWorksheet(wb, 'Figure V metabolite')      # V loadings restricted to metabolites
addWorksheet(wb, 'Figure V microbiome')      # V loadings restricted to microbiome
addWorksheet(wb, 'Figure U')                 # U loadings (instrument effects)
addWorksheet(wb, 'Figure C')                 # C matrix (covariate-adjusted effects)
addWorksheet(wb, 'Causal Effect')            # Causal effects for all phenotypes
# addWorksheet(wb, 'Figure V (summary V1)')
# addWorksheet(wb, 'Figure V (summary V2)')
# addWorksheet(wb, 'Figure V (summary V3)')
# addWorksheet(wb, 'Figure V (summary V4)')
# addWorksheet(wb, 'Figure V (summary V5)')
addWorksheet(wb, 'Metabolites_meta')         # GNPS annotation metadata
addWorksheet(wb, 'microbiome_meta')          # Microbiome taxonomic annotation

# ----------------------------------------------------------
# Write each data object to the corresponding sheet
# ----------------------------------------------------------
writeData(wb, sheet = 'Figure U', Uest, rowNames = TRUE)
writeData(wb, sheet = 'Figure V', Vest, rowNames = TRUE)

# Separate metabolite and microbiome components in V
writeData(wb, sheet = 'Figure V metabolite', Vest[1:nrow(ColAnot), ], rowNames = TRUE)
writeData(wb, sheet = 'Figure V microbiome', Vest[-(1:nrow(ColAnot)), ], rowNames = TRUE)

# Write additional model outputs
writeData(wb, sheet = 'Figure C', Cest, rowNames = TRUE)
writeData(wb, sheet = 'Causal Effect', causal_effect, rowNames = TRUE)

# (Optional summaries retained for reproducibility)
# writeData(wb, sheet = 'Figure V (summary V1)', dtf1, rowNames = TRUE)
# writeData(wb, sheet = 'Figure V (summary V2)', dtf2, rowNames = TRUE)
# writeData(wb, sheet = 'Figure V (summary V3)', dtf3, rowNames = TRUE)
# writeData(wb, sheet = 'Figure V (summary V4)', dtf2, rowNames = TRUE)
# writeData(wb, sheet = 'Figure V (summary V5)', dtf3, rowNames = TRUE)

# Metadata sheets
writeData(wb, sheet = 'Metabolites_meta', ColAnot, rowNames = TRUE)
writeData(wb, sheet = 'microbiome_meta', tax_anot, rowNames = TRUE)

# ----------------------------------------------------------
# Save workbook to disk
# ----------------------------------------------------------
saveWorkbook(wb, fmicrobiome, overwrite = TRUE)

# ==========================================================
# End of Export Section
# ==========================================================




# ==========================================================
# Section: Differential Abundance Analysis (Microbiome)
# ==========================================================

rm(list = ls())

# ----------------------------------------------------------
# Cluster run configurations (retained for reproducibility)
# ----------------------------------------------------------
# Example submission commands for HPC/LSF environments:
# bsub -W 24:00 -q medium -n 1 -M 32 -R 'rusage[mem=32]' \
#      -o out/out_2 -e err/err_2 "Rscript model_fit_gofar_joint_otu.R Family 123"
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' \
#      -o out/out_1 -e err/err_1 "Rscript model_fit_gofar_mbolites.R Family 123"
# bsub -W 24:00 -q medium -n 6 -M 64 -R 'rusage[mem=64]' \
#      -o out/out_3 -e err/err_3 "Rscript model_fit_gofar_mbiom.R Family 123"

# ==========================================================
# Load dependencies
# ==========================================================
library(phyloseq)
library(gofar)
library(nbfar)
library(microbiome)
library(magrittr)
library(tidyverse)
library(DESeq2)

# ----------------------------------------------------------
# Load preprocessed microbiome–metabolite data
# ----------------------------------------------------------
dfol <- 'HCCAnalysis/data/MDAnderson/'
load(file.path(dfol, 'processed_data_metabolites_otu.rds'))

# (Optional) argument parsing for batch runs
# args <- commandArgs(trailingOnly = TRUE)
# setting <- as.character(args[1])
# seed_id <- as.numeric(args[2])

# ==========================================================
# Step 1: Prepare data for DESeq2 analysis
# ==========================================================

# Clean phenotype variable names (remove spaces)
colnames(phenotypes) <- gsub(' ', '_', colnames(phenotypes))
selVars <- colnames(phenotypes)

# Prepare object containers
out_list_res <- vector('list', length = length(selVars))

# Use the microbiome data at OTU level
microbiome_use <- microbiome_raw_meta  # could aggregate later at 'Family' level if needed

# Filter out low-variance taxa
x_filt <- subset_taxa(microbiome_use, apply(otu_table(microbiome_use), 1, sd) > 3)
sum(taxa_sums(x_filt) == 0)

# Attach phenotype and control variable metadata
sample_data(x_filt) <- phenotypes %>% cbind(control_var, .)

# Sanity check: taxa ordering
all(colnames(microbiome_use_clr) == taxa_names(x_filt))

# ==========================================================
# Step 2: Differential abundance testing for each phenotype
# ==========================================================
for (i in seq_along(out_list_res)) {
  # Create DESeq2 object using phenotype i as outcome
  diagdds <- phyloseq_to_deseq2(
    x_filt,
    as.formula(paste('~ age + weight + Strain +', selVars[i]))
  )
  
  # Run DESeq2 with Wald test and local dispersion fit
  diagdds <- DESeq(diagdds, test = "Wald", fitType = "local", sfType = "poscounts")
  
  # Extract results and annotate
  out_list_res[[i]] <- results(diagdds, cooksCutoff = FALSE) %>%
    data.frame() %>%
    cbind(Phenotype = selVars[i], .) %>%
    cbind(tax_anot, .)
  
  # Remove missing adjusted p-values
  out_list_res[[i]] <- out_list_res[[i]][!is.na(out_list_res[[i]]$padj), ]
}

# Combine all phenotype-level results into one dataframe
out_deseq <- lapply(out_list_res, function(x) rownames_to_column(x, 'OTU_id'))
res_all <- do.call('rbind', out_deseq)

# ==========================================================
# Step 3: Save results
# ==========================================================
write.csv(res_all, file = 'HCCAnalysis/output/deseq_microbiome.csv', row.names = FALSE)



