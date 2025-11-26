# helper for NULL names
`%||%` <- function(a, b) if (!is.null(a)) a else b

estimate_alpha_naive <- function(X, Y, Q = NULL,
                                 p_adj_method = "BH",
                                 alpha_level = 0.05) {
  # Ensure proper structures
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  if (!is.null(Q)) {
    Q <- as.data.frame(Q)
  }
  
  n <- nrow(X)
  d <- ncol(X)
  
  if (length(Y) != n) {
    stop("Y must have the same number of rows as X.")
  }
  if (!is.null(Q) && nrow(Q) != n) {
    stop("Q must have the same number of rows as X.")
  }
  
  # Storage for results
  est   <- numeric(d)
  se    <- numeric(d)
  pval  <- numeric(d)
  
  for (j in seq_len(d)) {
    df <- data.frame(
      y = Y,
      x = X[, j]
    )
    if (!is.null(Q)) {
      df <- cbind(df, Q)
    }
    
    fit <- lm(y ~ ., data = df)
    sm  <- summary(fit)
    
    # Coefficient for x (first predictor after intercept)
    # name is exactly "x" because we called the column x
    coef_x <- sm$coefficients["x", ]
    
    est[j]  <- coef_x["Estimate"]
    se[j]   <- coef_x["Std. Error"]
    pval[j] <- coef_x["Pr(>|t|)"]
  }
  
  # Multiple testing correction across all X_j
  p_adj <- p.adjust(pval, method = p_adj_method)
  
  # Zero out non-significant effects
  est_shrunk <- ifelse(p_adj < alpha_level, est, 0)
  
  # Build result table
  res <- data.frame(
    variable    = colnames(X) %||% paste0("X", seq_len(d)),
    alpha_hat   = est,
    se          = se,
    p_value     = pval,
    p_adj       = p_adj,
    alpha_hat_shrunk = est_shrunk,
    significant = p_adj < alpha_level
  )
  
  res
}




estimate_alpha_naive_deseq2 <- function(X, Q, y, alpha_level = 0.05) {
  # X: n x d count matrix (samples x features)
  # Q: n x c matrix/data.frame of observed confounders
  # y: length-n vector (numeric or factor) – outcome / exposure of interest
  # alpha_level: FDR threshold for significance (on padj)
  #
  # Returns:
  #  - alpha_naive: length-d vector of estimated effects (0 if not significant)
  #  - pval: length-d vector of raw p-values
  #  - padj: length-d vector of adjusted p-values
  #  - res_df: full DESeq2 result table (with feature IDs)
  
  if (!is.matrix(X)) X <- as.matrix(X)
  
  n <- nrow(X)
  d <- ncol(X)
  
  # -----------------------------
  # 1. Build DESeq2 dataset
  # -----------------------------
  # Confounders Q as data.frame with clean names
  Q_df <- as.data.frame(Q)
  
  # Outcome / exposure of interest
  col_df <- cbind(Q_df, y = y)
  rownames(col_df) <- paste0("sample_", seq_len(n))
  
  # Count matrix must be features x samples
  countData <- t(round(X))
  rownames(countData) <- paste0("feat_", seq_len(d))
  colnames(countData) <- rownames(col_df)
  
  # Build design: ~ Q1 + Q2 + ... + y
  q_terms <- colnames(Q_df)
  design_formula <- as.formula(
    paste("~", paste(c(q_terms, "y"), collapse = " + "))
  )
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = countData,
    colData   = col_df,
    design    = design_formula
  )
  
  # -----------------------------
  # 2. Fit NB GLM via DESeq2
  # -----------------------------
  dds <- DESeq2::DESeq(
    dds,
    test    = "Wald",
    fitType = "local",
    sfType  = "poscounts"
  )
  
  # Identify the coefficient corresponding to 'y'
  coef_names <- DESeq2::resultsNames(dds)
  y_coef_name <- grep("^y", coef_names, value = TRUE)[1]
  if (is.na(y_coef_name)) {
    stop("Could not find a coefficient corresponding to 'y' in resultsNames(dds).")
  }
  
  res <- DESeq2::results(dds, name = y_coef_name, cooksCutoff = FALSE)
  res_df <- as.data.frame(res)
  res_df$feature_id <- rownames(res_df)
  
  # -----------------------------
  # 3. Extract effects and significance
  # -----------------------------
  # Use log2FoldChange as the effect of y on each feature
  alpha_hat <- res_df$log2FoldChange
  pval <- res_df$pvalue
  padj <- res_df$padj
  
  # Replace NA padj with 1 (definitely not significant)
  padj_clean <- ifelse(is.na(padj), 1, padj)
  
  # Set non-significant effects to 0
  alpha_naive <- alpha_hat
  alpha_naive[padj_clean > alpha_level] <- 0
  
  # -----------------------------
  # 4. Return structured output
  # -----------------------------
  list(
    alpha_naive = alpha_naive,
    pval        = pval,
    padj        = padj,
    res_df      = res_df
  )
}


