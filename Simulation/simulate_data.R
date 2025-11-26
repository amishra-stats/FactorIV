simulate_X_given_U <- function(U, n, Xsigma) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Please install MASS for mvrnorm()")
  }
  
  basis_vec <- function(x) {
    if (diff(dim(x)) < 0) x <- t(x)
    qr_x <- qr(x)
    k <- qr.Q(qr_x) %*% qr.R(qr_x)[, 1:qr_x$rank]
    k[abs(k) < 1e-6] <- 0
    idx <- integer(qr_x$rank)
    for (i in seq_len(qr_x$rank)) {
      idx[i] <- which(apply(x, 2, function(col, ref) sum(abs(col - ref)), k[, i]) < 1e-6)[1]
    }
    list(ind = idx, vec = x[, idx, drop = FALSE])
  }
  
  p <- nrow(U)
  r <- ncol(U)
  # Complement basis
  U_t <- diag(max(dim(U)))
  U_t <- U_t[, -basis_vec(U)$ind, drop = FALSE]
  P <- cbind(U, U_t)
  
  UtXsUt <- t(U_t) %*% Xsigma %*% U_t
  UtXsU  <- t(U_t) %*% Xsigma %*% U
  UXsU   <- t(U)   %*% Xsigma %*% U
  UXsUinv <- solve(UXsU)
  sigma_X2 <- UtXsUt - UtXsU %*% UXsUinv %*% t(UtXsU)
  sigma_X2 <- (sigma_X2 + t(sigma_X2)) / 2  # symmetrize
  
  X1 <- matrix(rnorm(n * r), nrow = r, ncol = n)
  mean_X2 <- UtXsU %*% UXsUinv %*% X1
  X2_noise <- t(MASS::mvrnorm(ncol(mean_X2), rep(0, nrow(mean_X2)), sigma_X2))
  X2 <- mean_X2 + X2_noise
  
  X <- t(solve(t(P)) %*% rbind(X1, X2))
  X
}



simulate_iv_data <- function(
    n = 200, d = 100, m = 20, c = 3, k = 3,
    diagonal_elements = c(3, 4, 5),
    tem = 0.5,                    # correlation parameter for correlated error
    snr_X = 3, snr_Y = 3,         # target SNRs
    disp = 2,                     # NB dispersion
    family = "Gaussian",          # "Gaussian" or "NB"
    endog_mode = "correlated_error", # "confounder" or "correlated_error"
    seed = 123) {
  set.seed(seed)
  library(MASS)
  
  # -------------------------------------------
  # 1. Latent structure and parameters
  # -------------------------------------------
  D <- diag(diagonal_elements)
  
  # Z <- matrix(rnorm(n * m), n, m)
  
  # Sparse βU, βV
  betaU <- matrix(runif(m * k, -1, 1), m, k)
  betaU[abs(betaU) < 0.75] <- 0
  # betaU <- sweep(betaU, 2, sqrt(colSums(betaU^2)), FUN = "/")
  
  # Instruments
  # set.seed(123)
  # p <- 20; r <- 3; n <- 200
  # U <- matrix(rnorm(p * r), p, r); U <- qr.Q(qr(U))  # orthonormal basis
  Xsigma <- 0.5^abs(outer(1:m, 1:m, "-"))
  Z <- simulate_X_given_U(betaU, n, Xsigma)
  
  betaV <- matrix(runif(d * k, -1, 1), d, k)
  betaV[abs(betaV) < 0.75] <- 0
  betaV <- sweep(betaV, 2, sqrt(colSums(betaV^2)), FUN = "/")
  betaV <- t(betaV)
  
  beta <- betaU %*% D %*% betaV
  U <- Z %*% betaU
  
  # Confounders
  Q <- matrix(rnorm(n * c), n, c)
  gam <- matrix(rnorm(c * d, sd = 0.5), c, d)
  eta <- matrix(rnorm(c * 1, sd = 0.5), c, 1)
  
  # True causal effects
  # alpha <- matrix(rnorm(d * 1), d, 1)
  # alpha[abs(alpha) < 1] <- 0
  alpha <- matrix(runif(k * 1, 2, 4), k, 1)
  alpha <- crossprod(betaV, alpha)
  
  # -------------------------------------------
  # 2. Endogeneity mechanisms
  # -------------------------------------------
  if (endog_mode == "confounder") {
    # Shared latent confounder
    Uc <- matrix(rnorm(n * 1), n, 1)
    lambda_X <- sqrt(tem)
    lambda_Y <- sqrt(tem)
    
    # First-stage signal and noise
    temp <- lambda_X*sample(c(0,1),d,replace = T, prob = c(0.95,0.05))
    signal_X <- Z %*% beta + Q %*% gam +  Uc %*% t(temp)
    eps_raw <- matrix(rnorm(n * d), n, d)
    scale_to_snr <- function(signal, noise, snr_target) {
      var_signal <- mean(signal^2)
      var_noise  <- mean(noise^2)
      noise * sqrt(var_signal / (snr_target * var_noise))
    }
    eps <- scale_to_snr(signal_X, eps_raw, snr_X)
    
    
    # Generate X (Gaussian or NB)
    if (tolower(family) == "nb" || tolower(family) == "negative_binomial") {
      LP <- -0.5 + signal_X + 0*eps
      mu <- 10*exp(LP)
      G <- matrix(rgamma(n * d, shape = disp, scale = as.vector(mu) / disp), n, d)
      X <- matrix(rpois(n * d, lambda = as.vector(G)), n, d)
    } else {
      LP <- signal_X + eps
      X <- LP
    }
    
    # Second-stage
    # LF <- Z %*% betaU %*% D
    # signal_Y <- LF %*% alpha + Q %*% eta + lambda_Y * Uc
    signal_Y <- signal_X %*% alpha + Q %*% eta + lambda_Y * Uc
    delta_raw <- rnorm(n, 0, 1)
    delta <- scale_to_snr(signal_Y, delta_raw, snr_Y)
    Y <- signal_Y + delta
    
  } else if (endog_mode == "correlated_error") {
    # Correlated error mode
    e1 <- matrix(rnorm(n * 1), n, 1)
    e2 <- matrix(rnorm(n * 1), n, 1)
    e3 <- matrix(rnorm(n * d), n, d)
    
    eps_raw <- sqrt(tem) * e1 %*% matrix(1, 1, d) + sqrt(1 - tem) * e3
    delta_raw <- sqrt(tem) * e1 + sqrt(1 - tem) * e2
    
    signal_X <- Z %*% beta + Q %*% gam
    scale_to_snr <- function(signal, noise, snr_target) {
      var_signal <- mean(signal^2)
      var_noise  <- mean(noise^2)
      noise * sqrt(var_signal / (snr_target * var_noise))
    }
    eps <- scale_to_snr(signal_X, eps_raw, snr_X)
    LP <- signal_X + eps
    
    # Generate X (Gaussian or NB)
    if (tolower(family) == "nb" || tolower(family) == "negative_binomial") {
      mu <- 10*exp(LP)
      G <- matrix(rgamma(n * d, shape = disp, scale = as.vector(mu) / disp), n, d)
      X <- matrix(rpois(n * d, lambda = as.vector(G)), n, d)
      # LF <- Z %*% betaU %*% D
      # signal_Y <- LF %*% alpha + Q %*% eta
      signal_Y <- signal_X %*% alpha + Q %*% eta
      delta <- scale_to_snr(signal_Y, delta_raw, snr_Y)
      Y <- signal_Y + delta
    } else {
      X <- LP
      # LF <- Z %*% betaU %*% D
      # signal_Y <- LF %*% alpha + Q %*% eta
      signal_Y <- signal_X %*% alpha + Q %*% eta
      delta <- scale_to_snr(signal_Y, delta_raw, snr_Y)
      Y <- signal_Y + delta
    }
    

    
    Uc <- NULL
  } else {
    stop("endog_mode must be either 'confounder' or 'correlated_error'")
  }
  
  # -------------------------------------------
  # 3. Compute empirical SNRs
  # -------------------------------------------
  snr_X_emp <- mean(signal_X^2) / mean(eps^2)
  snr_Y_emp <- mean(signal_Y^2) / mean(delta^2)
  
  # -------------------------------------------
  # 4. Diagnostics
  # -------------------------------------------
  cat("Simulation Summary\n",
      "-------------------\n",
      sprintf("Family: %s | Endog. mode: %s\n", family, endog_mode),
      sprintf("n=%d, d=%d, m=%d, c=%d, k=%d\n", n, d, m, c, k),
      sprintf("Target SNR_X=%.2f → achieved %.2f\n", snr_X, snr_X_emp),
      sprintf("Target SNR_Y=%.2f → achieved %.2f\n", snr_Y, snr_Y_emp),
      sprintf("tem (endogeneity strength)=%.2f\n", tem))
  
  # -------------------------------------------
  # 5. Return structured output
  # -------------------------------------------
  list(
    Z = Z, Q = Q, X = X, Y = Y,
    beta = beta, alpha = alpha,
    gam = gam, eta = eta,
    betaU = betaU, betaV = betaV, D = D,
    eps = eps, delta = delta,
    signal_X = signal_X, signal_Y = signal_Y,
    snr_X_emp = snr_X_emp,
    snr_Y_emp = snr_Y_emp,
    Uc = Uc,
    config = list(
      n = n, d = d, m = m, c = c, k = k,
      snr_X = snr_X, snr_Y = snr_Y, tem = tem,
      family = family, endog_mode = endog_mode
    )
  )
}

# helper for NULL names
`%||%` <- function(a, b) if (!is.null(a)) a else b

