# Simulation Study for Factor-IV Framework

Here we provide details of simulation pipeline for evaluating **Factor-IV**, a supervised latent-factor instrumental variable (IV) method for recovering causal effects in high-dimensional biological systems.

The workflow combines:

- Data generation via a flexible IV simulation engine (`simulate_data.R`)
- Analysis and comparison across multiple scenarios (`SimulationStudy.Rmd`)
- Supporting function for output analysis (`simulate_analysis.R`)
- Storage of fitted models in `output/`
- Publication-ready figures in `plots/`

The goal is to systematically compare **Factor-IV** against naïve methods that ignore endogeneity, under a variety of realistic data-generating mechanisms.

---

## 1. Folder Structure

A typical layout of the `Simulation` folder is:

```text
Simulation/
│
├── simulate_data.R           # Data-generation engine for IV simulations
├── simulate_analysis.R       # Suuporting function for output analysis
├── SimulationStudy.Rmd       # Main simulation workflow & analysis
│
├── output/                   # Saved simulation objects and metrics
│     ├── Gaussian_confounder_spreduced.rds            # Gaussian setting with unobserved confounder
│     ├── Gaussian_correlated_error_spreduced.rds      # Gaussian setting with correlated errror 
│     ├── NB_confounder_spreduced.rds                  # Negative Binomial setting with unobserved confounder
│     ├── NB_correlated_error_spreduced.rds            # Negative Binomial setting with correlated errror 
│
└── plots/                    # Plots for comparison and diagnostics
      ├── Gaussian_confounder_alpha_comparison.pdf
      |── Gaussian_confounder_model_fit_true.pdf
      ├── Gaussian_correlated_error_alpha_comparison.pdf
      |── Gaussian_correlated_error_model_fit_true.pdf
      ├── NB_confounder_alpha_comparison.pdf
      |── NB_confounder_model_fit_true.pdf
      ├── NB_correlated_error_alpha_comparison.pdf
      |── NB_correlated_error_model_fit_true.pdf
```


---

## 2. Data-Generation Engine (`simulate_data.R`)

The file `simulate_data.R` defines a function (e.g. `simulate_iv_data()`) that generates synthetic data under a high-dimensional IV framework. The core objects are:

- **Instruments**: $Z \in \mathbb{R}^{n \times m}$
- **High-dimensional exposures** (e.g. microbes, metabolites): $X \in \mathbb{R}^{n \times d}$
- **Observed confounders**: $Q \in \mathbb{R}^{n \times c}$
- **Latent factors**: $U \in \mathbb{R}^{n \times k}$ encoded via sparse loadings
- **Outcome**: $Y \in \mathbb{R}^n$

### 2.1 Latent Factor & First-Stage Model

A supervised factor structure is imposed via sparse matrices $\( U \)$ and $\( V \)$ and a diagonal matrix $\( D \)$:

```math
X = Z U D V^\top + Q \Gamma + \varepsilon,
```

where:

- $U \in \mathbb{R}^{m \times k}$ and $V \in \mathbb{R}^{d \times k}$ are sparse and column-normalized,
- $ D $ controls the strength of each latent factor,
- $\Gamma \in \mathbb{R}^{c \times d}$ describes the contribution of observed confounders $Q$,
- $\varepsilon$ is noise scaled to achieve a prescribed **signal-to-noise ratio (SNR)** in the first stage.

### 2.2 Outcome Model

The causal signal is carried through the instrument–factor path:

```math
\text{LF} = Z U D, \qquad
Y = \text{LF} \alpha + Q \eta + \delta,
```

where:

- $\alpha \in \mathbb{R}^k$ are the **true causal effects of the latent factors** (sparse or dense, depending on the setting),
- $\( \eta \in \mathbb{R}^c \)$ encodes the effects of confounders,
- $\( \delta \)$ is outcome noise scaled to match target $\( \text{SNR}_Y \)$.

Note that in the high-dimensional setting, the “true” effect of each individual covariate in $\( X \)$ can be thought of as embedded in $\( V \alpha \)$, though the simulation is parameterized at the factor level.

### 2.3 Endogeneity Mechanisms

Two sources of endogeneity can be toggled via arguments:

1. **Correlated Error (`endog_mode = "correlated_error"`)**

   A shared error component drives both $\( X \)$ and $\( Y \)$:

```math
   \varepsilon = \sqrt{\text{tem}}\, e_1 \mathbf{1}_d^\top + \sqrt{1-\text{tem}}\, e_3,
   \quad
   \delta = \sqrt{\text{tem}}\, e_1 + \sqrt{1-\text{tem}}\, e_2,
```

   where $\( e_1 \)$ and $\( e_2 \)$ are independent standard normal random variables, and $\( e_3 \)$ is $\( e_1 \)$ plus $\( e_2 \)$.

   where `tem` controls the strength of the correlation (endogeneity).

2. **Latent Confounder (`endog_mode = "confounder"`)**

   A shared unobserved confounder $\( U_c \)$ directly affects both $\( X \)$ and $\( Y \)$:

  ```math
   X = Z\beta + Q\Gamma + \lambda_X U_c + \varepsilon,
   \quad
   Y = \text{LF} \alpha + Q\eta + \lambda_Y U_c + \delta.
  ```

   Here `tem` (or derived scaling parameters) controls how strongly $\( U_c \)$ induces endogeneity.

### 2.4 SNR Control

The function rescales the noise terms so that:

  ```math
\text{SNR}_X
\approx 
\frac{\mathbb{E}[\text{signal}_X^2]}{\mathbb{E}[\varepsilon^2]}, \qquad
\text{SNR}_Y
\approx
\frac{\mathbb{E}[\text{signal}_Y^2]}{\mathbb{E}[\delta^2]}
```

match user-specified target values `snr_X` and `snr_Y`.  
Empirical SNRs (computed from the simulated sample) are returned for diagnostics.

### 2.5 Distributional Family (Gaussian vs NB)

- **Gaussian**: `family = "Gaussian"`  
  $\( X \)$ is continuous, using the linear model as written above.

- **Negative Binomial**: `family = "NB"`  
  Each entry of $\(X\)$ is generated from an NB model:

```math
  \mu = \exp(LP), \quad
  X_{ij} \sim \text{NB}(\mu_{ij}, \text{disp}),
  ```

  where `disp` controls overdispersion and `LP` is the linear predictor.

The returned object typically includes:

- `Z, Q, X, Y`
- `betaU, betaV, D, beta`
- `alpha, gam, eta`
- `eps, delta`
- `signal_X, signal_Y`
- empirical `snr_X_emp`, `snr_Y_emp`
- configuration list: dimensions, SNR, family, endogeneity mode, etc.

---

## 3. Simulation Study (`SimulationStudy.Rmd`)

`SimulationStudy.Rmd` orchestrates full  experiments:

### 3.1 Study Design

For each scenario, the RMarkdown typically varies:

- **Dimensions**: $\(n, d, m, c, k\)$
- **Endogeneity strength**: `tem`
- **SNR levels**: `snr_X`, `snr_Y`
- **Family**: `"Gaussian"` vs `"NB"`
- **Endogeneity mode**: `"correlated_error"` vs `"confounder"`

A loop over combinations of these settings:

1. Calls `simulate_iv_data(...)` (from `simulate_data.R`)
2. Fits different estimators of the causal effect
3. Computes performance metrics
4. Stores fitted objects and summary statistics in `output/`
5. Generates comparison plots in `plots/`

### 3.2 Naïve Estimator (Ignoring Instruments)

The naïve approach regresses **each exposure on the outcome and confounders**, effectively ignoring the IV structure:

- For Gaussian $\(X\)$:  
  Linear regression / multiple regression:
  ```math
  X_j \sim Y + Q.
  ```

- For NB/count $\(X\)$:  
  DESeq2 vignette-style models:
    ```math
  X_j \sim Y + Q
  ```
  with NB likelihood and dispersion estimation.

From these models, we extract:

- naive coefficient $\( \hat{\alpha}^{\text{naive}}_j \)$
- p-values and FDR-adjusted p-values
- set coefficients to 0 if not significant (adjusted p-value > threshold).

### 3.3 Factor-IV Estimator

The Factor-IV pipeline typically proceeds through:

1. **Supervised factor model** (e.g., GOFAR or related methods):  
     ```math
   X \approx U D V^\top
   ```
   where $\(U\)$ is linked to $\(Z\)$ and $\(Q\)$.

2. **Instrumental-variables outcome model** using the learned latent factors (`U_hat`) and instruments (`Z`):

   - 2SLS / GMM / IV-probit depending on the design.  
   - Extract latent-factor-level estimates $\( \hat{\alpha}^{\text{IV}} \)$.

3. Optionally, transform back to feature space using $\(V\)$ for per-feature causal interpretation.

### 3.4 Performance Metrics

For each scenario, the simulation computes:

- **Mean Squared Error (MSE)** of α-estimates:
    ```math
  \text{MSE}(\hat{\alpha}) = \frac{1}{k}\sum_{i=1}^k (\hat{\alpha}_i - \alpha_i)^2
  ```
  (or in feature space, using $\(V\alpha\)$ if desired),

- **R²** from regressing estimated α on true α:
  - Factor-IV vs True
  - Naïve vs True

- **Support recovery** summaries (TP, FP, FN) if α is sparse.

These metrics are saved in the `output/` folder for later aggregation.

---

## 4. Plots (`plots/`)

The `plots/` folder contains figures that visually compare methods across scenarios.

Typical plots include:

### 4.1 Scatterplots: True vs Estimated Effects

Faceted comparison of **Factor-IV** and **Naïve OLS/GLM**:

- X-axis: True α  
- Y-axis: Estimated α  
- Per-method panels (Factor-IV, Naïve)  
- R² and MSE annotated in titles or subtitles.

These show how closely each estimator recovers the ground-truth signal.

### 4.2 SNR Sensitivity

Plots illustrating how performance changes as:

- `snr_X` and `snr_Y` vary,
- or as `tem` increases/decreases (endogeneity strength).

Examples:

- MSE vs SNR curves
- Heatmaps over `(snr_X, snr_Y)` combinations.

### 4.3 Endogeneity Effect

Diagnostic plots showing:

- correlation between a representative $\( X_j \)$ and $\(Y\)$,
- correlation between a representative $\( Z_j \)$ and$ \(Y\)$,
- to visually confirm endogeneity and the relevance of instruments.

---

## 5. Running the Entire Simulation

To reproduce the simulation study:

1. Ensure the working directory is the `Simulation` folder.
2. Open `SimulationStudy.Rmd` in RStudio (or call from terminal):
   ```r
   rmarkdown::render("SimulationStudy.Rmd")
   ```
3. This will:
   - source `simulate_data.R`, `simulate_analysis.R`,
   - generate data for each scenario,
   - fit models (naïve & Factor-IV),
   - save results to `output/`,
   - generate plots in `plots/`.

---

## 6. Interpreting Results

### 6.1 When Factor-IV Should Outperform Naïve Methods

- **Non-zero endogeneity** (`tem` moderately large),
- Strong instrument–latent factor relationship,
- Moderate or high noise in $\(X\)$ or $\(Y\)$,
- Sparse signals in $\( \alpha \)$ and structured loading matrices.

Under these settings, naïve methods are biased, while Factor-IV recovers a consistent estimate (up to finite-sample error and factor estimation error).

### 6.2 When Naïve Methods May Look Competitive

- When endogeneity is very weak (`tem` near 0),
- When instruments are weak or nearly irrelevant,
- In very low-noise settings with simple confounding.

Simulation results across SNR and `tem` help illustrate these transitions.

---

## 7. Extending the Simulation

You can:

- Modify `simulate_data.R` to:
  - change sparsity of `betaU`, `betaV`, `alpha`,
  - alter distributions of `Z` and `Q`,
  - add additional outcome types (binary, survival, etc.).

- Extend `SimulationStudy.Rmd` to:
  - add more methods (e.g., ridge IV, Lasso-2SLS),
  - include power and Type I error analyses,
  - generate summary tables for inclusion in the manuscript.

---

## 8. Reproducibility and Random Seeds

All simulations are set up with explicit `seed` parameters inside the data-generation function to ensure reproducibility. You can fix or vary the seed per scenario to:

- reproduce exact results reported in the manuscript, or
- perform additional Monte Carlo runs for robustness.

---

## 9. Contact / Citation

If you use this simulation framework in your own work or build upon it, please cite the associated manuscript:

> **Factor Instrumental Variables Analysis: A New Approach for Uncovering Cause-and-Effect Relationships in Disease from Multi-Omics Data.**

and acknowledge the simulation framework from this repository.
