# ============================================================
# Residual Diagnostics for Logistic Regression Models
# Each predictor vs DEATH_EVENT
# Checks: Linearity, Independence, Normality, Constant Variance
# ============================================================

library(tidyverse)

df <- read.csv('/Users/kaisarthik/Library/CloudStorage/OneDrive-purdue.edu/25-26/STAT 512/STAT-512-Final-Project/heart_failure_clinical_records_dataset 2.csv')

# ============================================================
# Helper: produce 4 diagnostic plots + binned residual plot
# for a given glm model
# ============================================================
run_diagnostics <- function(model, title_prefix) {
  
  par(mfrow = c(2, 2))
  
  # --- 1. Residuals vs Fitted (Linearity + Constant Variance) ---
  plot(fitted(model), residuals(model, type = "deviance"),
       xlab = "Fitted Values (Predicted Probabilities)",
       ylab = "Deviance Residuals",
       main = paste(title_prefix, "— Residuals vs Fitted"),
       pch = 20, col = rgb(0, 0, 0, 0.4))
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(fitted(model), residuals(model, type = "deviance")),
        col = "blue", lwd = 2)
  # Interpretation: The loess curve should stay near 0. 
  # Systematic patterns suggest non-linearity.
  # Fanning/funneling suggests non-constant variance.
  
  # --- 2. Q-Q Plot (Normality of Residuals) ---
  qqnorm(residuals(model, type = "deviance"),
         main = paste(title_prefix, "— Normal Q-Q Plot"),
         pch = 20, col = rgb(0, 0, 0, 0.4))
  qqline(residuals(model, type = "deviance"), col = "red", lwd = 2)
  # Interpretation: Points should follow the line.
  # Note: normality of residuals is not strictly required for
  # logistic regression (binary outcomes produce non-normal
  # residuals by nature), but this is shown for completeness.
  
  # --- 3. Scale-Location (Constant Variance) ---
  std_resid <- rstandard(model, type = "deviance")
  plot(fitted(model), sqrt(abs(std_resid)),
       xlab = "Fitted Values",
       ylab = expression(sqrt("|Standardized Residuals|")),
       main = paste(title_prefix, "— Scale-Location"),
       pch = 20, col = rgb(0, 0, 0, 0.4))
  lines(lowess(fitted(model), sqrt(abs(std_resid))),
        col = "blue", lwd = 2)
  # Interpretation: A roughly horizontal loess line indicates
  # constant variance. Trends indicate heteroscedasticity.
  
  # --- 4. Cook's Distance (Influential Observations) ---
  plot(cooks.distance(model), type = "h",
       main = paste(title_prefix, "— Cook's Distance"),
       ylab = "Cook's Distance", xlab = "Observation Index",
       col = rgb(0, 0, 0, 0.5))
  abline(h = 4 / length(cooks.distance(model)), col = "red", lty = 2, lwd = 2)
  # Interpretation: Points above the dashed line (4/n threshold)
  # are potentially influential observations.
  
  par(mfrow = c(1, 1))
}

# ============================================================
# Helper: Binned Residual Plot (preferred linearity check
# for logistic regression with continuous predictors)
# ============================================================
binned_residual_plot <- function(model, predictor_vals, predictor_name, n_bins = 30) {
  resid_vals <- residuals(model, type = "response")
  
  bins <- cut(predictor_vals, breaks = n_bins)
  bin_df <- data.frame(predictor = predictor_vals, resid = resid_vals, bin = bins)
  
  bin_summary <- bin_df %>%
    group_by(bin) %>%
    summarise(
      mean_x = mean(predictor, na.rm = TRUE),
      mean_resid = mean(resid, na.rm = TRUE),
      se = sd(resid, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  plot(bin_summary$mean_x, bin_summary$mean_resid,
       xlab = predictor_name, ylab = "Mean Residual",
       main = paste("Binned Residual Plot —", predictor_name),
       pch = 19, col = "steelblue", cex = 1.3)
  abline(h = 0, col = "red", lwd = 2)
  # Add approximate 95% bounds
  n <- nrow(bin_df)
  abline(h = c(-2 / sqrt(n), 2 / sqrt(n)), col = "grey50", lty = 2)
  legend("topright", legend = c("Bin means", "Zero line", "±2/√n bounds"),
         col = c("steelblue", "red", "grey50"), lty = c(NA, 1, 2),
         pch = c(19, NA, NA), cex = 0.8)
  # Interpretation: Most bin means should fall within the grey
  # bounds. A clear trend indicates non-linearity in the
  # predictor-logit relationship.
}

# ============================================================
# Helper: Independence check — residuals vs row index
# ============================================================
independence_plot <- function(model, title_prefix) {
  plot(1:length(residuals(model)), residuals(model, type = "deviance"),
       xlab = "Observation Order", ylab = "Deviance Residuals",
       main = paste(title_prefix, "— Residuals vs Order (Independence)"),
       pch = 20, col = rgb(0, 0, 0, 0.4))
  abline(h = 0, col = "red", lwd = 2)
  lines(lowess(1:length(residuals(model)), residuals(model, type = "deviance")),
        col = "blue", lwd = 2)
  # Interpretation: No trend or pattern should be visible.
  # A systematic pattern (wave, drift) suggests observations
  # are not independent (e.g., time-based correlation).
}


# ============================================================
# 1. AGE vs DEATH_EVENT (Sai)
# ============================================================
cat("\n========== AGE vs DEATH_EVENT ==========\n")
model_age <- glm(DEATH_EVENT ~ age, data = df, family = binomial)
summary(model_age)

run_diagnostics(model_age, "Age")
binned_residual_plot(model_age, df$age, "Age")
independence_plot(model_age, "Age")


# ============================================================
# 2. SEX vs DEATH_EVENT (Jake)
# ============================================================
cat("\n========== SEX vs DEATH_EVENT ==========\n")
model_sex <- glm(DEATH_EVENT ~ sex, data = df, family = binomial)
summary(model_sex)

run_diagnostics(model_sex, "Sex")
# Note: Linearity/binned residual plot is not meaningful for a
# binary predictor. The standard 4-panel diagnostics suffice.
independence_plot(model_sex, "Sex")


# ============================================================
# 3. PLATELETS vs DEATH_EVENT (Jad)
# ============================================================
cat("\n========== PLATELETS vs DEATH_EVENT ==========\n")
model_plt <- glm(DEATH_EVENT ~ platelets, data = df, family = binomial)
summary(model_plt)

run_diagnostics(model_plt, "Platelets")
binned_residual_plot(model_plt, df$platelets, "Platelets")
independence_plot(model_plt, "Platelets")


# ============================================================
# 4. SERUM CREATININE vs DEATH_EVENT (Alex/Jinwook)
# ============================================================
cat("\n========== SERUM CREATININE vs DEATH_EVENT ==========\n")
model_sc <- glm(DEATH_EVENT ~ serum_creatinine, data = df, family = binomial)
summary(model_sc)

run_diagnostics(model_sc, "Serum Creatinine")
binned_residual_plot(model_sc, df$serum_creatinine, "Serum Creatinine")
independence_plot(model_sc, "Serum Creatinine")


# ============================================================
# 5. EJECTION FRACTION vs DEATH_EVENT (Lohit)
# ============================================================
cat("\n========== EJECTION FRACTION vs DEATH_EVENT ==========\n")
model_ef <- glm(DEATH_EVENT ~ ejection_fraction, data = df, family = binomial)
summary(model_ef)

run_diagnostics(model_ef, "Ejection Fraction")
binned_residual_plot(model_ef, df$ejection_fraction, "Ejection Fraction")
independence_plot(model_ef, "Ejection Fraction")


# ============================================================
# INTERPRETATION GUIDE (printed to console)
# ============================================================
cat("
============================================================
HOW TO INTERPRET EACH DIAGNOSTIC PLOT
============================================================

1. RESIDUALS vs FITTED (Linearity + Constant Variance)
   - Loess line near zero → linearity holds
   - Random scatter → constant variance holds
   - Fanning / funneling → heteroscedasticity

2. NORMAL Q-Q PLOT (Normality)
   - Points on the line → residuals approx. normal
   - Note: for logistic regression, residuals are inherently
     non-normal due to the binary outcome. Heavy tails or
     S-shapes are expected and not necessarily a concern.

3. SCALE-LOCATION (Constant Variance)
   - Horizontal loess → homoscedasticity
   - Upward/downward trend → variance changes with fitted values

4. COOK'S DISTANCE (Influential Points)
   - Points above 4/n line are potentially influential
   - Investigate if removing them changes conclusions

5. BINNED RESIDUAL PLOT (Linearity for logistic regression)
   - Most points within the +/- 2/sqrt(n) bounds → good fit
   - Systematic curve → non-linear relationship; consider
     adding polynomial or spline terms

6. RESIDUALS vs ORDER (Independence)
   - No trend or wave pattern → independence holds
   - Systematic pattern → possible autocorrelation
")