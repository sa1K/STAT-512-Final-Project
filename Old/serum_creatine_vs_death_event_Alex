---
title: "Multiple Logistic Regression: Effect of Serum Creatinine on DEATH_EVENT"
output: html_document
---

# Setup
df <- read.csv("/Users/alexjung/Downloads/heart_failure_clinical_records_dataset 2.csv")

# Full Multiple Logistic Regression


sc_model <- glm(
  DEATH_EVENT ~ serum_creatinine + age + anaemia +
    creatinine_phosphokinase + diabetes +
    ejection_fraction + high_blood_pressure +
    platelets + serum_sodium + sex + smoking,
  data = df,
  family = binomial(link = "logit")
)

summary(sc_model)

# Null Model for Likelihood Ratio Test
null_model <- glm(DEATH_EVENT ~ 1, data = df, family = binomial)

# Compare Null vs Full Model
anova(null_model, sc_model, test = "Chisq")

# McFadden Pseudo R^2
1 - (sc_model$deviance / sc_model$null.deviance)

# Odds Ratios + 95% CI
exp(cbind(OR = coef(sc_model), confint(sc_model)))

# Reduced Model (significant predictors only)
reduced_sc_model <- glm(
  DEATH_EVENT ~ age + ejection_fraction + serum_creatinine,
  data = df,
  family = binomial(link = "logit")
)

summary(reduced_sc_model)

# Compare Reduced vs Full
anova(reduced_sc_model, sc_model, test = "Chisq")

# AIC / BIC Comparison
AIC(sc_model, reduced_sc_model)
BIC(sc_model, reduced_sc_model)
