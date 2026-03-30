# ============================================================
# Heart Failure EDA
# Variables: Age, Sex, Platelets, Serum Creatinine, 
#            Ejection Fraction vs DEATH_EVENT
# ============================================================

library(tidyverse)

# --- Read Data -----------------------------------------------
df <- read.csv('/Users/kaisarthik/Library/CloudStorage/OneDrive-purdue.edu/25-26/STAT 512/STAT-512-Final-Project/heart_failure_clinical_records_dataset 2.csv')

df$death_label <- factor(df$DEATH_EVENT, levels = c(0, 1), labels = c("Survived", "Died"))
df$sex_label   <- factor(df$sex, levels = c(0, 1), labels = c("Female", "Male"))

# --- Dataset Overview ----------------------------------------
str(df)
summary(df[, c("age", "sex", "platelets", "serum_creatinine", 
               "ejection_fraction", "DEATH_EVENT")])
colSums(is.na(df))  # check missing values
table(df$DEATH_EVENT)

# --- Summary Stats by Death Event ----------------------------
df %>%
  group_by(death_label) %>%
  summarise(
    across(c(age, platelets, serum_creatinine, ejection_fraction),
           list(mean = mean, sd = sd, median = median, min = min, max = max),
           .names = "{.col}_{.fn}")
  ) %>%
  print(width = Inf)

# --- Sex vs Death Event Crosstab -----------------------------
table(df$sex_label, df$death_label)
prop.table(table(df$sex_label, df$death_label), margin = 1)

# ============================================================
# PLOTS
# ============================================================

# Color palette
cols <- c("Survived" = "#2ecc71", "Died" = "#e74c3c")

# --- 1. Age vs Death Event (Sai) ----------------------------
ggplot(df, aes(x = death_label, y = age, fill = death_label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Age vs Death Event", x = "", y = "Age") +
  theme_minimal() +
  theme(legend.position = "none")

# Overlapping density
ggplot(df, aes(x = age, fill = death_label)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Age Distribution by Outcome", x = "Age", y = "Density", fill = "Outcome") +
  theme_minimal()

# --- 2. Sex vs Death Event (Jake) ---------------------------
sex_pct <- df %>%
  count(sex_label, death_label) %>%
  group_by(sex_label) %>%
  mutate(pct = n / sum(n) * 100)

ggplot(sex_pct, aes(x = sex_label, y = pct, fill = death_label)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", pct)), 
            position = position_dodge(0.6), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Sex vs Death Event", x = "", y = "Percentage (%)", fill = "Outcome") +
  theme_minimal()

# --- 3. Platelets vs Death Event (Jad) ----------------------
ggplot(df, aes(x = death_label, y = platelets, fill = death_label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Platelets vs Death Event", x = "", y = "Platelets (kiloplatelets/mL)") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(df, aes(x = platelets, fill = death_label)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Platelets Distribution by Outcome", x = "Platelets", y = "Density", fill = "Outcome") +
  theme_minimal()

# --- 4. Serum Creatinine vs Death Event (Alex/Jinwook) ------
ggplot(df, aes(x = death_label, y = serum_creatinine, fill = death_label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Serum Creatinine vs Death Event", x = "", y = "Serum Creatinine (mg/dL)") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(df, aes(x = serum_creatinine, fill = death_label)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Serum Creatinine Distribution by Outcome", 
       x = "Serum Creatinine (mg/dL)", y = "Density", fill = "Outcome") +
  theme_minimal()

# --- 5. Ejection Fraction vs Death Event (Lohit) -----------
ggplot(df, aes(x = death_label, y = ejection_fraction, fill = death_label)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 1.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Ejection Fraction vs Death Event", x = "", y = "Ejection Fraction (%)") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(df, aes(x = ejection_fraction, fill = death_label)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = cols) +
  labs(title = "Ejection Fraction Distribution by Outcome", 
       x = "Ejection Fraction (%)", y = "Density", fill = "Outcome") +
  theme_minimal()

# --- 6. Correlation Heatmap ---------------------------------
# install.packages("corrplot") if needed
library(corrplot)

cor_mat <- cor(df[, c("age", "platelets", "serum_creatinine", 
                      "ejection_fraction", "DEATH_EVENT")])
corrplot(cor_mat, method = "color", type = "upper", addCoef.col = "black",
         tl.col = "black", tl.srt = 45, title = "Correlation Heatmap",
         mar = c(0, 0, 2, 0))

# ============================================================
# PRELIMINARY STATISTICAL TESTS
# ============================================================

# Age: two-sample t-test
t.test(age ~ DEATH_EVENT, data = df)

# Sex: chi-squared test
chisq.test(table(df$sex, df$DEATH_EVENT))

# Platelets: two-sample t-test
t.test(platelets ~ DEATH_EVENT, data = df)

# Serum Creatinine: two-sample t-test
t.test(serum_creatinine ~ DEATH_EVENT, data = df)

# Ejection Fraction: two-sample t-test
t.test(ejection_fraction ~ DEATH_EVENT, data = df)