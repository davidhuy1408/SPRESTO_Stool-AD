# ==============================================================================
# Longitudinal Analysis: GLMM for Binary Outcome (Control vs AD)
# ==============================================================================

# 1. Load necessary libraries
library(lme4)
library(dplyr)
library(broom.mixed)

# 2. Load Data
# ------------------------------------------------------------------------------
data <- read.csv("your_data.csv")

# 3. Define Variable Names
# ------------------------------------------------------------------------------
outcome_col <- "Group" 
predictors_to_test <- c("Bifidobacterium", "Lactobacillus", "IL_6", "Butyrate")
time_col <- "TimePoint" # Important for longitudinal analysis
subject_id_col <- "SubjectID"
covariate_1 <- "Mat_History_Allergy"
covariate_2 <- "Mat_Antibiotics"
covariate_3 <- "Mode_of_Delivery"

# 4. Data Preparation
# ------------------------------------------------------------------------------
data[[outcome_col]] <- as.factor(data[[outcome_col]])
data[[subject_id_col]] <- as.factor(data[[subject_id_col]])
data[[covariate_1]] <- as.factor(data[[covariate_1]])
data[[covariate_2]] <- as.factor(data[[covariate_2]])
data[[covariate_3]] <- as.factor(data[[covariate_3]])

results_table <- data.frame()

# 5. Run GLMM Loop
# ------------------------------------------------------------------------------
for (predictor in predictors_to_test) {

  print(paste("Analyzing predictor:", predictor))
  
  formula_str <- paste(outcome_col, "~", predictor, "+", 
                       time_col, "+",
                       covariate_1, "+", 
                       covariate_2, "+", 
                       covariate_3, "+", 
                       "(1|", subject_id_col, ")")
  
  f <- as.formula(formula_str)
  
  tryCatch({
    model <- glmer(f, data = data, family = binomial, nAGQ = 0)
    
    model_res <- tidy(model, exponentiate = TRUE, conf.int = TRUE) %>% 
      filter(term == predictor) %>% 
      mutate(Predictor_Variable = predictor) %>%
      select(Predictor_Variable, term, estimate, std.error, p.value, conf.low, conf.high)
    
    colnames(model_res)[colnames(model_res) == "estimate"] <- "Odds_Ratio"
    
    results_table <- rbind(results_table, model_res)
    
  }, error = function(e) {
    message(paste("Error fitting model for:", predictor))
    message(e)
  })
}

# 6. Apply Benjamini–Hochberg FDR Correction
# ------------------------------------------------------------------------------
if (nrow(results_table) > 0) {
  
  results_table$FDR_adj_p <- p.adjust(results_table$p.value, method = "BH")
  
  results_table$Significant_10pct_FDR <- ifelse(results_table$FDR_adj_p < 0.10, "Yes", "No")
  
  # 7. View and Save Results
  # ------------------------------------------------------------------------------
  print(results_table)
  