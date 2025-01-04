##################################
##########   STAT 409   ##########
####  Data Science Practicum  ####
##################################
####### LOGISTIC REGRESSION ######
###########  LOGLOSS  ############
##################################

##  NOTE: Run 'Data_Cleaning_and_Integration.R' first

# Load necessary libraries
library(dplyr)
library(caret)
library(car)
library(foreach)
library(doParallel)

################################################################################

# Define dataset and response variable
data <- data_time_merge
response_var <- "amyloid_binary"

# Additional predictors to always include in the model
additional_predictors <-
  c("age_at_appointment",
    "time_diff",
    "sex",
    "Calculated_Consensus_dx")

################################################################################

# Function to generate all valid predictor combinations
generate_combinations <- function(predictors, must_include = NULL) {
  # Generate all possible combinations of predictors
  predictor_combinations <-
    unlist(lapply(1:length(predictors), function(x) {
      combn(predictors, x, simplify = FALSE)
    }), recursive = FALSE)
  
  # Filter combinations to include at least one plasma predictor
  valid_combinations <- Filter(function(comb) {
    any(comb %in% plasma_predictors)
  }, predictor_combinations)
  
  # If a must_include variable is specified, ensure it is in each combination
  if (!is.null(must_include)) {
    valid_combinations <-
      lapply(valid_combinations, function(comb)
        c(comb, must_include))
  }
  
  return(valid_combinations)
}

################################################################################

# Function to calculate log-loss for logistic regression using LOOCV
calculate_logloss_loocv <-
  function(data, response_var, predictors) {
    n <- nrow(data)
    logloss_values <- numeric(n)
    
    for (i in 1:n) {
      train_data <- data[-i, ]
      test_data <- data[i, , drop = FALSE]
      
      if (length(predictors) == 1) {
        formula_str <- paste(response_var, "~", predictors)
      } else {
        formula_str <-
          paste(response_var, "~", paste(predictors, collapse = " + "))
      }
      
      # Debugging: Print the formula being used
      print(paste("Trying formula:", formula_str))
      
      fit <-
        eval(parse(
          text = paste("glm(", formula_str, ", 
                       data = train_data, family = binomial)")
        ))
      
      prediction <- predict(fit, test_data, type = "response")
      actual <- test_data[[response_var]]
      
      logloss_values[i] <-
        -(actual * log(prediction) + (1 - actual) * log(1 - prediction))
    }
    
    logloss <- mean(logloss_values)
    return(logloss)
  }

################################################################################

# Function to perform LOOCV for predictor combos & calculate log loss using glm
perform_logloss <- function(combinations, data, response_var) {
  # Initialize a data frame to store the results
  results <-
    data.frame(
      model = character(),
      logloss = numeric(),
      stringsAsFactors = FALSE
    )
  
  # Perform LOOCV for each combination of predictors in parallel
  results <-
    foreach(comb = combinations,
            .combine = rbind,
            .packages = c("dplyr")) %dopar% {
              logloss <- calculate_logloss_loocv(data, response_var, comb)
              
              formula_str <- if (length(comb) == 1) {
                paste(response_var, "~", comb)
              } else {
                paste(response_var, "~", paste(comb, collapse = " + "))
              }
              
              # Return the model formula and its log loss as a data frame row
              data.frame(model = formula_str,
                         logloss = logloss,
                         stringsAsFactors = FALSE)
            }
  return(results)
}

################################################################################

# Set up parallel backend to use the available cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Export the necessary variables and functions to the cluster
clusterExport(
  cl,
  varlist = c(
    "data",
    "response_var",
    "plasma_predictors",
    "additional_predictors",
    "calculate_logloss_loocv"
  )
)

# Generate combinations for each scenario
combinations_apoe_ternary <-
  generate_combinations(c(plasma_predictors, additional_predictors), 
                        "apoe_ternary")
combinations_apoe_binary <-
  generate_combinations(c(plasma_predictors, additional_predictors), 
                        "apoe_binary")
combinations_neither <-
  generate_combinations(c(plasma_predictors, additional_predictors))

################################################################################

# Perform log loss using LOOCV for models with 'apoe_ternary'
results_logreg_tern <-
  perform_logloss(combinations_apoe_ternary, data, response_var)

# Find the best model based on the lowest log loss using LOOCV
best_model_logreg_tern <-
  results_logreg_tern %>% filter(logloss == min(logloss))

# Print the best model
print("Best model with apoe_ternary:")
print(best_model_logreg_tern)

# Display top 10 models with the lowest log loss scores
top_10_logreg_tern <-
  results_logreg_tern %>% arrange(logloss) %>% head(10)
print("Top 10 models with the lowest log loss scores (apoe_ternary):")
print(top_10_logreg_tern)

################################################################################

# Perform log loss for models with 'apoe_binary'
results_logreg_bin <-
  perform_logloss(combinations_apoe_binary, data, response_var)

# Find the best model based on the lowest log loss
best_model_logreg_bin <-
  results_logreg_bin %>% filter(logloss == min(logloss))

# Print the best model
print("Best model with apoe_binary:")
print(best_model_logreg_bin)

# Display top 10 models with the lowest log loss scores
top_10_logreg_bin <-
  results_logreg_bin %>% arrange(logloss) %>% head(10)
print("Top 10 models with the lowest log loss scores (apoe_binary):")
print(top_10_logreg_bin)

################################################################################

# Perform log loss for models without 'apoe_binary' or 'apoe_ternary'
results_logreg_neither <-
  perform_logloss(combinations_neither, data, response_var)

# Find the best model based on the lowest log loss
best_model_logreg_neither <-
  results_logreg_neither %>% filter(logloss == min(logloss))

# Print the best model
print("Best model without apoe:")
print(best_model_logreg_neither)

# Display top 10 models with the lowest log loss scores
top_10_logreg_neither <-
  results_logreg_neither %>% arrange(logloss) %>% head(10)
print("Top 10 models with the lowest log loss scores (no apoe):")
print(top_10_logreg_neither)

################################################################################

# Stop the parallel backend
stopCluster(cl)

# Print the best models
print("Best model with apoe_ternary:")
print(best_model_logreg_tern)

print("Best model with apoe_binary:")
print(best_model_logreg_bin)

print("Best model without apoe:")
print(best_model_logreg_neither)

################################################################################

# Write results to .csv
# Combine all results
all_results <-
  rbind(results_logreg_tern,
        results_logreg_bin,
        results_logreg_neither)

# Save combined results to CSV
write.csv(all_results, "logreg_logloss_results.csv", row.names = FALSE)
