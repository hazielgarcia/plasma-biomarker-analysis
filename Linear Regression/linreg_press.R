##################################
##########   STAT 409   ##########
####  Data Science Practicum  ####
##################################
########### SKELETON  ############
####### LINEAR REGRESSION ########
##################################

##  NOTE: Run 'Data_Cleaning_and_Integration.R' first

# Load necessary libraries
library(dplyr)
library(caret)
library(qpcR)
library(car)
library(doParallel)
library(foreach)

################################################################################

# Define dataset and response variable
data <- data_time_merge
response_var <- "Amyloid_PET"

# Additional predictors to always include in the model
additional_predictors <- 
  c("age_at_appointment", "time_diff", "sex", "Calculated_Consensus_dx")

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

# Function to perform LOOCV for all predictor combos and apoe variable using lm
perform_press <- function(combinations, data, response_var) {
  # Initialize a data frame to store the results
  results <-
    data.frame(model = character(),
               press = numeric(),
               stringsAsFactors = FALSE)
  
  
  # Perform LOOCV for each combination of predictors in parallel
  results <-
    foreach(comb = combinations,
            .combine = rbind,
            .packages = c("qpcR", "dplyr")) %dopar% {
              if (length(comb) == 1) {
                formula_str <- paste(response_var, "~", comb)
              } else {
                formula_str <-
                  paste(response_var, "~", paste(comb, collapse = " + "))
              }
              
              # Debugging: Print the formula being used
              print(paste("Trying formula:", formula_str))
              
              fit <-
                {
                  eval(parse(text = 
                               paste("lm(", formula_str, ", data = data)")))
                }
              press <- PRESS(fit)
              
              # Return the model formula and its PRESS stat as a data frame row
              data.frame(model = formula_str,
                         press = press$stat,
                         stringsAsFactors = FALSE)
            }
  return(results)
}


################################################################################

# Set up parallel backend to use the available cores
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Export the necessary variables to the cluster
clusterExport(cl,
              varlist = c(
                "data",
                "response_var",
                "plasma_predictors",
                "additional_predictors"
              ))

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

# Perform PRESS for models with 'apoe_ternary'
results_linreg_tern <- 
  perform_press(combinations_apoe_ternary, data, response_var)

# Find the best model based on the lowest PRESS
best_model_linreg_tern <- results_linreg_tern %>% filter(press == min(press))

# Print the best model
print("Best model with apoe_ternary:")
print(best_model_linreg_tern)

# Display top 10 models with the lowest PRESS scores
top_10_linreg_tern <- results_linreg_tern %>% arrange(press) %>% head(10)
print("Top 10 models with the lowest PRESS scores (apoe_ternary):")
print(top_10_linreg_tern)

################################################################################

# Perform PRESS for models with 'apoe_binary'
results_linreg_bin <- 
  perform_press(combinations_apoe_binary, data, response_var)

# Find the best model based on the lowest PRESS
best_model_linreg_bin <- results_linreg_bin %>% filter(press == min(press))

# Print the best model
print("Best model with apoe_binary:")
print(best_model_linreg_bin)


# Display top 10 models with the lowest PRESS scores
top_10_linreg_bin <- results_linreg_bin %>% arrange(press) %>% head(10)
print("Top 10 models with the lowest PRESS scores (apoe_binary):")
print(top_10_linreg_bin)


################################################################################

# Perform PRESS for models without  'apoe_binary' or 'apoe_ternary'
results_linreg_neither <- 
  perform_press(combinations_neither, data, response_var)

# Find the best model based on the lowest PRESS
best_model_linreg_neither <- 
  results_linreg_neither %>% filter(press == min(press))

# Print the best model
print("Best model without apoe:")
print(best_model_linreg_neither)


# Display top 10 models with the lowest PRESS scores
top_10_linreg_neither <- results_linreg_neither %>% arrange(press) %>% head(10)
print("Top 10 models with the lowest PRESS scores (no apoe):")
print(top_10_linreg_neither)

################################################################################

# Stop the parallel backend
stopCluster(cl)

# Print the best models
print("Best model with apoe_ternary:")
print(best_model_linreg_tern)

print("Best model with apoe_binary:")
print(best_model_linreg_bin)

print("Best model without apoe:")
print(best_model_linreg_neither)

################################################################################

# Write results to .csv
# Combine all results
all_results <- 
  rbind(results_linreg_tern, results_linreg_bin, results_linreg_neither)

# Save combined results to CSV
write.csv(all_results, "linreg_press_results.csv", row.names = FALSE)
