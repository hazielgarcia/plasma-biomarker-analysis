##################################
##########   STAT 409   ##########
####  Data Science Practicum  ####
##################################
###### LINEAR REGRESSION #########
##################################


##  NOTE: Run 'Data_Cleaning_and_Integration.R' first

################################################################################
# Load necessary libraries
library(dplyr)
library(caret)
library(ggplot2)
library(tidyr)
library(knitr)
library(kableExtra)
library(MASS)

################################################################################
# Define dataset and predictors
data <- data_time_merge

# Response Variable
response_var <- "Amyloid_PET"
# Ensure the response variable is numeric
data[[response_var]] <- as.numeric(as.character(data[[response_var]]))

# Combinations of Predictors found with PRESS
best_formula_full <- 
  "Amyloid_PET ~ ptau181 + pTau217 + GFAP + NFL + Ab_ratio + sex 
                 + Calculated_Consensus_dx + apoe_ternary"
best_formula_reduced <- 
  "Amyloid_PET ~ pTau217 + GFAP + NFL + Ab_ratio + sex"

################################################################################
######################  CHECKING THE 4 ASSUMPTIONS  ############################
################################################################################

# Function to plot the model with summary
plot_model <- function(model) {
  print(summary(model))
  par(mfrow = c(3, 2))
  for (i in 1:6) {
    plot(model,
         which = i)
  }
  par(mfrow = c(1, 1))
}

# Create initial models for analysis
full_model <-
  eval(parse(text = paste("lm(", best_formula_full, ", data = data)")))
reduced_model <-
  eval(parse(text = paste("lm(", best_formula_reduced, ", data = data)")))


# Linearity
# Full
plot(
  model.matrix(full_model) %*% coef(full_model),
  data$Amyloid_PET,
  xlab = "β * x",
  ylab = "y",
  main = "Scatterplot of y vs. β * x"
)
abline(a = 0, b = 1, col = "red")

# Reduced
# Full
plot(
  model.matrix(reduced_model) %*% coef(reduced_model),
  data$Amyloid_PET,
  xlab = "β * x",
  ylab = "y",
  main = "Scatterplot of y vs. β * x"
)
abline(a = 0, b = 1, col = "red")

# Relatively descent linearity is achieved


plot_model(full_model)
plot_model(reduced_model)


# Log transformation of the response variable
data$log_Amyloid_PET <- log(data$Amyloid_PET)

# Full Model
log_model_full <-
  lm(
    log_Amyloid_PET ~ ptau181 + pTau217 + GFAP + NFL + Ab_ratio + sex
    + Calculated_Consensus_dx + apoe_ternary,
    data = data
  )

# Reduced Model
log_model_reduced <-
  lm(
    log_Amyloid_PET ~ pTau217 + GFAP + NFL + Ab_ratio + sex,
    data = data
  )

plot_model(log_model_full)
plot_model(log_model_reduced)



# Finding the optimal transformation using Box-Cox
# Full model BoxCox
boxcox_model_full <-
  boxcox(
    lm(
      Amyloid_PET ~ ptau181 + pTau217 + GFAP + NFL + Ab_ratio + sex
      + Calculated_Consensus_dx + apoe_ternary,
      data = data
    )
  )
best_lambda_full <-
  boxcox_model_full$x[which.max(boxcox_model_full$y)]
data$transformed_Amyloid_PET_full <-
  (data$Amyloid_PET ^ best_lambda_full - 1) / best_lambda_full
transformed_model_full <-
  lm(
    transformed_Amyloid_PET_full ~ ptau181 + pTau217 + GFAP + NFL + Ab_ratio
    + sex + Calculated_Consensus_dx + apoe_ternary,
    data = data
  )

# Reduced Model BoxCox
boxcox_model_reduced <-
  boxcox(lm(Amyloid_PET ~ pTau217 + GFAP + NFL + Ab_ratio + sex,
            data = data))
best_lambda_reduced <-
  boxcox_model_reduced$x[which.max(boxcox_model_reduced$y)]
data$transformed_Amyloid_PET_reduced <-
  (data$Amyloid_PET ^ best_lambda_reduced - 1) / best_lambda_reduced
transformed_model_reduced <-
  lm(
    transformed_Amyloid_PET_reduced ~ pTau217 + GFAP + NFL + Ab_ratio + sex,
    data = data
  )

plot_model(transformed_model_full)
plot_model(transformed_model_reduced)




# Calculate Cook's distances
cooksd_full <- cooks.distance(log_model_full)
cooksd_reduced <- cooks.distance(log_model_reduced)

# Identify influential points
influential_points_full <-
  as.numeric(names(cooksd_full)[(cooksd_full > (4 / nrow(data)))])
influential_points_reduced <-
  as.numeric(names(cooksd_reduced)[(cooksd_reduced > (4 / nrow(data)))])

inf_pts_full <- 
  which(rownames(data) %in% as.character(influential_points_full))
inf_pts_reduced <- 
  which(rownames(data) %in% as.character(influential_points_reduced))


str(inf_pts_full)
str(inf_pts_reduced)

# Remove influential points
data_clean_full <- data[-inf_pts_full, ]
data_clean_reduced <- data[-inf_pts_reduced, ]

str(data_clean_full)
str(data_clean_reduced)

# Re-fit the models
# Full
clean_model_full <-
  lm(
    log_Amyloid_PET ~ ptau181 + pTau217 + GFAP + NFL + Ab_ratio + sex
    + Calculated_Consensus_dx + apoe_ternary,
    data = data_clean_full
  )
# Reduced
clean_model_reduced <-
  lm(
    log_Amyloid_PET ~ pTau217 + GFAP + NFL + Ab_ratio + sex,
    data = data_clean_reduced
  )

plot_model(clean_model_full)
plot_model(clean_model_reduced)


################################################################################
# Re-Define dataset and predictors
data_full <- data_clean_full
data_reduced <- data_clean_reduced

# Response Variable
response_var <- "log_Amyloid_PET"
# Ensure the response variable is numeric
data_full[[response_var]] <- 
  as.numeric(as.character(data_full[[response_var]]))
data_reduced[[response_var]] <- 
  as.numeric(as.character(data_reduced[[response_var]]))


# Combinations of Predictors found with PRESS
best_formula_full <- 
  "log_Amyloid_PET ~ ptau181 + pTau217 + GFAP + NFL + Ab_ratio + sex 
                 + Calculated_Consensus_dx + apoe_ternary"
best_formula_reduced <- 
  "log_Amyloid_PET ~ pTau217 + GFAP + NFL + Ab_ratio + sex"



################################################################################
##########################  FUNCTION DEFINITIONS  ##############################
################################################################################

# Function to perform LOOCV and gather metrics for the best model
perform_loocv_and_gather_metrics <-
  function(data, response_var, formula_str) {
    # Ensure the data is a data frame
    if (!is.data.frame(data)) {
      stop("The data is not a data frame.")
    }
    
    # Print the number of rows to debug
    print(paste("Number of rows in data:", nrow(data)))
    
    # Print the formula to debug
    print(paste("Using formula:", formula_str))
    
    # Initialize vectors to store predictions and actual values
    predictions <- numeric(nrow(data))
    true_labels <- numeric(nrow(data))
    
    # Ensure the response variable is numeric
    data[[response_var]] <- as.numeric(data[[response_var]])
    
    # Perform LOOCV
    for (i in 1:nrow(data)) {
      train_data <- data[-i, ]
      test_data <- data[i, , drop = FALSE]
    
      fit <-
        eval(parse(text = paste("lm(", formula_str, ", data = train_data)")))
      
      prediction <- predict(fit, test_data)
      
      predictions[i] <- prediction
      true_labels[i] <- test_data[[response_var]]
    }
    
    # Calculate and return metrics
    mse <- mean((predictions - true_labels) ^ 2)
    rmse <- sqrt(mse)
    mae <- mean(abs(predictions - true_labels))
    r_squared <- summary(lm(predictions ~ true_labels))$r.squared
    adj_r_squared <-
      summary(lm(predictions ~ true_labels))$adj.r.squared
    
    return(
      list(
        mse = mse,
        rmse = rmse,
        mae = mae,
        r_squared = r_squared,
        adj_r_squared = adj_r_squared,
        predictions = predictions,
        actual_values = true_labels
      )
    )
  }


################################################################################

# Function to plot variable importance
plot_variable_importance <- function(fit) {
  coef_df <- as.data.frame(coef(summary(fit)))
  coef_df$Variable <- rownames(coef_df)
  coef_df <-
    coef_df[order(abs(coef_df$Estimate), decreasing = TRUE), ]
  
  ggplot(coef_df, aes(x = reorder(Variable, abs(Estimate)), 
                      y = Estimate)) +
    geom_bar(stat = "identity", 
             fill = "#97C835") +
    coord_flip() +
    labs(title = "Variable Importance", 
         x = "Variables", 
         y = "Coefficient Estimate") +
    theme_minimal()
}

################################################################################

# Function to plot prediction errors
plot_prediction_errors <- function(predictions, true_labels) {
  error_df <-
    data.frame(Predictions = predictions, TrueValues = true_labels)
  ggplot(error_df, aes(x = TrueValues, y = Predictions)) +
    geom_point(color = "#97C835") +
    geom_abline(linetype = "dashed", color = "gray") +
    labs(title = "Prediction Errors", x = "True Values", y = "Predictions") +
    theme_minimal()
}

################################################################################

# Function to analyze predictions above and below threshold
analyze_threshold <-
  function(predictions_full,
           true_labels_full,
           predictions_reduced,
           true_labels_reduced,
           threshold = log(1.17)) {
    # Initialize lists to store results
    results <- list()
    
    # Define function to calculate metrics
    calculate_metrics <-
      function(predictions, true_labels) {
        mse <- mean((predictions - true_labels) ^ 2)
        rmse <- sqrt(mse)
        mae <- mean(abs(predictions - true_labels))
        return(list(
          mse = mse,
          rmse = rmse,
          mae = mae
        ))
      }
    
    # Calculate metrics for below and above threshold for full model
    below_threshold_full <-
      calculate_metrics(predictions_full[true_labels_full < threshold],
                        true_labels_full[true_labels_full < threshold])
    above_threshold_full <-
      calculate_metrics(predictions_full[true_labels_full >= threshold],
                        true_labels_full[true_labels_full >= threshold])
    
    # Calculate metrics for below and above threshold for reduced model
    below_threshold_reduced <-
      calculate_metrics(predictions_reduced[true_labels_reduced < threshold],
                        true_labels_reduced[true_labels_reduced < threshold])
    above_threshold_reduced <-
      calculate_metrics(predictions_reduced[true_labels_reduced >= threshold],
                        true_labels_reduced[true_labels_reduced >= threshold])
    
    # Print metrics for verification
    print("Full Model Below Threshold Metrics:")
    print(below_threshold_full)
    print("Full Model Above Threshold Metrics:")
    print(above_threshold_full)
    print("Reduced Model Below Threshold Metrics:")
    print(below_threshold_reduced)
    print("Reduced Model Above Threshold Metrics:")
    print(above_threshold_reduced)
    
    # Combine results into a data frame
    results_df <- data.frame(
      Metric = rep(c("MSE", "RMSE", "MAE"), 4),
      Value = c(
        below_threshold_full$mse,
        below_threshold_full$rmse,
        below_threshold_full$mae,
        above_threshold_full$mse,
        above_threshold_full$rmse,
        above_threshold_full$mae,
        below_threshold_reduced$mse,
        below_threshold_reduced$rmse,
        below_threshold_reduced$mae,
        above_threshold_reduced$mse,
        above_threshold_reduced$rmse,
        above_threshold_reduced$mae
      ),
      Model = rep(c("Full", "Reduced"), each = 6),
      Threshold = rep(
        c("Below Threshold", "Above Threshold"),
        each = 3,
        times = 2
      )
    )
    
    return(results_df)
  }


# Function to create and display a summary table
create_summary_table1 <-
  function(threshold_analysis_df) {
    kable(
      threshold_analysis_df,
      format = "html",
      table.attr = "class='table table-bordered'",
      caption = "Threshold Analysis of Model Performance Metrics"
    ) %>%
      kable_styling(full_width = FALSE, position = "center")
  }


# Function to plot threshold analysis results
plot_threshold_analysis1 <-
  function(threshold_analysis_df) {
    ggplot(threshold_analysis_df,
           aes(
             x = interaction(Threshold, Metric),
             y = Value,
             fill = Model
           )) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(values = c("Full" = "#4B611F",
                                   "Reduced" = "tomato")) +
      labs(title = "Threshold Analysis of Model Performance",
           x = "Metric (Threshold)",
           y = "Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

# Function to plot threshold analysis results with facets
plot_threshold_analysis2 <- 
  function(threshold_analysis_df) {
    ggplot(threshold_analysis_df, aes(x = Metric, 
                                      y = Value,
                                      fill = Model)) +
      geom_bar(stat = "identity", 
               position = position_dodge()) +
      facet_wrap(~ Threshold, scales = "free_x") +
      scale_fill_manual(values = c("Full" = "#4B611F", 
                                   "Reduced" = "tomato")) +
      labs(title = "Threshold Analysis of Model Performance", 
           x = "Metric", 
           y = "Value") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }


################################################################################

# Function to plot overall model performance
plot_overall_performance <-
  function(metrics) {
    performance_df <- data.frame(
      Metric = c("MSE", "RMSE", "MAE", "R-squared", "Adjusted R-squared"),
      Value = c(
        metrics$mse,
        metrics$rmse,
        metrics$mae,
        metrics$r_squared,
        metrics$adj_r_squared
      )
    )
    ggplot(performance_df, aes(x = Metric, y = Value)) +
      geom_bar(stat = "identity",
               fill = "#97C835") +
      labs(title = "Overall Model Performance",
           x = "Metric",
           y = "Value") +
      theme_minimal()
  }


################################################################################

# Function to plot side-by-side comparison of metrics
plot_comparison <-
  function(metrics_full, metrics_reduced) {
    # Creating a dataframe for comparison
    comparison_df <- data.frame(
      Metric = c("MSE", "RMSE", "MAE", "R-squared", "Adjusted R-squared"),
      Full = c(
        metrics_full$mse,
        metrics_full$rmse,
        metrics_full$mae,
        metrics_full$r_squared,
        metrics_full$adj_r_squared
      ),
      Reduced = c(
        metrics_reduced$mse,
        metrics_reduced$rmse,
        metrics_reduced$mae,
        metrics_reduced$r_squared,
        metrics_reduced$adj_r_squared
      )
    )
    
    # Melting the dataframe for ggplot
    comparison_df_long <-
      tidyr::gather(comparison_df,
                    key = "Model",
                    value = "Value",
                    Full,
                    Reduced)
    
    # Plotting the comparison
    ggplot(comparison_df_long, aes(x = Metric,
                                   y = Value,
                                   fill = Model)) +
      geom_bar(stat = "identity",
               position = position_dodge()) +
      scale_fill_manual(values = c("Full" = "#4B611F",
                                   "Reduced" = "tomato")) +
      labs(title = "Comparison of Model Performance Metrics",
           x = "Metric",
           y = "Value") +
      theme_minimal()
  }


# Function to create a summary table
create_summary_table2 <-
  function(metrics_full, metrics_reduced) {
    summary_df <- data.frame(
      Metric = c("MSE", "RMSE", "MAE", "R-squared", "Adjusted R-squared"),
      Full = c(
        metrics_full$mse,
        metrics_full$rmse,
        metrics_full$mae,
        metrics_full$r_squared,
        metrics_full$adj_r_squared
      ),
      Reduced = c(
        metrics_reduced$mse,
        metrics_reduced$rmse,
        metrics_reduced$mae,
        metrics_reduced$r_squared,
        metrics_reduced$adj_r_squared
      )
    )
    return(summary_df)
  }

################################################################################
################################################################################

# Fit the model using all data for variable importance plot
fit_full <- 
  eval(parse(text = paste("lm(", best_formula_full, ", 
                          data = data_clean_full)")))
fit_reduced <- 
  eval(parse(text = paste("lm(", best_formula_reduced, ", 
                          data = data_clean_reduced)")))


# Plot variable importance for the full model
plot_variable_importance(fit_full)
# Plot variable importance for the reduced model
plot_variable_importance(fit_reduced)

# Calculate AIC and BIC to compare 2 models
cat("AIC for Full Model:", AIC(clean_model_full), "\n")
cat("AIC for Reduced Model:", AIC(clean_model_reduced), "\n")

cat("BIC for Full Model:", BIC(clean_model_full), "\n")
cat("BIC for Reduced Model:", BIC(clean_model_reduced), "\n")

# Perform LOOCV and gather metrics for the full model
metrics_full <-
  perform_loocv_and_gather_metrics(data_clean_full, 
                                   response_var, 
                                   best_formula_full)

# Perform LOOCV and gather metrics for the reduced model
metrics_reduced <-
  perform_loocv_and_gather_metrics(data_clean_reduced, 
                                   response_var, 
                                   best_formula_reduced)

# Print metrics for both models
print("Metrics for the best linear regression model (Full):")
print(metrics_full)

print("Metrics for the best linear regression model (Reduced):")
print(metrics_reduced)

# Plot overall model performance for the full model
plot_overall_performance(metrics_full)

# Plot overall model performance for the reduced model
plot_overall_performance(metrics_reduced)

# Plot prediction errors for the full model
plot_prediction_errors(metrics_full$predictions,
                       metrics_full$actual_values)

# Plot prediction errors for the reduced model
plot_prediction_errors(metrics_reduced$predictions, 
                       metrics_reduced$actual_values)


# Analyze predictions above and below the threshold for both models
threshold_analysis_df <- 
  analyze_threshold(
  metrics_full$predictions, 
  metrics_full$actual_values,
  metrics_reduced$predictions, 
  metrics_reduced$actual_values
)

# Reorder the dataframe by Model and Threshold
threshold_analysis_df <- threshold_analysis_df %>%
  arrange(Model, Threshold)

create_summary_table1(threshold_analysis_df)
# 2 options for how you want to display
plot_threshold_analysis1(threshold_analysis_df)
plot_threshold_analysis2(threshold_analysis_df)

# Plot comparison of metrics for full and reduced models
plot_comparison(metrics_full, metrics_reduced)


# Create and print summary table of metrics
summary_table <-
  create_summary_table2(metrics_full, metrics_reduced)
print("Summary Table of Metrics:")
print(summary_table)


# Align the datasets for Wilcox test 
# (different lengths due to removing different influential points)


# Verify the lengths of the predictions
length_full_predictions <- length(metrics_full$predictions)
length_reduced_predictions <- length(metrics_reduced$predictions)

cat("Length of Full Model Predictions:",
    length_full_predictions,
    "\n")
cat("Length of Reduced Model Predictions:",
    length_reduced_predictions,
    "\n")

# If lengths are unequal, align the data sets
if (length_full_predictions != length_reduced_predictions) {
  cat("The lengths of predictions are not equal. Aligning data sets...\n")
  
  # Get the indices of the overlapping rows based on the unique IDs
  overlapping_indices <- intersect(rownames(data_clean_full),
                                   rownames(data_clean_reduced))
  
  # Filter both datasets to include only overlapping rows
  data_clean_full <- data_clean_full[overlapping_indices,]
  data_clean_reduced <- data_clean_reduced[overlapping_indices,]
  
  # Recompute the predictions
  metrics_full <-
    perform_loocv_and_gather_metrics(data_clean_full,
                                     response_var,
                                     best_formula_full)
  metrics_reduced <-
    perform_loocv_and_gather_metrics(data_clean_reduced,
                                     response_var,
                                     best_formula_reduced)
}

# Verify the lengths again
length_full_predictions <- length(metrics_full$predictions)
length_reduced_predictions <- length(metrics_reduced$predictions)

cat("Length of Full Model Predictions after alignment:",
    length_full_predictions,
    "\n")
cat(
  "Length of Reduced Model Predictions after alignment:",
  length_reduced_predictions,
  "\n"
)

# Perform Wilcoxon signed-rank test for MSE if lengths are equal
if (length_full_predictions == length_reduced_predictions) {
  wilcox_test_mse <-
    wilcox.test(metrics_full$predictions,
                metrics_reduced$predictions,
                paired = TRUE)
  
  print(wilcox_test_mse)
} else {
  cat(
    "The lengths of predictions are still not equal. 
    Cannot perform Wilcoxon signed-rank test.\n"
  )
}

