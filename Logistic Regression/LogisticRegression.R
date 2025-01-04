##################################
##########   STAT 409   ##########
####  Data Science Practicum  ####
##################################
###### LOGISTIC REGRESSION #######
##################################


##  NOTE: Run 'Data_Cleaning_and_Integration.R' first


# Load necessary libraries
library(dplyr)
library(caret)
library(ROCR)
library(ggplot2)
library(scales)

# Function to perform LOOCV and gather metrics for the best model
perform_loocv_and_gather_metrics <-
  function(data, response_var, formula_str) {
    
    # Print the number of rows to debug
    print(paste("Number of rows in data:", nrow(data)))
    
    # Initialize vectors to store predictions and actual values
    predictions <- numeric(nrow(data))
    pred_probs <- numeric(nrow(data))
    true_labels <- numeric(nrow(data))
    confusion_matrices <- list()
    
    # Ensure the response variable is numeric
    data[[response_var]] <- as.numeric(data[[response_var]])
    
    # Perform LOOCV
    for (i in 1:nrow(data)) {
      train_data <- data[-i, ]
      test_data <- data[i, , drop = FALSE]
      
      fit <-
        eval(parse(
          text = paste("glm(", formula_str, ", 
                       data = train_data, family = binomial)")
        ))
      prediction <- predict(fit, test_data, type = "response")
      predicted_class <- ifelse(prediction > 0.5, 1, 0)
      
      predictions[i] <- predicted_class
      pred_probs[i] <- prediction
      true_labels[i] <- test_data[[response_var]]
      
      # Print predictions and probabilities for debugging
      print(paste("Prediction for row", i, ":", prediction))
      
      # Create a 2x2 confusion matrix
      confusion_matrix <-
        matrix(
          0,
          nrow = 2,
          ncol = 2,
          dimnames = list(
            "Predicted" = c("0", "1"),
            "Actual" = c("0", "1")
          )
        )
      confusion_matrix[as.character(predicted_class), 
                       as.character(test_data[[response_var]])] <-  1
      confusion_matrices[[i]] <- confusion_matrix
    }
    
    # Calculate and return metrics
    accuracy <- mean(predictions == true_labels)
    sensitivity <-
      sum(predictions == 1 & true_labels == 1) / sum(true_labels == 1)
    specificity <-
      sum(predictions == 0 & true_labels == 0) / sum(true_labels == 0)
    precision <-
      sum(predictions == 1 & true_labels == 1) / sum(predictions == 1)
    f1_score <-
      2 * (precision * sensitivity) / (precision + sensitivity)
    
    # Calculate AUC and ROC
    pred <- prediction(pred_probs, true_labels)
    perf <- performance(pred, measure = "tpr", "fpr")
    auc <- performance(pred, measure = "auc")@y.values[[1]]
    
    # Calculate Log-Loss
    logloss <-
      -(mean(
        true_labels * log(pred_probs) + (1 - true_labels) * log(1 - pred_probs)
      ))
    
    # Calculate Brier Score
    brier <- mean((pred_probs - true_labels) ^ 2)
    
    aggregated_confusion_matrix <- Reduce("+", confusion_matrices)
   
    return(
      list(
        accuracy = accuracy,
        sensitivity = sensitivity,
        specificity = specificity,
        precision = precision,
        f1_score = f1_score,
        auc = auc,
        logloss = logloss,
        brier = brier,
        predictions = predictions,
        pred_probs = pred_probs,
        actual_values = true_labels,
        aggregated_confusion_matrix = aggregated_confusion_matrix,
        roc_curve = perf
      )
    )
  }

# Plot aggregated ROC curve
plot_aggregated_roc_curve <- function(roc_curve, auc) {
  # Plot the ROC curve
  plot(roc_curve, col = "red", lwd = 2, main = "Aggregated ROC Curve")
  # Diagonal reference line
  abline(a = 0, b = 1, col = "#97C835")
  # Add AUC value to the plot
  legend("bottomright", 
         legend = sprintf("AUC = %.3f", auc), 
         col = "red", lwd = 2)
}

################################################################################
################################################################################

# Define dataset and predictors
data <- data_time_merge

# Response Variable
response_var <- "amyloid_binary"
# Ensure the response variable is numeric (0/1)
data[[response_var]] <- as.numeric(as.character(data[[response_var]]))

# Combinations of Predictors found with LogLoss
# LogLoss:  0.157409889
best_formula_full <- 
  "amyloid_binary ~ ptau231 + ptau181 + pTau217 + GFAP + NFL + Ab_ratio 
                    + age_at_appointment + apoe_ternary"

# LogLoss:  0.1692471
best_formula_reduced <- 
  "amyloid_binary ~ ptau181 + pTau217 + GFAP + NFL + Ab_ratio"

################################################################################


# Function to plot variable importance
plot_variable_importance <- function(fit) {
  library(ggplot2)
  
  coef_df <- as.data.frame(coef(summary(fit)))
  coef_df$Variable <- rownames(coef_df)
  coef_df <-
    coef_df[order(abs(coef_df$Estimate), decreasing = TRUE),]
  
  ggplot(coef_df, aes(x = reorder(Variable, abs(Estimate)), 
                      y = Estimate)) +
    geom_bar(stat = "identity", fill = "#97C835") +
    coord_flip() +
    labs(title = "Variable Coefficients", 
         x = "Variables", 
         y = "Coefficient Estimate") +
    theme_minimal()
}

# Fit the model using all data for variable importance plot
fit_full <-
  eval(parse(
    text = paste("glm(", best_formula_full, ", 
                       data = data, family = binomial)")
  ))
plot_variable_importance(fit_full)

# Fit the model using all data for variable importance plot
fit_reduced <-
  eval(parse(
    text = paste("glm(", best_formula_reduced, ", 
                       data = data, family = binomial)")
  ))
plot_variable_importance(fit_reduced)

################################################################################

anova(fit_reduced, fit_full, test = "Chisq")

################################################################################

# Perform LOOCV and gather metrics for the best model considering ALL predictors
metrics_full <- 
  perform_loocv_and_gather_metrics(data, response_var, best_formula_full)

# Print metrics
print("Metrics for the best logistic regression model (Full):")
print(metrics_full)

# Print aggregated confusion matrix
print("Aggregated Confusion Matrix (Full):")
print(metrics_full$aggregated_confusion_matrix)

# Plot the aggregated ROC curve
plot_aggregated_roc_curve(metrics_full$roc_curve, metrics_full$auc)

################################################################################

# Perform LOOCV and gather metrics for the best model w/ only reduced predictors
metrics_reduced <- 
  perform_loocv_and_gather_metrics(data, response_var, best_formula_reduced)

# Print metrics
print("Metrics for the best logistic regression model (Reduced):")
print(metrics_reduced)

# Print aggregated confusion matrix
print("Aggregated Confusion Matrix (Reduced):")
print(metrics_reduced$aggregated_confusion_matrix)

# Plot the aggregated ROC curve
plot_aggregated_roc_curve(metrics_reduced$roc_curve, metrics_reduced$auc)

################################################################################

# Function to plot histogram of prediction probabilities
plot_prediction_probabilities <- function(pred_probs, title) {
  # Print the type and first few values of pred_probs for debugging
  print(paste("Type of pred_probs:", typeof(pred_probs)))
  print(paste("First few values of pred_probs:", head(pred_probs)))

  # Plot the histogram
  hist(
    pred_probs,
    breaks = 30,
    main = title,
    xlab = "Prediction Probabilities",
    col = "#97C835"
  )
}

# Plot prediction probabilities for the full model
plot_prediction_probabilities(metrics_full$pred_probs, 
                              "Prediction Probabilities (Full)")

# Plot prediction probabilities for the reduced model
plot_prediction_probabilities(metrics_reduced$pred_probs, 
                              "Prediction Probabilities (Reduced)")

# Compare distributions of prediction probabilities
boxplot(as.numeric(metrics_full$pred_probs), 
        as.numeric(metrics_reduced$pred_probs), 
        names = c("Full", "Reduced"), 
        main = "Prediction Probabilities Comparison", 
        col = c("#4B611F", "tomato"))


################################################################################
##########################  COMPARING THE TWO MODELS  ##########################
################################################################################

# Function to plot ROC curves for both models
plot_roc_curves <-
  function(roc_curve_full,
           auc_full,
           roc_curve_reduced,
           auc_reduced) {
    plot(
      roc_curve_full,
      col = "#4B611F",
      lwd = 2,
      main = "ROC Curves Comparison"
    )
    plot(roc_curve_reduced,
         col = "tomato",
         lwd = 2,
         add = TRUE)
    legend(
      "bottomright",
      legend = c(paste(
        "Full AUC =", round(auc_full, 3)
      ),
      paste("Reduced AUC =", round(auc_reduced, 3))),
      col = c("#4B611F", "tomato"),
      lwd = 2
    )
    abline(a = 0, b = 1, col = "grey")
  }

# Plot ROC curves
plot_roc_curves(
  metrics_full$roc_curve,
  metrics_full$auc,
  metrics_reduced$roc_curve,
  metrics_reduced$auc
)

################################################################################

# Function to plot Precision-Recall curves for both models
plot_precision_recall_curves <-
  function(pred_probs_full,
           true_labels_full,
           pred_probs_reduced,
           true_labels_reduced) {
    # Full model
    pred_full <-
      prediction(pred_probs_full, true_labels_full)
    perf_full <- performance(pred_full, "prec", "rec")
    
    # Reduced model
    pred_reduced <- prediction(pred_probs_reduced, true_labels_reduced)
    perf_reduced <- performance(pred_reduced, "prec", "rec")
    
    # Plot the curves
    plot(perf_full,
         col = "#4B611F",
         lwd = 2,
         main = "Precision-Recall Curves Comparison")
    plot(perf_reduced,
         col = "tomato",
         lwd = 2,
         add = TRUE)
    legend(
      "bottomleft",
      legend = c("Full", "Reduced"),
      col = c("#4B611F", "tomato"),
      lwd = 2
    )
  }

# Plot Precision-Recall curves
plot_precision_recall_curves(
  metrics_full$pred_probs,
  metrics_full$actual_values,
  metrics_reduced$pred_probs,
  metrics_reduced$actual_values
)

################################################################################

# Function to plot calibration curves
plot_calibration_curves <- function(pred_probs_full, true_labels_full,
                                    pred_probs_reduced, true_labels_reduced) {
  
  # Convert factors to characters for compatibility
  true_labels_full <- as.factor(true_labels_full)
  true_labels_reduced <- as.factor(true_labels_reduced)
  
  # Calibration for full model
  cal_full <- calibration(
    true_labels_full ~ (1 - pred_probs_full),
    data = data.frame(pred_probs_full, true_labels_full)
  )
  cal_full_df <- as.data.frame(cal_full$data)
  
  # Calibration for reduced model
  cal_reduced <- calibration(
    true_labels_reduced ~ (1 - pred_probs_reduced),
    data = data.frame(pred_probs_reduced, true_labels_reduced)
  )
  cal_reduced_df <- as.data.frame(cal_reduced$data)
  
  # Plotting the calibration curves using ggplot2
  ggplot() +
    geom_abline(intercept = 0, 
                slope = 1, 
                color = "darkgray",
                linewidth = 0.5) +
    geom_line(
      data = cal_full_df,
      aes(x = midpoint/100, y = Percent / 100, color = "Full"),
      linewidth = 1.2
    ) +
    geom_line(
      data = cal_reduced_df,
      aes(x = midpoint/100, y = Percent / 100, color = "Reduced"),
      linewidth = 1.2
    ) +
    geom_smooth(
      data = cal_full_df,
      aes(x = midpoint/100, y = Percent / 100, color = "Full"),
      method = "loess",
      se = FALSE,
      linewidth = 1,
      linetype = "dotted"
    ) +
    geom_smooth(
      data = cal_reduced_df,
      aes(x = midpoint/100, y = Percent / 100, color = "Reduced"),
      method = "loess",
      se = FALSE,
      linewidth = 1,
      linetype = "dotted"
    ) +
    scale_color_manual(
      name = "Model",
      values = c("Full" = "#4B611F", "Reduced" = "tomato")
    ) +
    theme_minimal() +
    labs(x = "Predicted Probability", y = "Observed Probability", 
         title = "Calibration Curves") +
    theme(legend.position = "bottom") +
    coord_cartesian(ylim = c(0, 1))
}

# Example usage
plot_calibration_curves(
  metrics_full$pred_probs,
  metrics_full$actual_values,
  metrics_reduced$pred_probs,
  metrics_reduced$actual_values
)

################################################################################

# Function to plot metrics at different thresholds
plot_threshold_metrics <- function(pred_probs, true_labels) {
  library(ggplot2)
  library(caret)
  
  thresholds <- seq(0, 1, by = 0.01)
  metrics <- data.frame(
    Threshold = thresholds,
    Accuracy = numeric(length(thresholds)),
    Precision = numeric(length(thresholds)),
    Recall = numeric(length(thresholds)),
    F1 = numeric(length(thresholds))
  )
  
  true_labels <- factor(true_labels, levels = c(0, 1))
  
  for (i in seq_along(thresholds)) {
    threshold <- thresholds[i]
    predicted <-
      factor(ifelse(pred_probs >= threshold, 1, 0), levels = c(0, 1))
    
    confusion <- confusionMatrix(predicted, true_labels)
    metrics$Accuracy[i] <- confusion$overall['Accuracy']
    metrics$Precision[i] <- confusion$byClass['Pos Pred Value']
    metrics$Recall[i] <- confusion$byClass['Sensitivity']
    metrics$F1[i] <- confusion$byClass['F1']
  }
  
  metrics <- metrics %>%
    tidyr::drop_na()
  
  ggplot(metrics, aes(x = Threshold)) +
    geom_line(aes(y = Accuracy, color = "Accuracy"), linewidth = 1) +
    geom_line(aes(y = Precision, color = "Precision"), linewidth = 1) +
    geom_line(aes(y = Recall, color = "Recall"), linewidth = 1) +
    geom_line(aes(y = F1, color = "F1"), linewidth = 1) +
    labs(title = "Metrics at Different Thresholds", y = "Metric Value") +
    scale_color_manual(
      name = "Metric",
      values = c(
        "Accuracy" = "blue",
        "Precision" = "green",
        "Recall" = "red",
        "F1" = "purple"
      )
    ) +
    theme_minimal()
}

# Plot threshold metrics for the full model
plot_threshold_metrics(metrics_full$pred_probs,
                       metrics_full$actual_values)

# Plot threshold metrics for the reduced model
plot_threshold_metrics(metrics_reduced$pred_probs, 
                       metrics_reduced$actual_values)


################################################################################

