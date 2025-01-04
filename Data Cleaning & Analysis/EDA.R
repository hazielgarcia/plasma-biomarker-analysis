##################################
##########   STAT 409   ##########
####  Data Science Practicum  ####
##################################
###########################################################
###########  EXPLORATORY DATA ANALYSIS (EDA)  #############
###########################################################

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(corrplot)
library(patchwork)
library(cowplot)


# View the first few rows of the dataset
head(data_time_merge)
# Summary statistics
summary(data_time_merge)
# Structure of the dataset
str(data_time_merge)
# Names of dataset
names(data_time_merge)
# View entire dataset in other window
#View(data_time_merge)


## QUALITATIVE DATA
# Proportion of amyloid_binary in dataset
prop.table(table(data_time_merge$amyloid_binary))

categorical_vars <- 
  names(data_time_merge) == "amyloid_binary" |
  names(data_time_merge) == "apoe_ternary" |
  names(data_time_merge) == "sex" |
  names(data_time_merge) == "race_primary" |
  names(data_time_merge) == "hispanic" |
  names(data_time_merge) == "Calculated_Consensus_dx"

cat_data <- data_time_merge[, categorical_vars]


# CALCULATE RELATIVE FREQUENCIES
rel_freq_list <- lapply(names(cat_data), function(var) {
  data_to_count <- cat_data %>%
    dplyr::select(dplyr::all_of(var), amyloid_binary) %>%
    dplyr::count(!!rlang::sym(var), amyloid_binary) %>%
    dplyr::group_by(amyloid_binary) %>%
    dplyr::mutate(Total = sum(n)) %>%
    dplyr::ungroup() %>%
    # Compute the relative frequencies
    dplyr::mutate(Rel_Freq = n / Total * 100) %>%
    # Reshape the data using pivot_wider
    tidyr::pivot_wider(
      names_from = amyloid_binary,
      values_from = c(n, Rel_Freq),
      names_glue = "{.value}_{amyloid_binary}"
    ) %>%
    dplyr::mutate(Variable = var)
}) %>%
  dplyr::bind_rows() %>%
  dplyr::select(Variable, everything())

# View the frequency list in a pop-up window
#View(rel_freq_list)

# CALCULATE TOTAL PERCENTAGES FOR CATEGORIES OF EACH VARIABLE
# Initialize empty list to store category percentage results
category_percentages_list <- list()

# Loop through each categorical variable to calculate percentages
for (var in names(cat_data)) {
  # Calculate counts and percentages
  percentages_df <- cat_data %>%
    dplyr::select(dplyr::all_of(var)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(var))) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = 'drop') %>%
    dplyr::mutate(Percentage = Count / sum(Count) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Variable = var,
                  Category = dplyr::across(dplyr::all_of(var))) %>%
    dplyr::select(Variable, Count, Percentage, Category)
  
  # Store the dataframe in the list
  category_percentages_list[[var]] <- percentages_df
}

# Combine all dataframes into one
all_category_percentages <-
  do.call(dplyr::bind_rows, category_percentages_list[])

# View the combined results
#View(all_category_percentages)


## QUALITATIVE
### Graphs of Qualitative Variables

# List to hold plots
plots_list <- list()
str(data_time_merge)

# Loop through categorical variables and create ggplots
for (var in names(data_time_merge[categorical_vars])) {
  # For the first variable, keep the legend
  if (var == names(data_time_merge[categorical_vars])[1]) {
    p <-
      ggplot(data_time_merge, aes_string(x = var, fill = "amyloid_binary")) +
      geom_bar(position = position_stack(reverse = TRUE)) + # Reverse stack order
      scale_fill_manual(labels = c("0", "1"),
                        values = c("#DC310F", "#97C835")) +
      labs(x = var, y = "Count") +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      theme_minimal()
  } else {
    # For all other variables, remove the legend
    p <-
      ggplot(data_time_merge, aes_string(x = var, fill = "amyloid_binary")) +
      geom_bar(position = position_stack(reverse = TRUE)) + # Reverse stack order
      scale_fill_manual(labels = c("Negative", "Positive"),
                        values = c("#DC310F", "#97C835")) +
      labs(x = var, y = "Count") +
      theme_minimal() +
      theme(legend.position = "none") +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
  }
  plot(p)
  # Add the plot to the list
  plots_list[[var]] <- p
}

# Create a blank plot to balance visualization
#blank_plot <- ggplot() +
# theme_void() +
#  theme(
#   plot.background = element_blank(),
#  panel.grid = element_blank(),
# panel.background = element_blank()
#)

# Position to insert the blank plot
#pos_insert_blank <- length(plots_list) - 2

# Insert the blank plot into the list at the specific position
#plots_list <- append(plots_list, list(blank_plot), pos_insert_blank)

# Arrange the plots on a grid, with only the first plot having a legend
grid.arrange(grobs = plots_list, ncol = 3)


## QUANTITATIVE
# Set up variables for numeric values
numeric_data <- 
  names(data_time_merge) == "Amyloid_PET" |
  names(data_time_merge) == "age_at_appointment" |
  names(data_time_merge) == "time_diff" |
  names(data_time_merge) == "ptau231" |
  names(data_time_merge) == "ptau181" |
  names(data_time_merge) == "pTau217" |
  names(data_time_merge) == "GFAP" |
  names(data_time_merge) == "NFL" |
  names(data_time_merge) == "Ab_ratio" 

continuous_vars <- c("Amyloid_PET", "age_at_appointment", "time_diff", 
                     "ptau231", "ptau181", "pTau217", 
                     "GFAP", "NFL", "Ab_ratio" )

### GRAPHS - Continuous VARIABLES
# Clear plot_list
plots_list <- list()
# Hist & boxplots of continuous
for (var in continuous_vars) {
  # Histogram for the variable
  p_hist <- ggplot(data_time_merge, aes_string(x = var)) +
    geom_histogram(
      aes(y = ..density..),
      position = "identity",
      #alpha = 0.5,
      binwidth = diff(range(data_time_merge[[var]], na.rm = TRUE)) / 30,
      fill = "#97C835",
      color = "black"
    ) +
    geom_density(alpha = .5,
                 fill = "#FFCCCC",
                 size = 0.5) +
    labs(title = paste("Histogram of", var),
         x = var,
         y = "Density") +
    theme_minimal()
  
  p_box <-
    ggplot(data_time_merge, aes_string(y = var, fill = "amyloid_binary")) +
    geom_boxplot(aes(x = factor(amyloid_binary))) +
    scale_fill_manual(values = c("#DC310F", "#97C835")) +
    labs(
      title = paste("Boxplot of", var, "by Amyloid Binary Values"),
      x = "",
      y = var
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  # Combine the plots side by side for the current variable
  combined_plot <- arrangeGrob(p_hist, p_box, ncol = 2)
  
  # Add the combined plot to the list
  plots_list[[var]] <- combined_plot
}

# Display all combined plots together, each pair in a row
do.call(grid.arrange, c(plots_list[1:3], ncol = 1))
do.call(grid.arrange, c(plots_list[4:6], ncol = 1))
do.call(grid.arrange, c(plots_list[7:9], ncol = 1))

# Creating a Correlation Matrix 
numeric_corr <- cor(data_time_merge[numeric_data])

# Only keep the upper triangle (including the diagonal)
numeric_corr_upper <- upper.tri(numeric_corr)

# Replace the lower triangle and diagonal with NA
numeric_corr[!numeric_corr_upper] <- NA

# PLOTS
# Visualize the correlation matrix
corrplot(
  numeric_corr,
  method = "circle",
  type = "upper",
  tl.col = "black",
  tl.srt = 45,
  addCoef.col = "black",
  diag = FALSE,
  is.corr = TRUE
)



##  DATA VISUALIZATION
# Create scatter plots for each biomarker vs Amyloid PET
# List of biomarkers
biomarkers <- c("Ab_ratio","GFAP","NFL","ptau181", "pTau217", "ptau231") 
plot_list <- list()

for (biomarker in biomarkers) {
  p <- ggplot(data_time_merge, aes_string(x = biomarker, y = "Amyloid_PET")) +
    geom_point(color = "blue") +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  
    labs(title = paste("Scatter plot of", biomarker, "vs Amyloid PET"),
         x = biomarker, y = "Amyloid PET")
  
  plot_list[[biomarker]] <- p
}

# Combine scatter plots into a grid
grid.arrange(grobs = plot_list, ncol = 3)


# Checking for any trends with age
amyloid_v_age <- ggplot(data_time_merge, aes(x = age_at_appointment, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs Age at Appointment")

print(amyloid_v_age)



# Scatter plots for each predictor grouped by age
for (predictor in c(
  "Ab40",
  "Ab42",
  "GFAP",
  "NFL",
  "ptau181",
  "ptau231",
  "pTau217",
  "time_diff",
  "race_primary",
  "hispanic",
  "sex"
)) {
  p <-
    ggplot(
      data_time_merge,
      aes_string(x = predictor, y = "Amyloid_PET", color = "age_at_appointment")
    ) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste(
        "Scatter plot of",
        predictor,
        "vs Amyloid PET by Age At Appointment"
      ),
      x = predictor,
      y = "Amyloid PET"
    )
  
  plot_list[[predictor]] <- p
}




# Scatter plots for each predictor grouped by Sex
for (predictor in c(
  "Ab40",
  "Ab42",
  "GFAP",
  "NFL",
  "ptau181",
  "ptau231",
  "pTau217",
  "time_diff",
  "race_primary",
  "hispanic"
)) {
  x <-
    ggplot(data_time_merge,
           aes_string(x = predictor, y = "Amyloid_PET", color = "sex")) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(
      title = paste("Scatter plot of", predictor, "vs Amyloid PET by Sex"),
      x = predictor,
      y = "Amyloid PET"
    )
  
  plot_list[[predictor]] <- x
}


###### BRIAN  #####

# Exploring relationships between PET scores and plasma biomarkers
amyloid_v_ptau231 <- ggplot(data_time_merge, aes(x = ptau231, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs ptau231 Levels")
amyloid_v_ptau181 <- ggplot(data_time_merge, aes(x = ptau181, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs ptau181 Levels")
amyloid_v_ptau217 <- ggplot(data_time_merge, aes(x = pTau217, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs ptau217 Levels")
amyloid_v_ab40 <- ggplot(data_time_merge, aes(x = Ab40, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs Ab40 Levels")
amyloid_v_ab42 <- ggplot(data_time_merge, aes(x = Ab42, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs Ab42 Levels")
amyloid_v_gfap <- ggplot(data_time_merge, aes(x = GFAP, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs GFAP Levels")
amyloid_v_nfl <- ggplot(data_time_merge, aes(x = NFL, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs NFL Levels")

# Checking for any trends with age
amyloid_v_age <- ggplot(data_time_merge, aes(x = age_at_appointment, y = Amyloid_PET)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Amyloid PET vs Age at Appointment")


amyloid_v_age + amyloid_v_ptau181 + amyloid_v_ptau217 + amyloid_v_ptau231 + plot_annotation(title = 'Comparison between Amyloid Measurement and Plasma Biomarkers')
plot_layout(ncol = 2)

amyloid_v_nfl + amyloid_v_gfap + amyloid_v_ab40 + amyloid_v_ab42 + plot_annotation(title = 'Comparison between Amyloid Measurement and Plasma Biomarkers') +
  plot_layout(ncol = 2)



####################


##### Distribution Analysis #####
# Histograms for Amyloid PET Scores
ggplot(data_time_merge, aes(x = Amyloid_PET)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Amyloid PET Scores")

# Histograms for Plasma Biomarkers Variables
biomarkers <- c("Ab40", "Ab42", "GFAP", "NFL", 
                "ptau231", "ptau181", "pTau217")

# List to hold plots
plots_list <- list()

#  Loop through the variables and create histogram plots
for (var in biomarkers) {
  # Create the histogram plot
  p <- ggplot(data_time_merge, aes_string(x = var)) +
    geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
    ggtitle(paste("Distribution of", var)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Add the plot to the list
  plots_list[[var]] <- p
}

# Check if the number of plots is not a perfect square and adjust by adding blank plots
num_plots <- length(plots_list)
num_cols <- ceiling(sqrt(num_plots))  # Number of columns for a square layout

# Create a blank plot to balance visualization
blank_plot <- ggplot() +
  theme_void() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

# Position to insert the blank plot
pos_insert_blank <- length(plots_list) - 1

# Insert the blank plot into the list at the specific position
plots_list <- append(plots_list, list(blank_plot), pos_insert_blank)

# Arrange the plots on a grid
grid.arrange(grobs = plots_list, ncol = num_cols)



# Histogram for Age at Amyloid PET
amyloid_age_hist <- ggplot(data_time_merge, aes(x = age_at_appointment)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Age at Amyloid PET")
# Histogram for Age at Amyloid PET
plasma_age_hist <- ggplot(data_time_merge, aes(x = age_at_acquisition)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Age at Plasma Test")
# Histogram for Age at Amyloid PET
time_diff_hist <- ggplot(data_time_merge, aes(x = time_diff)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Age Differences")

(amyloid_age_hist | plasma_age_hist | time_diff_hist)

# Histogram for sex
ggplot(data_time_merge, aes(x = sex)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Sex")
# Histogram for apoe_e1
apoe1_hist <- ggplot(data_time_merge, aes(x = apoe_e1)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Apoe_e1")
# Histogram for apoe_e2
apoe2_hist <- ggplot(data_time_merge, aes(x = apoe_e2)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Apoe_e2")

(apoe1_hist | apoe2_hist)

# Histogram for Age at Amyloid PET
race_hist <- ggplot(data_time_merge, aes(x = race_primary)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Race")
# Histogram for Age at Amyloid PET
hisp_hist <- ggplot(data_time_merge, aes(x = hispanic)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Hispanic / Non-Hispanic")

(race_hist | hisp_hist)

##### Distribution Analysis #####
# Histograms for Amyloid PET Scores
ggplot(data_time_merge, aes(x = Amyloid_PET)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Amyloid PET Scores")

# Histograms for Plasma Biomarkers Variables
biomarkers <- c("Ab40", "Ab42", "GFAP", "NFL", 
               "ptau231", "ptau181", "pTau217")

# List to hold plots
plots_list <- list()

#  Loop through the variables and create histogram plots
for (var in biomarkers) {
  # Create the histogram plot
  p <- ggplot(data_time_merge, aes_string(x = var)) +
    geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
    ggtitle(paste("Distribution of", var)) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Add the plot to the list
  plots_list[[var]] <- p
}

# Check if the number of plots is not a perfect square and adjust by adding blank plots
num_plots <- length(plots_list)
num_cols <- ceiling(sqrt(num_plots))  # Number of columns for a square layout

# Create a blank plot to balance visualization
blank_plot <- ggplot() +
  theme_void() +
  theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank()
  )

# Position to insert the blank plot
pos_insert_blank <- length(plots_list) - 1

# Insert the blank plot into the list at the specific position
plots_list <- append(plots_list, list(blank_plot), pos_insert_blank)

# Arrange the plots on a grid
grid.arrange(grobs = plots_list, ncol = num_cols)



# Histogram for Age at Amyloid PET
amyloid_age_hist <- ggplot(data_time_merge, aes(x = age_at_appointment)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Age at Amyloid PET")
# Histogram for Age at Amyloid PET
plasma_age_hist <- ggplot(data_time_merge, aes(x = age_at_acquisition)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Age at Plasma Test")
# Histogram for Age at Amyloid PET
time_diff_hist <- ggplot(data_time_merge, aes(x = time_diff)) +
  geom_histogram(bins = 30, fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Age Differences")

(amyloid_age_hist | plasma_age_hist | time_diff_hist)

# Histogram for sex
ggplot(data_time_merge, aes(x = sex)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Sex")
# Histogram for apoe_e1
apoe1_hist <- ggplot(data_time_merge, aes(x = apoe_e1)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Apoe_e1")
# Histogram for apoe_e2
apoe2_hist <- ggplot(data_time_merge, aes(x = apoe_e2)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Apoe_e2")

(apoe1_hist | apoe2_hist)

# Histogram for Age at Amyloid PET
race_hist <- ggplot(data_time_merge, aes(x = race_primary)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Race")
# Histogram for Age at Amyloid PET
hisp_hist <- ggplot(data_time_merge, aes(x = hispanic)) +
  geom_histogram(stat = "count", fill = "cornflowerblue", color = "black") +
  ggtitle("Distribution of Hispanic / Non-Hispanic")

(race_hist | hisp_hist)

##### Relationship Analysis #####
numeric_var <- c("age_at_appointment", "time_diff", "age_at_acquisition", 
                 "ptau231", "ptau181", "pTau217", "Ab40", "Ab42", "GFAP", "NFL", 
                 "apoe_e1", "apoe_e2")



