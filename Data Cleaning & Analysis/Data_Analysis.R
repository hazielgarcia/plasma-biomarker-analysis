##################################
##########   STAT 409   ##########
####  Data Science Practicum  ####
##################################
########## DATA CLEANING #########
##################################

# Clear environment
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(corrplot)
library(patchwork)

# Load the datasets
participants <- read.csv("participants_anon.csv")
pib_dvr_aal <- read.csv("pib_dvr_aal_anon.csv") #Not for primary use
amyloid_pet <- read.csv("Amyloid_PET_anon.csv")
plasma <- read.csv("plasma_simoa_gothenburg_anon.csv")
cog_status <- read.csv("consensus_diagnosis_anon.csv")


# Display the structure of the datasets to understand their variables
# Generate summary statistics for each dataset

## PIB_DVR_AAL Data
str(amyloid_pet)
str(pib_dvr_aal)  # For comparison to amyloid_pet
#View(pib_dvr_aal)

## Amyloid PET Data
summary(amyloid_pet)
#View(amyloid_pet)
count(distinct(amyloid_pet["subject_id_number"])) # Number of dist. subjects
count(distinct(pib_dvr_aal["subject_id_number"])) # Check against pib data


## Participants (Demographic Data)
str(participants)
summary(participants)
#View(participants)
count(distinct(participants["subject_id_number"])) # Number of Distinct Subjects
# Analyze categorical columns with data
table(participants["apoe_e1"])
table(participants["apoe_e2"])
table(participants["sex"])
table(participants["race_primary"])
table(participants["hispanic"])
table(participants["ed_years"])


## Plasma Biomarker Data
str(plasma)
summary(plasma)
#View(plasma)
count(distinct(plasma["subject_id_number"])) # Number of dist. subjects

## Cognitive Status at time of visit data
str(cog_status)
# View(cog_status)
summary(cog_status)
table(cog_status["Calculated_Consensus_dx"])

# Calculate the proportion of missing data in each column of each dataset
sapply(participants, function(x)
  sum(is.na(x)) / length(x))
sapply(plasma, function(x)
  sum(is.na(x)) / length(x))
sapply(amyloid_pet, function(x)
  sum(is.na(x)) / length(x))
sapply(cog_status, function(x)
  sum(is.na(x)) / length(x))


### Creating plots of biomarkers over time ###
# Re-load plasma (if necessary)
plasma <- read.csv("plasma_simoa_gothenburg_anon.csv")

# Ensure all biomarker columns are numeric and sort data by subject and time
plasma_clean <- plasma %>%
  mutate(
    ptau181 = as.numeric(as.character(ptau181)),
    Ab42_40_ratio = Ab42 / Ab40  # Creating the Ab42/Ab40 ratio column
  ) %>%
  select(
    subject_id_number,
    age_at_acquisition,
    ptau231,
    ptau181,
    pTau217,
    Ab42_40_ratio,  # Include the new ratio column
    GFAP,
    NFL
  ) %>%
  arrange(subject_id_number, age_at_acquisition) %>%  # Sorting by subject & age
  na.omit()

# Filtering subjects with more than one test and
# Calculating the age difference from the first test
multi_visit_data <- plasma_clean %>%
  group_by(subject_id_number) %>%
  filter(n() > 1) %>%
  mutate(age_diff = age_at_acquisition - first(age_at_acquisition)) %>%
  mutate(across(c(
    ptau231, ptau181, pTau217, Ab42_40_ratio, GFAP, NFL
  ), ~ . - first(.))) %>%
  ungroup()

# Generate individual plots for each biomarker
biomarkers <-
  c("ptau231", "ptau181", "pTau217", "Ab42_40_ratio", "GFAP", "NFL")

plot_list <- lapply(biomarkers, function(biomarker) {
  ggplot(
    multi_visit_data,
    aes(
      x = age_diff,
      y = !!sym(biomarker),
      group = subject_id_number,
      color = as.factor(subject_id_number)
    )
  ) +
    geom_line(show.legend = FALSE) +
    geom_point(show.legend = FALSE) +
    labs(
      title = paste("Changes in", biomarker, "Over Time Relative to First Test"),
      x = "Years Since First Test",
      y = paste("Change in", biomarker)
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
})

# Print the plots
plot_list
plasma <- read.csv("plasma_simoa_gothenburg_anon.csv")


# Make Histograms of variables to visually analyze pre-cleaned distribution
# Amyloid PET Data
init_amyloid_plot <-
  ggplot(amyloid_pet %>% filter(is.finite(Amyloid_PET)), aes(x = Amyloid_PET)) +
  geom_histogram(bins = 30,
                 fill = "blue",
                 color = "black") +
  ggtitle("Initial Distribution of Amyloid values")

# Participants Data
init_ed_plot <-
  ggplot(participants %>% filter(is.finite(ed_years)), aes(x = ed_years)) +
  geom_histogram(bins = 30,
                 fill = "red",
                 color = "black") +
  ggtitle("Initial Distribution of ed_years")

init_apoe1_plot <-
  ggplot(participants %>% filter(is.finite(apoe_e1)), aes(x = apoe_e1)) +
  geom_histogram(bins = 30,
                 fill = "green",
                 color = "black") +
  ggtitle("Initial Distribution of apoe_e1 values")

init_apoe2_plot <-
  ggplot(participants %>% filter(is.finite(apoe_e2)), aes(x = apoe_e2)) +
  geom_histogram(bins = 30,
                 fill = "orange",
                 color = "black") +
  ggtitle("Initial Distribution of apoe_e2 values")

# Plasma Data
init_gfap_plot <-
  ggplot(plasma %>% filter(is.finite(GFAP)), aes(x = GFAP)) +
  geom_histogram(bins = 30,
                 fill = "pink",
                 color = "black") +
  ggtitle("Initial Distribution of GFAP values")

init_nfl_plot <-
  ggplot(plasma %>% filter(is.finite(NFL)), aes(x = NFL)) +
  geom_histogram(bins = 30,
                 fill = "grey",
                 color = "black") +
  ggtitle("Initial Distribution of NFL values")

init_ab42_plot <-
  ggplot(plasma %>% filter(is.finite(Ab42)), aes(x = Ab42)) +
  geom_histogram(bins = 30,
                 fill = "lightblue",
                 color = "black") +
  ggtitle("Initial Distribution of AB42 values")

init_ab40_plot <-
  ggplot(plasma %>% filter(is.finite(Ab40)), aes(x = Ab40)) +
  geom_histogram(bins = 30,
                 fill = "orangered",
                 color = "black") +
  ggtitle("Initial Distribution of AB40 values")

init_ptau231_plot <-
  ggplot(plasma %>% filter(is.finite(ptau231)), aes(x = ptau231)) +
  geom_histogram(bins = 30,
                 fill = "cyan",
                 color = "black") +
  ggtitle("Initial Distribution of ptau231 values")

plasma$ptau181 <- as.numeric(as.character(plasma$ptau181))
init_ptau181_plot <-
  ggplot(plasma %>% filter(is.finite(ptau181)), aes(x = ptau181)) +
  geom_histogram(bins = 30,
                 fill = "lightgreen",
                 color = "black") +
  ggtitle("Initial Distribution of ptau181 values")

init_ptau217_plot <-
  ggplot(plasma %>% filter(is.finite(pTau217)), aes(x = pTau217)) +
  geom_histogram(bins = 30,
                 fill = "grey",
                 color = "black") +
  ggtitle("Initial Distribution of ptau217 values")

# Combine and display plots
(init_amyloid_plot + init_ed_plot + init_apoe1_plot + init_apoe2_plot) + plot_layout(ncol = 2)
(
  init_gfap_plot + init_nfl_plot + init_ptau181_plot + init_ptau231_plot + init_ptau217_plot
) + plot_layout(ncol = 2)



### Visualizations for cognitive status:
# Create a bar plot for cognitive status distribution
ggplot(cog_status, aes(x = Calculated_Consensus_dx)) +
  geom_bar(fill = "skyblue") +
  theme_minimal() +
  labs(title = "Distribution of Cognitive Status Categories",
       x = "Cognitive Status",
       y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


