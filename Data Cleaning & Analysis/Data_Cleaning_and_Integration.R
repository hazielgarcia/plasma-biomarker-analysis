##################################
##########   STAT 409   ##########
####  Data Science Practicum  ####
##################################
## DATA CLEANING & INTEGRATION  ##
##################################

### FINAL CLEANED DATA after this has run ###
# 'amyloid_pet'
# 'participants'
# 'plasma'

### FINAL MERGED DATASET ##
# 'data_time_merge'


# <<<  YOUR WORKING DIRECTORY HERE  >>>

# Clear environment
rm(list = ls())

# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(corrplot)

# Load the datasets
participants <- read.csv("participants_anon.csv")
pib_dvr_aal <- read.csv("pib_dvr_aal_anon.csv") #Not for primary use
amyloid_pet <- read.csv("Amyloid_PET_anon.csv")
plasma <- read.csv("plasma_simoa_gothenburg_anon.csv")
cog_status <- read.csv("consensus_diagnosis_anon.csv")


# Display the structure of the datasets &Generate summary statistics
## Participants (Demographic Data)
#View(participants)
str(participants)
summary(participants)

count(distinct(participants["subject_id_number"])) # Number of Distinct Subjects
table(participants["race_primary"])
table(participants["hispanic"])


## Plasma Biomarker Data
#View(plasma)
str(plasma)
summary(plasma)

count(distinct(plasma["subject_id_number"])) # Number of dist. subjects


## Amyloid PET Data
# str(pib_dvr_aal)  # For comparison to amyloid_pet
#View(amyloid_pet)
str(amyloid_pet)
summary(amyloid_pet)

count(distinct(amyloid_pet["subject_id_number"])) # Number of dist. subjects
count(distinct(pib_dvr_aal["subject_id_number"])) # Check against pib data


## Cognitive Status
#View(cog_status)
str(cog_status)
summary(cog_status)

count(distinct(cog_status["subject_id_number"])) # Check against pib data
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



#####################################################
#################  DATA CLEANING  ###################
#####################################################
##### Dropping NA Data as Discussed with Client #####


#####  PARTICIPANTS Data Cleaning  #####
participants <- read.csv("participants_anon.csv")
str(participants)
summary(participants)

# Seems as if duplicate row in 'participants' verify & remove if true
# Checking for duplicates based on 'subject_id_number', which should be unique
id_duplicates <-
  participants[duplicated(participants$subject_id_number) |
                 duplicated(participants$subject_id_number, fromLast = TRUE),]
print(id_duplicates)

# Removing exact duplicated line
participants <- participants[!duplicated(participants[, -1]), ]

# Verify by checking the number of rows and looking for duplicates again
print(nrow(participants))
print(any(duplicated(participants)))
print(any(duplicated(participants$subject_id_number)))

# Drop 'X' and NA columns
participants <-
  participants %>% select(-c(X, race_secondary, race_tertiary, hisp_orig))

# Drop 'ed_years' column per client advice
participants <-
  participants %>% select(-c(ed_years))

# Verify Results
str(participants)
summary(participants)

# Remaining NA's in 'apoe_e1', & 'apoe_e2'
# Select columns that are relevant or that you know contain NA values
# Apply a function across rows to count NAs
na_count_per_row <-
  apply(participants[, c("apoe_e1", "apoe_e2")], 1, function(x)
    sum(is.na(x)))

# View rows where the count of NAs is greater than 0
participants[na_count_per_row > 0,]

# Table to show how many rows have exactly 1 NA, 2 NAs, etc.
table(na_count_per_row)

# Remove rows with 'apoe_e1' & 'apoe_e2' with NA
participants <-
  participants[!is.na(participants$apoe_e1) &
                 !is.na(participants$apoe_e2), ]
# Adding 'apoe_ternary' --> int: 0, 1, or 2
##  Column for number of '4' values in apoe_e1' & 'apoe_e2'
participants <- participants %>%
  mutate(apoe_ternary = (apoe_e1 == 4) + (apoe_e2 == 4))
# Adding 'apoe_binary' --> int: 0 or 1
##  Column for if any '4' values in apoe_e1' & 'apoe_e2'
participants <- participants %>%
  mutate(apoe_binary = ifelse(apoe_e1 == 4 | apoe_e2 == 4, 1, 0))


# Check if issues with blanks & NA's in 'race_primary' & 'hispanic'
# have been resolved with other row removals
table(participants["race_primary"])
table(participants["hispanic"])
# There are no longer any blank or "Unknown" values in these columns

# Check if all instances where race_primary is "Other" have hispanic as "yes"
check <- participants %>%
  filter(race_primary == "Other") %>%
  summarize(all_hispanic_yes = all(hispanic == "yes"))

# Print the result
print(check)
# ALL INSTANCES MATCH. WE WILL NOT REMOVE 'hispanic' FROM THE DATASET,
# BUT IT WILL NOT BE INCLUDED IN OUR MODELS TO AVOID SINGULARITY ISSUES.

# Convert categorical columns to factors
participants$sex <- as.factor(participants$sex)
participants$race_primary <- as.factor(participants$race_primary)
participants$hispanic <- as.factor(participants$hispanic)
participants$apoe_e1 <- as.factor(participants$apoe_e1)
participants$apoe_e2 <- as.factor(participants$apoe_e2)
participants$apoe_ternary <- as.factor(participants$apoe_ternary)
participants$apoe_binary <- as.factor(participants$apoe_binary)


# Verify Results
str(participants)
summary(participants)
count(distinct(participants["subject_id_number"])) # Number of Distinct Subjects
#View(participants)



##### PLASMA Data Cleaning #####
plasma <- read.csv("plasma_simoa_gothenburg_anon.csv")
str(plasma)
summary(plasma)
count(distinct(plasma["subject_id_number"]))

## Coerce 'ptau181' from character to numeric ##
plasma$ptau181 <- suppressWarnings(as.numeric(as.character(plasma$ptau181)))
##  Remove column 'X'  ##
plasma <- plasma %>% select(-c(X))
# Verify Results
summary(plasma)

# Checking relationship between 'ptau231' and 'ptau231_text'
# Check if all non-NA values are the same in 'ptau231' and 'ptau231_text'
all_equal <-
  all(plasma$ptau231[!is.na(plasma$ptau231) &
                       !is.na(plasma$ptau231_text)] ==
        plasma$ptau231_text[!is.na(plasma$ptau231) &
                              !is.na(plasma$ptau231_text)])
# Print the result
print(all_equal)

# The Non-NA values all match
# 'ptua_231' has 1 NA, while 'ptau231_text' has 2
# Check the NA in 'ptau_231' is also in 'ptau_231_text'
# Identify rows where 'ptau231' is NA
na_in_ptau231 <- is.na(plasma$ptau231)
# Check if the same rows are NA in 'ptau231_text'
(matching_na <-
    all(na_in_ptau231 == (
      na_in_ptau231 & is.na(plasma$ptau231_text)
    )))
# This NA is in both columns
# It seems dropping 'ptau231_text' will lose no information

## Drop 'ptau231_text'  ##
if (matching_na) {
  plasma$ptau231_text <-
    NULL  # Drop 'ptau231_text' from the dataframe
} else {
  print("NA values do not align perfectly. Reevaluate before dropping.")
}

# Verify Results
str(plasma)
summary(plasma)

##  Remove column 'pTau217_cv'  ##
# Column contains coeff. variation produced by previous experiments
# This information is not germane to our study
plasma <- plasma %>% select(-c(pTau217_cv))
# Verify Results
str(plasma)
summary(plasma)

## Dropping Rows With NA
# Count the number of NA values in each row of the plasma data frame
na_count_per_row <- rowSums(is.na(plasma))
# Display the count of NA values for rows that have at least one NA
na_count_per_row[na_count_per_row > 0]
# Filter to see rows with more than 1 NA
plasma_with_multiple_nas <- plasma[na_count_per_row > 1,]
# View these specific rows
print(plasma_with_multiple_nas)
# Get a summary table of NA counts per row
table(na_count_per_row)

# Drop the row with 6 NA's & Check Data
plasma <- plasma[rowSums(is.na(plasma)) != 6,]
summary(plasma)

# This Leaves NA values in:
#  'ptau181':6  & 'pTau217':6

# Re-check the distribution of NAs
na_count_per_row <- rowSums(is.na(plasma))
table(na_count_per_row)

# Remove rows with NAs in 'ptau181' & 'ptau217'; Check Data
plasma <- plasma[!is.na(plasma$pTau217),]
plasma <- plasma[!is.na(plasma$ptau181),]
summary(plasma)
str(plasma)
#View(plasma)

# Found possible duplicate rows for subjects
# Identify indices of all duplicated rows
dup_indices <-
  which(duplicated(plasma) | duplicated(plasma, fromLast = TRUE))

# Extract the duplicated rows
duplicates <- plasma[dup_indices,]

# Print the duplicated rows along with their duplicates
print(duplicates)

# Duplicates found for subject 1; Removing duplicate rows
plasma <- plasma %>% filter(!duplicated(.))

# Create the 'Ab_ratio' column
plasma <- plasma %>%
  mutate(Ab_ratio = Ab42 / Ab40)

# Verify results
str(plasma)
summary(plasma)
count(distinct(plasma["subject_id_number"])) # Number of Distinct Subjects



#####  Amyloid PET Data Cleaning  #####
amyloid_pet <- read.csv("Amyloid_PET_anon.csv")
str(amyloid_pet)
summary(amyloid_pet)
count(distinct(amyloid_pet["subject_id_number"])) # Number of Distinct Subjects

# Removing rows with NAs
amyloid_pet <- na.omit(amyloid_pet)
# Remove column 'X'
amyloid_pet <- amyloid_pet %>% select(-c(X))

# Adding 'amyloid_binary' column with a threshold of 1.17
amyloid_pet$amyloid_binary <-
  ifelse(amyloid_pet$Amyloid_PET > 1.17, 1, 0)


# Checking the results
table(amyloid_pet$amyloid_binary)

# Verify Results
str(amyloid_pet)
summary(amyloid_pet)
count(distinct(amyloid_pet["subject_id_number"])) # Number of Distinct Subjects



#####  Cognitive Status PET Data Cleaning  #####
cog_status <- read.csv("consensus_diagnosis_anon.csv")
str(cog_status)
summary(cog_status)
count(distinct(cog_status["subject_id_number"])) # Number of Distinct Subjects

# Remove column 'X'
cog_status <- cog_status %>% select(-c(X))

table(cog_status["Calculated_Consensus_dx"])

# Group categories of 'Calculated_Consensus_dx' per client suggestions
cog_status <- cog_status %>%
  mutate(
    Calculated_Consensus_dx = case_when(
      Calculated_Consensus_dx == 
        "Cog_Unimpaired_Stable" ~ "Cog_Unimpaired_Stable",
      Calculated_Consensus_dx %in% 
        c("Cog_Unimpaired_Declining", "Impaired_Not_MCI") 
      ~ "Unimpaired_Declining_or_Impaired",
      Calculated_Consensus_dx %in% 
        c("Clinical_MCI", "Dementia") ~ "MCI_or_Dementia",
      Calculated_Consensus_dx == 
        "No_Diagnosis_Calculated" ~ "No_Diagnosis_Calculated",
      TRUE ~ as.character(Calculated_Consensus_dx)
    )
  )

# Convert 'Calculated_Consensus_dx' to factor
cog_status$Calculated_Consensus_dx <-
  as.factor(cog_status$Calculated_Consensus_dx)

# View the updated table
table(cog_status$Calculated_Consensus_dx)

# Verify Results
str(cog_status)
summary(cog_status)
count(distinct(cog_status["subject_id_number"])) # Number of Distinct Subjects


##################################
#######  DATA INTEGRATION  #######
########### TIME-MERGE ###########
##################################

### TIME-ALIGNED MERGE - Combining by nearest ages in tests  ###
# Ensure only subjects with both tests are considered
common_ids <-
  intersect(
    intersect(amyloid_pet$subject_id_number, plasma$subject_id_number),
    participants$subject_id_number
  )

str(common_ids)

# Filter both datasets to include only these common subjects
filtered_amyloid_pet <-
  amyloid_pet[amyloid_pet$subject_id_number %in% common_ids,]
filtered_plasma <-
  plasma[plasma$subject_id_number %in% common_ids,]

# Define function to find nearest plasma record for each PET record, by subject
nearest_plasma <- function(pet_row, plasma_df) {
  # Subset the plasma dataframe for the current PET row's subject ID
  subject_plasma_df <-
    plasma_df[plasma_df$subject_id_number == pet_row$subject_id_number,]
  # Calculate the signed time differences
  time_diffs_signed <-
    pet_row$age_at_appointment - subject_plasma_df$age_at_acquisition
  # Calculate the absolute time differences
  time_diffs_abs <- abs(time_diffs_signed)
  # Find the index of the minimum time difference
  closest_index <- which.min(time_diffs_abs)
  # Select the closest plasma record
  closest_plasma <- subject_plasma_df[closest_index,]
  # Ensure that closest_plasma is a flat dataframe
  closest_plasma <- data.frame(t(unlist(closest_plasma)))
  # Combine the PET and selected plasma record
  combined_row <-
    c(pet_row, closest_plasma[setdiff(names(closest_plasma), 
                                      "subject_id_number")])
  # Convert combined vector to dataframe
  combined_df <- data.frame(t(unlist(combined_row)))
  # Assign column names appropriately
  colnames(combined_df) <-
    c(names(pet_row), setdiff(names(closest_plasma), "subject_id_number"))
  # Add the time difference as a new column
  combined_df$time_diff <- time_diffs_signed[closest_index]
  return(combined_df)
}

# Apply function to each PET record using lapply
aligned_data <-
  do.call(rbind, lapply(1:nrow(filtered_amyloid_pet), function(i)
    nearest_plasma(filtered_amyloid_pet[i,], filtered_plasma)))

# Convert the result into a proper dataframe and ensure types are correct
aligned_data <- as.data.frame(aligned_data)
aligned_data <- type.convert(aligned_data, as.is = TRUE)

# Incorporating demographic data
data_time_merge <-
  left_join(aligned_data, participants, by = "subject_id_number")

# View the structure and summary of the final integrated data
# View(data_time_merge)
str(data_time_merge)
summary(data_time_merge)

# Number of Distinct Subjects
count(distinct(data_time_merge["subject_id_number"])) 

# Append cognitive status information
# Ensure the column names match for merging
# Rename VisNo to visno in cog_status
cog_status <- cog_status %>%
  rename(visno = VisNo)

# Merge cog_status with the existing integrated dataset 
data_time_merge <-
  left_join(data_time_merge,
            cog_status,
            by = c("subject_id_number", "visno"))

# View the structure and summary of the final integrated data 
str(data_time_merge)
summary(data_time_merge)

# Number of distinct subjects 
count(distinct(data_time_merge["subject_id_number"]))


# Filter the data to eliminate rows where abs(time_diff) > 2
data_time_merge <- data_time_merge[abs(data_time_merge$time_diff) <= 2, ]
table(data_time_merge[,"Calculated_Consensus_dx"])

# Convert response variable, 'amyloid_binary', to factor
#data_time_merge$amyloid_binary <- as.factor(data_time_merge$amyloid_binary)

# Display the structure of the filtered data to confirm the changes
str(data_time_merge)
summary(data_time_merge)
#View(data_time_merge)
# Number of distinct subjects 
count(distinct(data_time_merge["subject_id_number"]))

# List of plasma biomarkers predictor variables
plasma_predictors <- 
  c("ptau231", "ptau181", "pTau217", "GFAP", "NFL", "Ab_ratio")
# Standardize plasma predictors
data_time_merge[plasma_predictors] <- scale(data_time_merge[plasma_predictors])