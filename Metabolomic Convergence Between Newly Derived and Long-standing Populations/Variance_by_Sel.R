# set working directory
setwd("C:/Users/huber/Box/Phillips Lab/Projects/A and C Metabolomic Trajectory Project/Linear Models/")

# load required packages
library(car)
library(dplyr)

# Clear the global environment, but not the plots
rm(list = ls(all.names = TRUE), envir = .GlobalEnv)

#read in data table
data <- read.csv("A-C Normalized Metabolomics Time Series Data.csv", header =T)

#subset data by selection regime
data_A <- subset(data, SimpleSel =="A-Type")
data_C <- subset(data, SimpleSel =="C-Type")


# Levene's test to look at differences in variance

#### A-type

#identify days
days_A <- unique(data_A$Day)

# Initialize a data frame to store results
results_A <- data.frame(Day = character(),
                      Metabolite = character(),
                      Sel_with_larger_variance = character(),
                      P_Value = numeric(),
                      stringsAsFactors = FALSE)

# Loop over each unique day
for (day in days_A) {
  # Subset the main data frame for the current day
  subset_data_A <- data_A[data_A$Day == day, ]
  
  # Loop over each column from 8 to the last column
  for (col_index in 8:ncol(subset_data_A)) {
    metabolite_name_A <- colnames(subset_data_A)[col_index]  # Get metabolite name
    
    # Perform Levene's Test on the current column
    test_result_A <- leveneTest(subset_data_A[[col_index]] ~ subset_data_A$Sel, data = subset_data_A)
    
    # Extract the p-value and store it
    p_value_A <- test_result_A$`Pr(>F)`[1]  # Extracting the p-value from the test result
    
      # Calculate variances for each Sel level
      variances_A <- tapply(subset_data_A[[col_index]], subset_data_A$Sel, var)
      
      # Identify which Sel value has larger variance
      larger_variance_sel_A <- names(variances_A)[which.max(variances_A)]
      
      # Store results for this instance
      results_A <- rbind(results_A, data.frame(Day = day,
                                           Metabolite = metabolite_name_A,
                                           Sel_with_larger_variance_A = larger_variance_sel_A,
                                           P_Value = p_value_A))
    }
  }

# Print the results
print(results_A)



#### C-type

#identify days
days_C <- unique(data_C$Day)

# Initialize a data frame to store results
results_C <- data.frame(Day = character(),
                        Metabolite = character(),
                        Sel_with_larger_variance = character(),
                        P_Value = numeric(),
                        stringsAsFactors = FALSE)

# Loop over each unique day
for (day in days_C) {
  # Subset the main data frame for the current day
  subset_data_C <- data_C[data_C$Day == day, ]
  
  # Loop over each column from 8 to the last column
  for (col_index in 8:ncol(subset_data_C)) {
    metabolite_name_C <- colnames(subset_data_C)[col_index]  # Get metabolite name
    
    # Perform Levene's Test on the current column
    test_result_C <- leveneTest(subset_data_C[[col_index]] ~ subset_data_C$Sel, data = subset_data_C)
    
    # Extract the p-value and store it
    p_value_C <- test_result_C$`Pr(>F)`[1]  # Extracting the p-value from the test result
    
      # Calculate variances for each Sel level
      variances_C <- tapply(subset_data_C[[col_index]], subset_data_C$Sel, var)
      
      # Identify which Sel value has larger variance
      larger_variance_sel_C <- names(variances_C)[which.max(variances_C)]
      
      # Store results for this instance
      results_C <- rbind(results_C, data.frame(Day = day,
                                               Metabolite = metabolite_name_C,
                                               Sel_with_larger_variance_C = larger_variance_sel_C,
                                               P_Value = p_value_C))
    }
  }

# Print the results
print(results_C)


#### combine the results

#rename column header so they can be combined
results_C <- results_C %>% 
  rename(Sel_with_larger_variance = Sel_with_larger_variance_C)

results_A <- results_A %>% 
  rename(Sel_with_larger_variance = Sel_with_larger_variance_A)


# Combine results_A and results_C into single DF
combined_results <- rbind(results_A, results_C)


#### Calculate 1 tailed p-value, testing the hypothesis that AO/nCO have higher variance than ACO/CO
combined_results <- combined_results %>%
  mutate(`1t_pval` = ifelse(Sel_with_larger_variance %in% c("ACO", "CO"),
                            (1 - P_Value) / 2,
                            ifelse(Sel_with_larger_variance %in% c("AO", "nCO"),
                                   P_Value / 2,
                                   NA)))  # NA if none of the conditions are met

#Check for NAs (should be zero)
sum(is.na(combined_results$`1t_pval`))



#### adjust for multiple tests (FDR)

# Create a new column to group AO with ACO, and CO with nCO (selection regime)
combined_results <- combined_results %>%
  mutate(Selection = case_when(
    Sel_with_larger_variance %in% c("AO", "ACO") ~ "A-type",
    Sel_with_larger_variance %in% c("CO", "nCO") ~ "C-type",
    TRUE ~ NA_character_
  ))


# Calculate FDR within each day and selection regime
combined_results <- combined_results %>%
  group_by(Day, Selection) %>%
  mutate(FDR = p.adjust(`1t_pval`, method = "fdr")) %>%
  ungroup()

#Check for sig FDR
sum(combined_results$FDR < 0.01, na.rm = TRUE)

# If there are >0, create a subset of combined_results with FDR < 0.01
combined_results_sigonly <- combined_results %>%
  filter(FDR < 0.01)


# Write results to file
#write.csv(x = combined_results, file = "Sig_Var_results.csv")
