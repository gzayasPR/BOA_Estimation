library(dplyr)
library(tidyr)
library(readr)

# Set the base directory
base_dir <- getwd()
chromosomes <- sprintf("%02d", 1:29)   # Adjust as per your data

# Read parameters from the CSV file
parameters_df <- read_csv('parameters.csv')
folds <- 0:4  # Assuming 5 folds


# Initialize empty lists to store the results
Angus_results <- list()
Brahman_results <- list()
combined_results <- list()
for (chr_num in chromosomes) {
  Angus.df <- data.frame()
  Brahman.df <- data.frame()
    for (fold in folds) {
      Angus_file <- sprintf("%s/Fold_%02d/Angus_Chr%s.txt", base_dir, fold, chr_num)
      Brahman_file <- sprintf("%s/Fold_%02d/Brahman_Chr%s.txt", base_dir, fold, chr_num)
      Angus.df <- rbind(Angus.df,read.table(Angus_file,header=T))
      Brahman.df <- rbind(Brahman.df,read.table(Brahman_file,header=T))}
  Angus_results[[chr_num]] <- Angus.df 
  Brahman.df[,-c(1:2)] <- 1 - Brahman.df[,-c(1:2)]
  Brahman_results[[chr_num]]  <- Brahman.df 
  combined_results[[chr_num]]  <- rbind(Brahman.df,Angus.df )
      }


# Function to calculate mean and sd, and create the formatted string
calculate_mean_sd <- function(df, colname) {
  mean_val <- mean(df[[colname]], na.rm = TRUE)
  sd_val <- sd(df[[colname]], na.rm = TRUE)
  return(sprintf("%.3f ± %.2f", mean_val, sd_val))
}

# Function to apply the mean and sd calculation across all columns for each dataframe in the list
combine_results <- function(results_list) {
  output_df <- data.frame(CHR = 1:29)
  for (chr_num in 1:29) {
    chr_df <- results_list[[chr_num]]
    for (colname in colnames(chr_df)[-c(1:2)]) { # Skip the first two columns if they are not hyperparameter data
      output_df[chr_num, colname] <- calculate_mean_sd(chr_df, colname)
    }
  }
  return(output_df)
}

# Combine the results for Brahman
Brahman_combined <- combine_results(Brahman_results)

# Combine the results for Angus
Angus_combined <- combine_results(Angus_results)

# Print out the combined results
print(Brahman_combined)
print(Angus_combined)

# After you've created your final combined DataFrame
readr::write_csv(Brahman_combined, "Brahman.CV.csv")
readr::write_csv(Angus_combined, "Angus.CV.csv")


library(dplyr)
library(tidyr)
library(purrr)

# Adjusted function to rank and return the top n parameter combinations for a given chromosome
get_best_parameters <- function(results_list, chr, top_n = 1) {
  # Extract the dataframe for the specified chromosome
  chr_df <- results_list[[as.character(chr)]]
  
  # Calculate the mean for each parameter across individuals and gather into long format
  parameters_summary <- chr_df %>%
    select(-PID, -IID) %>%  # Exclude non-parameter columns
    summarise_all(mean, na.rm = TRUE) %>%
    gather(key = "Parameter", value = "Mean") %>%
    mutate(CHR = chr) %>%
    # Convert Parameter back into separate Window and Hidden columns
    separate(Parameter, into = c("Window", "Hidden"), sep = "_") %>%
    # Arrange by descending mean values
    arrange(desc(Mean)) %>%
    # Ensure only the top n results are returned
    slice(1:top_n) %>%
    # Create new columns names based on rank and pivot wider
    mutate(Rank = row_number()) %>%
    pivot_wider(names_from = Rank, 
                names_glue = "{.value}{Rank}", 
                values_from = c(Window, Hidden, Mean))
  
  return(parameters_summary)
}
chromosomes <- sprintf("%02d", 1:29) 
# Function to apply get_best_parameters across all chromosomes
get_all_best_parameters <- function(results_list, top_n = 1) {
  # Apply get_best_parameters to each chromosome and store the results in a list
  results <- map_dfr(chromosomes,  ~get_best_parameters(results_list, .x, top_n))
  return(results)
}

# Example usage:
# This will get the best parameters for all chromosomes
best_params_Angus <- get_all_best_parameters(Angus_results, top_n = 3)
best_params_Brahman <- get_all_best_parameters(Brahman_results, top_n = 3)

readr::write_csv(best_params_Angus , "Angus.best_parmaeters.csv")
readr::write_csv(best_params_Brahman ,"Brahman.best_parmaeters.csv")


best_params_combined <- get_all_best_parameters(combined_results, top_n = 3)
readr::write_csv(best_params_combined ,"Combined.best_parameters.csv")



library(dplyr)
library(tidyr)
library(readr)

# Set the base directory
base_dir <- getwd()
chromosomes <- sprintf("%02d", 1:29)   # Adjust as per your data

# Read parameters from the CSV file
parameters_df <- read_csv('parameters.csv')
folds <- 0:4  # Assuming 5 folds


# Initialize empty lists to store the results
Angus_results <- list()
Brahman_results <- list()
combined_results <- list()
for (chr_num in chromosomes) {
  Angus.df <- data.frame()
  Brahman.df <- data.frame()
    for (fold in folds) {
      Angus_file <- sprintf("%s/Fold_%02d/Angus_Chr%s.txt", base_dir, fold, chr_num)
      Brahman_file <- sprintf("%s/Fold_%02d/Brahman_Chr%s.txt", base_dir, fold, chr_num)
      Angus.df <- rbind(Angus.df,read.table(Angus_file,header=T))
      Brahman.df <- rbind(Brahman.df,read.table(Brahman_file,header=T))}
  Angus_results[[chr_num]] <- Angus.df 
  Brahman.df[,-c(1:2)] <- 1 - Brahman.df[,-c(1:2)]
  Brahman_results[[chr_num]]  <- Brahman.df 
  combined_results[[chr_num]]  <- rbind(Brahman.df,Angus.df )
      }

all_combined_results <- bind_rows(combined_results)

library(dplyr)
library(tidyr)
library(stringr)

# Function to calculate mean, sd, and select the top n parameters
get_top_parameters <- function(data, top_n = 1) {
  # First, gather the data into a long format
  long_data <- data %>%
    gather(key = "Parameter", value = "Value", -PID, -IID)
  
  # Now, calculate the mean and sd for each parameter
  data_summary <- long_data %>%
    group_by(Parameter) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD = sd(Value, na.rm = TRUE)
    ) %>%
    ungroup()

  # Extract Window and Hidden from Parameter
  data_summary <- data_summary %>%
    mutate(
      Window = as.numeric(gsub("W(\\d+)_H\\d+", "\\1", Parameter)),
      Hidden = as.numeric(gsub("W\\d+_H(\\d+)", "\\1", Parameter))
    )

  # Arrange by descending Mean
  top_parameters <- data_summary %>%
    arrange(desc(Mean)) %>%
    slice_head(n = top_n)
  
  return(top_parameters)
}

# Example usage
# Get the best parameter across all chromosomes
best_parameter <- get_top_parameters(all_combined_results, top_n = 1)

# Get the top 3 parameters across all chromosomes
top_3_parameters <- get_top_parameters(all_combined_results, top_n = 3)

# Get the top 3 parameters across all chromosomes
all_parameters <- get_top_parameters(all_combined_results, top_n = nrow(parameters_df))

# Print results
print(best_parameter)
print(top_3_parameters)

readr::write_csv(all_parameters ,"Combined.ALL_parameters.csv")