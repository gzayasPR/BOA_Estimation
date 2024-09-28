 rm(list=ls()) 
library(dplyr)
library(readr)
library(tidyr)
# Set the base directory
base_dir <- getwd()  # Adjust as per your data location
chromosomes <- sprintf("%02d", 1:29)
folds <- 0:4  # Assuming 5 folds

# Initialize a dataframe to store combined results
all_results <- data.frame()
parameters <- read.csv("parameters.csv")

# Process data for each chromosome and fold
for (chr_num in chromosomes) {
  for (fold in folds) {
    file_pattern <- sprintf("*_Chr%s.txt", chr_num)
    fold_dir <- sprintf("%s/Fold_%02d/", base_dir, fold)
    file_paths <- list.files(path = fold_dir, pattern = file_pattern, full.names = TRUE)
    
    if (length(file_paths) == 0) next  # Skip if no files found
    
    F1.IDs <- read_table(sprintf("%s/F1.BO.ID", fold_dir ))
    F1.IDs$Group <- "F1"
    pop1.IDs <- read_table(sprintf("%s/pop1.BO.ID", fold_dir ))
    pop1.IDs$Group <- "pop1"
    pop2.IDs <- read_table(sprintf("%s/pop2.BO.ID", fold_dir ))
    pop2.IDs$Group <- "pop2"
    IDs <- rbind(F1.IDs,pop1.IDs,pop2.IDs)
    IDs$PID <- as.character(IDs$PID )
    IDs$IID <- as.character(IDs$IID )
    # Read and combine breed data from files
    df <- read_table(file_paths)
    df.merged <- merge(df,IDs,by = c("PID","IID"))
    # Calculate statistics for each group and parameter
    stats <- df.merged %>%
      gather(key = "parameter", value = "ancestry", -PID, -IID, -Group) %>%
      separate(parameter, into = c("window", "hidden_state"), sep = "_") %>%
      group_by(Group, window, hidden_state) %>%
      summarize(
        mean_ancestry = mean(ancestry, na.rm = TRUE),
        sd_ancestry = sd(ancestry, na.rm = TRUE),
        .groups = 'drop'
      )

    # Combine results
    all_results <- bind_rows(all_results, stats)
  }
}

all_results <- all_results %>%
  mutate(
    Parameters = paste(window, hidden_state, sep = "_"),
    Window = as.numeric(sub("W", "", window)),  # Extract numeric part from window
    Hidden_States = as.numeric(sub("H", "", hidden_state))  # Extract numeric part from hidden_state
  )


# Now pivot the data to widen it
all_results_wide <- all_results %>%
  pivot_wider(
    id_cols = c("Parameters", "Window", "Hidden_States"),
    names_from = Group,
    names_sep = "_",
    values_from = c(mean_ancestry, sd_ancestry),
    values_fn = list(mean_ancestry = mean, sd_ancestry = mean)  # Ensure to use appropriate aggregation if necessary
  )

# Rename the columns to match your desired output format
colnames(all_results_wide) <- gsub("mean_ancestry", "mean", colnames(all_results_wide))
colnames(all_results_wide) <- gsub("sd_ancestry", "sd", colnames(all_results_wide))


col_order <- c("Parameters", "Window", "Hidden_States","mean_pop1","sd_pop1","mean_pop2","sd_pop2","mean_F1","sd_F1")
all_results_wide <- all_results_wide[col_order]
# View the restructured data
print(all_results_wide)
# Write results to CSV
write_csv(all_results_wide, "final_hyperparameter_results.csv")

