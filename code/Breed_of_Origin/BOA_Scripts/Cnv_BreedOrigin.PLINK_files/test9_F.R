# Load necessary library
library(data.table)

# Input ID file (assuming only one .ID file in the directory)
ID.file <- list.files(pattern = "\\.ID$")
test.name <- sub('\\.ID$', '', basename(ID.file))

# Output .ped file name
out.ped <- paste(test.name, ".ped", sep = "")

# Get the working directory
path <- getwd()

# List of all genotype files (.geno.txt)
file.names <- list.files(path, pattern = "\\.geno\\.txt$", full.names = TRUE)

# Fast reading of genotype files using lapply and fread
geno_list <- lapply(file.names, fread, header = FALSE, sep = " ", data.table = FALSE)

# Combine all genotype data using do.call and cbind (efficient binding)
geno_combined <- do.call(cbind, geno_list)

# Fast reading of the ID file using fread
ID_data <- fread(ID.file, header = FALSE, sep = " ", data.table = FALSE)

# Combine the ID data with the genotype data
X2 <- cbind(ID_data, geno_combined)

# Write the combined data to the output .ped file
write.table(X2, out.ped, quote = FALSE, row.names = FALSE, col.names = FALSE)

# Output the completion message with the output file name
cat("Output saved to:", out.ped, "\n")