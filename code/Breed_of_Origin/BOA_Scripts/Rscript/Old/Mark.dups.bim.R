# Load necessary library
if (!require(data.table)) install.packages("data.table")
library(data.table)

# Read the BIM file
bim_file <- "UFID_UF_ILUMINA.ID_OCT23_Illumina.ID_ARS.bim"  # Replace with your BIM file path
bim_data <- fread(bim_file, header = FALSE, col.names = c('chr', 'snp', 'gd', 'pos', 'a1', 'a2'))

# Identifying and adjusting duplicate positions
bim_data[, pos := pos + seq_len(.N) - 1, by = .(chr, pos)]

# Write the modified BIM file
fwrite(bim_data, file = "UFID_UF_ILUMINA.ID_OCT23_Illumina.ID_ARS.bim", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)