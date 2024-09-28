# List of required packages
required_packages <- c("dplyr", "readr", "tidyr", "tidyverse", "vcfR", "data.table")

# Function to check and install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install missing packages
lapply(required_packages, install_if_missing)

# Load the libraries
lapply(required_packages, library, character.only = TRUE)
