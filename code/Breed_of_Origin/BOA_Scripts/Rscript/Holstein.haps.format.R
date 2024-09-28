library(vcfR)

# Load VCF file
vcf_file <- "Holstein.phased.vcf"  # Replace with your VCF file path
vcf <- read.vcfR(vcf_file)

# Extract genotype data
genotypes <- extract.gt(vcf, element = "GT")
genotypes[genotypes == "NA"] <- "?|?"
# Replace genotype codes (0 and 1) with actual alleles
replace_alleles <- function(genotypes, ref, alt) {
  sapply(seq_len(nrow(genotypes)), function(i) {
    gt <- genotypes[i, ]
    gt <- gsub("0[|/]", paste0(ref[i], "|"), gt)
    gt <- gsub("[|/]0", paste0("|", ref[i]), gt)
    gt <- gsub("1[|/]", paste0(alt[i], "|"), gt)
    gt <- gsub("[|/]1", paste0("|", alt[i]), gt)
    gt
  })
}

map <- as.data.frame(vcf@fix)
# Apply the function to all genotypes
haplotypes <- replace_alleles(genotypes, map$REF, map$ALT)

# Split the haplotypes into two rows per individual
split_haplotypes <- function(haplotypes) {
  t(apply(haplotypes, 2, function(x) {
    unlist(lapply(strsplit(x, "[|/]"), function(y) c(y[1], y[2])))
  }))
}

ID.vcf <- row.names(haplotypes)
# Repeat each element twice
repeated_vector <- rep(ID.vcf , each = 2)
# Create a new vector with .hap1 and .hap2 appended
new_vector <- mapply(function(x, y) paste0(x, ".hap", y), repeated_vector, rep(1:2, length(ID.vcf )))
# Flatten the vector since mapply returns a matrix
ID.haps <- as.vector(new_vector)


# Apply the split function
split_haps <- split_haplotypes(haplotypes)

haps <- data.frame(t(split_haps))

# Concatenate all columns into one column per row
comb_haps <- apply(haps , 1, function(row) paste(row, collapse = ""))


# Function to split by the last underscore
split_by_last_underscore <- function(x) {
  s <- sub("^(.*)_(.*)$", "\\1", x)
  e <- sub("^(.*)_(.*)$", "\\2", x)
  return(data.frame(start = s, end = e))
}

# Apply the function to the column
split_columns <- do.call(rbind, lapply(ID.haps, split_by_last_underscore))
comb_haps2 <- data.frame(split_columns ,comb_haps)
# Write the data to a file
write.table(comb_haps2 , file = "Holstein.formatted_haplotypes.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
