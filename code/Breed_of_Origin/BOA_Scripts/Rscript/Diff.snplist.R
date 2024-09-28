# Clear the workspace
rm(list = ls())

# Load arguments from the command line
args <- commandArgs(TRUE)
afd <- args[1]

library(tidyverse)
dataframe <- read.table("Pure.frq.strat", header = T)
df2 <- dataframe %>% group_by(CHR,SNP) %>% summarise(dif= abs(MAF[2] - MAF[1]) )
diff.snp <- df2[df2$dif > as.numeric(afd),]
write.table(diff.snp[2],"DIFF.snplist",col.names=F,row.name=F,quote=F)
