

raw.files <- list.files(pattern="*.raw$")

raw_reads <- c()
for (Window in raw.files){
    raw_reads[[Window]] <- read.table(Window, header=T)[,-c(1,3,4,5,6)]
}

concordance <- c()
for (Window1 in raw.files){
    i <- match(Window1,raw.files)
    for (Window2 in raw.files[-i]){
    j <- match(Window2,raw.files)   
    if (j > i) {
    windows_compare <- paste( Window1,Window2,sep="_vs_")  
    comparison <- raw_reads[[raw.files[i]]][-1] == raw_reads[[raw.files[j]]][-1]
    correct <- sum(comparison )
    total <- length(comparison )
    concordance[[windows_compare]] <- correct/total
    }else{}}}

concordance_df <- data.frame(
    Pairwise_Comparison = names(concordance),
    Concordance_Score = unlist(concordance))
print(concordance_df)
write.table(concordance_df,"Concordance.txt")