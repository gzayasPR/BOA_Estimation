X.hap <- cbind(read.table("hap.ID"),read.table("hap.ID"),read.table("hap.6"),sep=' ')
write.table(X.hap,"hap.7",row.names=F,col.names=F,quote =F)
