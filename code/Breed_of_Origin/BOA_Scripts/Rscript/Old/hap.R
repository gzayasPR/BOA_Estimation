  hap.2 <- t(read.table("hap.2"))
  hap.2.1 <- hap.2[-1,]
  hap.2.2 <- hap.2.1[-1,]
  write.table(hap.2.2,"hap.3",row.names=F, col.names=F, quote =F)