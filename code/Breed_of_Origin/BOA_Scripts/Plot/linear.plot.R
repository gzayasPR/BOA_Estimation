library(reshape2)
library(tidyverse)
library(ggplot2)
library(xlsx)
library(ggpubr)
# Get working directory and save into vector
path <- getwd()
# Read pattern from working directory and save file names into a vector
Chr.snp <- dir(path, pattern ="chr.*.pos$")
file.names <- dir(path, pattern ="*.BO$")
G1 <- read.xlsx2("Grouped.MAB.xlsx",sheetIndex = 1)[2]
G2 <- read.xlsx2("Grouped.MAB.xlsx",sheetIndex = 2)[2]
G3 <- read.xlsx2("Grouped.MAB.xlsx",sheetIndex = 3)[2]
G4 <- read.xlsx2("Grouped.MAB.xlsx",sheetIndex = 4)[2]
G5 <- read.xlsx2("Grouped.MAB.xlsx",sheetIndex = 5)[2]
G6 <- read.xlsx2("Grouped.MAB.xlsx",sheetIndex = 6)[2]
Lplots <- c()
for (i in 01:length(file.names)){
  x <- cbind(read.table("MAB.ID",col.names="IID"),read.table(file.names[i]))
  G1.x <- merge(x,G1,by="IID")         
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  G2.x <- merge(x,G2,by="IID")         
  G3.x <- merge(x,G3,by="IID")         
  G4.x <- merge(x,G4,by="IID")         
  G5.x <- merge(x,G5,by="IID")         
  G6.x <- merge(x,G6,by="IID")         
  Linearplot_FN <- paste("Chromosome",Chr.NO,"Linearplot.png")
  # Linear plots
  G1.x.line <- cbind( "Brahman" = colMeans(G1.x[-1]), "Angus"= 1 - colMeans(G1.x[-1]),"Position" = read.table(Chr.snp[i], header = F)/1000000)
  colnames(G1.x.line) <- c("Brahman","Angus","Position")
  G2.x.line <- cbind( "Brahman" = colMeans(G2.x[-1]), "Angus"= 1 - colMeans(G2.x[-1]),"Position" = read.table(Chr.snp[i], header = F)/1000000)
  colnames(G2.x.line) <- c("Brahman","Angus","Position")
  G3.x.line <- cbind( "Brahman" = colMeans(G3.x[-1]), "Angus"= 1 - colMeans(G3.x[-1]),"Position" = read.table(Chr.snp[i], header = F)/1000000)
  colnames(G3.x.line) <- c("Brahman","Angus","Position")
  G4.x.line <- cbind( "Brahman" = colMeans(G4.x[-1]), "Angus"= 1 - colMeans(G4.x[-1]),"Position" = read.table(Chr.snp[i], header = F)/1000000)
  colnames(G4.x.line) <- c("Brahman","Angus","Position")
  G5.x.line <- cbind( "Brahman" = colMeans(G5.x[-1]), "Angus"= 1 - colMeans(G5.x[-1]),"Position" = read.table(Chr.snp[i], header = F)/1000000)
  colnames(G5.x.line) <- c("Brahman","Angus","Position")
  G6.x.line <- cbind( "Brahman" = colMeans(G6.x[-1]), "Angus"= 1 - colMeans(G6.x[-1]),"Position" = read.table(Chr.snp[i], header = F)/1000000)
  colnames(G6.x.line) <- c("Brahman","Angus","Position")
  Groups.line <- cbind(G1.x.line,G2.x.line,G3.x.line,G4.x.line,G5.x.line,G6.x.line)
  Groups.line <- Groups.line[,c(3,1,2,4,5,7,8,10,11,13,14,16,17)]
  names(Groups.line) <- c("Position","Brahman1","Angus1","Brahman2","Angus2","Brahman3","Angus3","Brahman4","Angus4",
                          "Brahman5","Angus5","Brahman6","Angus6")
  l <- ggplot(Groups.line) + geom_line(aes(y=Angus1,x=Position,color="orange")) + 
    geom_line(aes(y=Angus2,x=Position,color="red"))+
    geom_line(aes(y=Angus3,x=Position,color="yellow")) +
    geom_line(aes(y=Angus4,x=Position,color="purple")) +
    geom_line(aes(y=Angus5,x=Position,color="green"))+
    geom_line(aes(y=Angus6,x=Position,color="blue"))+
    geom_hline(yintercept=.625,size = 2)+ 
    ylab("Angus percent Haplotype Origin") + scale_color_manual(name = "Markers per Haplotype ",
                                                                values= c("orange","red","yellow","purple","green","blue"),
                                                                labels = c("Brahman","1/4 Angus","Angus","1/2 Angus","3/4 Angus","Brangus")) + 
    xlab(paste(Chr.NO,"Position(Mb)"))
  Lplots[[Chr.NO]] <- l
}

# Lplots
png("MAB Chromosomes 1-15 Haplotypes Lplots.png",width=1500,height=900)
print(ggarrange(Lplots$'01',Lplots$'02',Lplots$'03',Lplots$'04',Lplots$'05',Lplots$'06',Lplots$'07',Lplots$'08',Lplots$'09',Lplots$'10',Lplots$'11',Lplots$'12',Lplots$'13',Lplots$'14',Lplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB Chromosmes 16-29 Haplotypes Lplots.png", width=1500, height=900)
print(ggarrange(Lplots$'16',Lplots$'17',Lplots$'18',Lplots$'19',Lplots$'20',Lplots$'21',Lplots$'22',Lplots$'23',Lplots$'24',Lplots$'25',Lplots$'26',Lplots$'27',Lplots$'28',Lplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
