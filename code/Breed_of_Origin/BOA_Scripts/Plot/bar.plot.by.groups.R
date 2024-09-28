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
#create a blank vector for barplots and line plots
G1.barplots <- c()
G2.barplots <- c()
G3.barplots <- c()
G4.barplots <- c()
G5.barplots <- c()
G6.barplots <- c()
Lplots <- c()
# Makes a line plot and barplot for each file in file.names and saves as an object into barplots and Lplots
for (i in 01:length(file.names)){
  x <- cbind(read.table("MAB.ID",col.names="IID"),read.table(file.names[i]))
  G1.x <- merge(x,G1,by="IID")         
  G1.x.mean <- rbind( "Brahman" = colMeans(G1.x[-1]), "Angus"= 1 - colMeans(G1.x[-1]))
  G1.y <- melt(G1.x.mean, measure.vars = c("Angus","Brahman"))
  snp <- read.table(Chr.snp[i], header = F) %>% .[rep(1:nrow(.),each=2),]/1000000
  x.area <- cbind(G1.y[,c(1,3)],snp)
  colnames(x.area) <- c("Breed","Percent","Position")
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  G1.p <- ggplot(x.area) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = F),alpha=.8)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=.625,size = 3)+ 
    ylim(0,1.01)+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  G1.barplots[[Chr.NO]] <- G1.p
  G2.x <- merge(x,G2,by="IID")         
  G2.x.mean <- rbind( "Brahman" = colMeans(G2.x[-1]), "Angus"= 1 - colMeans(G2.x[-1]))
  G2.y <- melt(G2.x.mean, measure.vars = c("Angus","Brahman"))
  snp <- read.table(Chr.snp[i], header = F) %>% .[rep(1:nrow(.),each=2),]/1000000
  x.area <- cbind(G2.y[,c(1,3)],snp)
  colnames(x.area) <- c("Breed","Percent","Position")
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  G2.p <- ggplot(x.area) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = F),alpha=.8)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=.625,size = 3)+ 
    ylim(0,1.01)+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  G2.barplots[[Chr.NO]] <- G2.p
  G3.x <- merge(x,G3,by="IID")         
  G3.x.mean <- rbind( "Brahman" = colMeans(G3.x[-1]), "Angus"= 1 - colMeans(G3.x[-1]))
  G3.y <- melt(G3.x.mean, measure.vars = c("Angus","Brahman"))
  snp <- read.table(Chr.snp[i], header = F) %>% .[rep(1:nrow(.),each=2),]/1000000
  x.area <- cbind(G3.y[,c(1,3)],snp)
  colnames(x.area) <- c("Breed","Percent","Position")
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  G3.p <- ggplot(x.area) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = F),alpha=.8)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=.625,size = 3)+ 
    ylim(0,1.01)+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  G3.barplots[[Chr.NO]] <- G3.p
  G4.x <- merge(x,G4,by="IID")         
  G4.x.mean <- rbind( "Brahman" = colMeans(G4.x[-1]), "Angus"= 1 - colMeans(G4.x[-1]))
  G4.y <- melt(G4.x.mean, measure.vars = c("Angus","Brahman"))
  snp <- read.table(Chr.snp[i], header = F) %>% .[rep(1:nrow(.),each=2),]/1000000
  x.area <- cbind(G4.y[,c(1,3)],snp)
  colnames(x.area) <- c("Breed","Percent","Position")
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  G4.p <- ggplot(x.area) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = F),alpha=.8)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=.625,size = 3)+ 
    ylim(0,1.01)+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  G4.barplots[[Chr.NO]] <- G4.p
  G5.x <- merge(x,G5,by="IID")         
  G5.x.mean <- rbind( "Brahman" = colMeans(G5.x[-1]), "Angus"= 1 - colMeans(G5.x[-1]))
  G5.y <- melt(G5.x.mean, measure.vars = c("Angus","Brahman"))
  snp <- read.table(Chr.snp[i], header = F) %>% .[rep(1:nrow(.),each=2),]/1000000
  x.area <- cbind(G5.y[,c(1,3)],snp)
  colnames(x.area) <- c("Breed","Percent","Position")
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  G5.p <- ggplot(x.area) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = F),alpha=.8)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=.625,size = 3)+ 
    ylim(0,1.01)+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  G5.barplots[[Chr.NO]] <- G5.p
  G6.x <- merge(x,G6,by="IID")         
  G6.x.mean <- rbind( "Brahman" = colMeans(G6.x[-1]), "Angus"= 1 - colMeans(G6.x[-1]))
  G6.y <- melt(G6.x.mean, measure.vars = c("Angus","Brahman"))
  snp <- read.table(Chr.snp[i], header = F) %>% .[rep(1:nrow(.),each=2),]/1000000
  x.area <- cbind(G6.y[,c(1,3)],snp)
  colnames(x.area) <- c("Breed","Percent","Position")
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  G6.p <- ggplot(x.area) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = F),alpha=.8)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=.625,size = 3)+ 
    ylim(0,1.01)+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  G6.barplots[[Chr.NO]] <- G6.p
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

#plot chromosomes together and save, # rename photos to what you want
# G1.barplots
png("MAB Angus.G1 Chromosomes 1-15 Haplotypes G1.barplots.png",width=1500,height=900)
print(ggarrange(G1.barplots$'01',G1.barplots$'02',G1.barplots$'03',G1.barplots$'04',G1.barplots$'05',G1.barplots$'06',G1.barplots$'07',G1.barplots$'08',G1.barplots$'09',G1.barplots$'10',G1.barplots$'11',G1.barplots$'12',G1.barplots$'13',G1.barplots$'14',G1.barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB Angus.G1 Chromosmes 16-29 Haplotypes G1.barplots.png", width=1500, height=900)
print(ggarrange(G1.barplots$'16',G1.barplots$'17',G1.barplots$'18',G1.barplots$'19',G1.barplots$'20',G1.barplots$'21',G1.barplots$'22',G1.barplots$'23',G1.barplots$'24',G1.barplots$'25',G1.barplots$'26',G1.barplots$'27',G1.barplots$'28',G1.barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
# G2.barplots
png("MAB 75.G2 Chromosomes 1-15 Haplotypes G2.barplots.png",width=1500,height=900)
print(ggarrange(G2.barplots$'01',G2.barplots$'02',G2.barplots$'03',G2.barplots$'04',G2.barplots$'05',G2.barplots$'06',G2.barplots$'07',G2.barplots$'08',G2.barplots$'09',G2.barplots$'10',G2.barplots$'11',G2.barplots$'12',G2.barplots$'13',G2.barplots$'14',G2.barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB 75.G2 Chromosmes 16-29 Haplotypes G2.barplots.png", width=1500, height=900)
print(ggarrange(G2.barplots$'16',G2.barplots$'17',G2.barplots$'18',G2.barplots$'19',G2.barplots$'20',G2.barplots$'21',G2.barplots$'22',G2.barplots$'23',G2.barplots$'24',G2.barplots$'25',G2.barplots$'26',G2.barplots$'27',G2.barplots$'28',G2.barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
# G3.barplots
png("MAB Brangus.G3 Chromosomes 1-15 Haplotypes G3.barplots.png",width=1500,height=900)
print(ggarrange(G3.barplots$'01',G3.barplots$'02',G3.barplots$'03',G3.barplots$'04',G3.barplots$'05',G3.barplots$'06',G3.barplots$'07',G3.barplots$'08',G3.barplots$'09',G3.barplots$'10',G3.barplots$'11',G3.barplots$'12',G3.barplots$'13',G3.barplots$'14',G3.barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB Brangus.G3 Chromosmes 16-29 Haplotypes G3.barplots.png", width=1500, height=900)
print(ggarrange(G3.barplots$'16',G3.barplots$'17',G3.barplots$'18',G3.barplots$'19',G3.barplots$'20',G3.barplots$'21',G3.barplots$'22',G3.barplots$'23',G3.barplots$'24',G3.barplots$'25',G3.barplots$'26',G3.barplots$'27',G3.barplots$'28',G3.barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
# G4.barplots
png("MAB 50.G4 Chromosomes 1-15 Haplotypes G4.barplots.png",width=1500,height=900)
print(ggarrange(G4.barplots$'01',G4.barplots$'02',G4.barplots$'03',G4.barplots$'04',G4.barplots$'05',G4.barplots$'06',G4.barplots$'07',G4.barplots$'08',G4.barplots$'09',G4.barplots$'10',G4.barplots$'11',G4.barplots$'12',G4.barplots$'13',G4.barplots$'14',G4.barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB 50.G4 Chromosmes 16-29 Haplotypes G4.barplots.png", width=1500, height=900)
print(ggarrange(G4.barplots$'16',G4.barplots$'17',G4.barplots$'18',G4.barplots$'19',G4.barplots$'20',G4.barplots$'21',G4.barplots$'22',G4.barplots$'23',G4.barplots$'24',G4.barplots$'25',G4.barplots$'26',G4.barplots$'27',G4.barplots$'28',G4.barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
# G5.barplots
png("MAB 25.G5 Chromosomes 1-15 Haplotypes G5.barplots.png",width=1500,height=900)
print(ggarrange(G5.barplots$'01',G5.barplots$'02',G5.barplots$'03',G5.barplots$'04',G5.barplots$'05',G5.barplots$'06',G5.barplots$'07',G5.barplots$'08',G5.barplots$'09',G5.barplots$'10',G5.barplots$'11',G5.barplots$'12',G5.barplots$'13',G5.barplots$'14',G5.barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB 25.G5 Chromosmes 16-29 Haplotypes G5.barplots.png", width=1500, height=900)
print(ggarrange(G5.barplots$'16',G5.barplots$'17',G5.barplots$'18',G5.barplots$'19',G5.barplots$'20',G5.barplots$'21',G5.barplots$'22',G5.barplots$'23',G5.barplots$'24',G5.barplots$'25',G5.barplots$'26',G5.barplots$'27',G5.barplots$'28',G5.barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
# G6.barplots
png("MAB Brahman.G6 Chromosomes 1-15 Haplotypes G6.barplots.png",width=1500,height=900)
print(ggarrange(G6.barplots$'01',G6.barplots$'02',G6.barplots$'03',G6.barplots$'04',G6.barplots$'05',G6.barplots$'06',G6.barplots$'07',G6.barplots$'08',G6.barplots$'09',G6.barplots$'10',G6.barplots$'11',G6.barplots$'12',G6.barplots$'13',G6.barplots$'14',G6.barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB Brahman.G6 Chromosmes 16-29 Haplotypes G6.barplots.png", width=1500, height=900)
print(ggarrange(G6.barplots$'16',G6.barplots$'17',G6.barplots$'18',G6.barplots$'19',G6.barplots$'20',G6.barplots$'21',G6.barplots$'22',G6.barplots$'23',G6.barplots$'24',G6.barplots$'25',G6.barplots$'26',G6.barplots$'27',G6.barplots$'28',G6.barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()

# Lplots
png("MAB Chromosomes 1-15 Haplotypes Lplots.png",width=1500,height=900)
print(ggarrange(Lplots$'01',Lplots$'02',Lplots$'03',Lplots$'04',Lplots$'05',Lplots$'06',Lplots$'07',Lplots$'08',Lplots$'09',Lplots$'10',Lplots$'11',Lplots$'12',Lplots$'13',Lplots$'14',Lplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png("MAB Chromosmes 16-29 Haplotypes Lplots.png", width=1500, height=900)
print(ggarrange(Lplots$'16',Lplots$'17',Lplots$'18',Lplots$'19',Lplots$'20',Lplots$'21',Lplots$'22',Lplots$'23',Lplots$'24',Lplots$'25',Lplots$'26',Lplots$'27',Lplots$'28',Lplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
