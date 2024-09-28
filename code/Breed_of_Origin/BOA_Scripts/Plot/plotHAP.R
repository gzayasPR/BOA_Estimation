library(reshape2)
library(tidyverse)
library(ggplot2)
# Get working directory and save into vector
path <- getwd()
# Read pattern from working directory and save file names into a vector
Chr.snp <- dir(path, pattern ="chr.*.pos$")
file.names <- dir(path, pattern ="*.BO$")
#create a blank vector for barplots and line plots
barplots <- c()
Lplots <- c()
# Makes a line plot and barplot for each file in file.names and saves as an object into barplots and Lplots
for (i in 01:length(file.names)){
  x <- read.table(file.names[i])
  x.mean <- rbind( "Angus" = colMeans(x), "Brahman"= 1 - colMeans(x))
  y <- melt(x.mean, measure.vars = c("Brahman","Angus"))
  snp <- read.table(Chr.snp[i], header = F) %>% .[rep(1:nrow(.),each=2),]/1000000
  x.area <- cbind(y[,c(1,3)],snp)
  colnames(x.area) <- c("Breed","Percent","Position")
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1)
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  p <- ggplot(x.area) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = F),alpha=.8)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=.5,size = 3)+
    ylim(0,1.01)+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  barplots[[Chr.NO]] <- p
  #Linear plots
  Linearplot_FN <- paste("Chromosome",Chr.NO,"Linearplot.png")
  # Linear plots
  x.line <- cbind( "Angus" = colMeans(x), "Brahman"= 1 - colMeans(x),"Position" = read.table(Chr.snp[i], header = F)/1000000)

  colnames(x.line) <- c("Angus","Brahman","Position")
  l <- ggplot(x.line) + geom_line(aes(y=Brahman,x=Position), col="orange") +
    geom_line(aes(y=Angus,x=Position), col="blue") + geom_hline(yintercept=.625,size = 2)+
    ylab("Haplotype Origin") + xlab(paste("Chromosome",Chr.NO,"Position(Mb)"))+ scale_fill_manual(values = c("orange","blue"))
  Lplots[[Chr.NO]] <- l
}

#plot chromosomes together and save, # rename photos to what you want
# barplots
name.txt <- read.table("NAME.txt")
png1_15.Bar <- paste(name.txt,"Chromosomes 1-15 Haplotypes Barplots.png",sep= " ")
png15_29.Bar  <- paste(name.txt,"Chromosomes 16-29 Haplotypes Barplots.png",sep= " ")
library(ggpubr)
png(png1_15.Bar ,width=1500,height=900)
print(ggarrange(barplots$'01',barplots$'02',barplots$'03',barplots$'04',barplots$'05',barplots$'06',barplots$'07',barplots$'08',barplots$'09',barplots$'10',barplots$'11',barplots$'12',barplots$'13',barplots$'14',barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png(png15_29.Bar , width=1500, height=900)
print(ggarrange(barplots$'16',barplots$'17',barplots$'18',barplots$'19',barplots$'20',barplots$'21',barplots$'22',barplots$'23',barplots$'24',barplots$'25',barplots$'26',barplots$'27',barplots$'28',barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
# linear plots
png1_15.line <- paste(name.txt,"Chromosomes 1-15 Haplotypes Linearplot.png",sep= " ")
png15_29.line  <- paste(name.txt,"Chromosomes 16-29 Haplotypes Linearplot.png",sep= " ")
png(png1_15.line ,width=1500,height=900)
print(ggarrange(Lplots$'01',Lplots$'02',Lplots$'03',Lplots$'04',Lplots$'05',Lplots$'06',Lplots$'07',Lplots$'08',Lplots$'09',Lplots$'10',Lplots$'11',Lplots$'12',Lplots$'13',Lplots$'14',Lplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()
png(png15_29.line, width=1500, height=900)
print(ggarrange(Lplots$'16',Lplots$'17',Lplots$'18',Lplots$'19',Lplots$'20',Lplots$'21',Lplots$'22',Lplots$'23',Lplots$'24',Lplots$'25',Lplots$'26',Lplots$'27',Lplots$'28',Lplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()
