#setwd("//ad.ufl.edu/ifas/ANS/Groups/MateescuLab/Mutibreed-Breed Angus-Brahman Project/Results/BreedComposition/Seminole_Haplotypes")
#getwd()
library(ggplot2)
path <- getwd()
out.file<-""

Chr.snp <- dir(path, pattern ="chr.*.pos$")
file.names <- dir(path, pattern ="lampld.std_ancestry.*.txt$")
barplots <- c()
Lplots <- c()
for (i in 1:length(file.names)){
  lamp <- read.table(file.names[i],  header=F)  
  lamp.m <- colMeans(lamp[, 2:ncol(lamp)])
  tmp <- data.frame(X=lamp.m, ind=rep(1:2))
  #lamp.df1 <- unstack(tmp, X~ind)
  snp <- read.table(Chr.snp[i], header = F)
  snp2 <- (snp[rep(1:nrow(snp),each=2),]/1000000)
  plot_matrix <- data.frame(tmp[1]*100,tmp[2],snp2)
  head(plot_matrix)
  plot_matrix2 <- data.frame(plot_matrix[,1],plot_matrix[,3],sapply(plot_matrix$ind,switch,'1'='Angus','2'='Brahman'))
  colnames(plot_matrix2) <- c("Percent","Position","Breed")
  head(plot_matrix2)
  Chr.NO1 <- sub("chr","",Chr.snp[i])
  Chr.NO <- sub(".pos","",Chr.NO1) 
  #Barplots
  Barplot_FN <- paste("Chromosome",Chr.NO,"Barplot.png")
  Barplot_FN
  p <- ggplot(plot_matrix2) + geom_area(aes(y = Percent, x = Position,fill=Breed),position = position_stack(vjust = 0.5,reverse = T),alpha=.6)+
    ylab("Haplotype Origin") +   xlab(paste("Chromosome",Chr.NO,"Position(Mb)")) + geom_hline(yintercept=62.5,size = 3)+ 
    coord_cartesian(ylim= c(0,101))+ theme_classic() + scale_fill_manual(values = c("blue","orange"))
  barplots[[Chr.NO]] <- p
  #Linear plots
  Linearplot_FN <- paste("Chromosome",Chr.NO,"Linearplot.png")
  # Linear plots
  lamp.df1 <- unstack(tmp, X~ind)
  linear_matrix <- data.frame(lamp.df1,Position = snp/1000000)
  names(linear_matrix) <- c("Angus","Brahman","Position")
  head(linear_matrix)
  l <- ggplot(linear_matrix) + geom_line(aes(y=Angus,x=Position), col="orange") + 
    geom_line(aes(y=Brahman,x=Position), col="blue") + geom_hline(yintercept=.625,size = 2)+ 
    ylab("Haplotype Origin") + xlab(paste("Chromosome",Chr.NO,"Position(Mb)"))+ scale_fill_manual(values = c("orange","blue"))
  Lplots[[Chr.NO]] <- l
}

#plots
#install.packages('ggpubr')
library(ggpubr)
png("Thermotolerance Chromosomes 1-15 Haplotypes Barplots.png",width=1500,height=900)
print(ggarrange(barplots$'01',barplots$'02',barplots$'03',barplots$'04',barplots$'05',barplots$'06',barplots$'07',barplots$'08',barplots$'09',barplots$'10',barplots$'11',barplots$'12',barplots$'13',barplots$'14',barplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()

png("Thermotolerance Chromosmes 16-29 Haplotypes Barplots.png", width=1500, height=900)
print(ggarrange(barplots$'16',barplots$'17',barplots$'18',barplots$'19',barplots$'20',barplots$'21',barplots$'22',barplots$'23',barplots$'24',barplots$'25',barplots$'26',barplots$'27',barplots$'28',barplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()


png("Thermotolerance Chromosomes 1-15 Haplotypes Linearplot.png",width=1500,height=900)
print(ggarrange(Lplots$'01',Lplots$'02',Lplots$'03',Lplots$'04',Lplots$'05',Lplots$'06',Lplots$'07',Lplots$'08',Lplots$'09',Lplots$'10',Lplots$'11',Lplots$'12',Lplots$'13',Lplots$'14',Lplots$'15',common.legend=TRUE,labels = c(1:15),ncol = 5, nrow = 3))
dev.off()

png("Thermotolerance Chromosmes 16-29 Haplotypes  Linearplot.png", width=1500, height=900)
print(ggarrange(Lplots$'16',Lplots$'17',Lplots$'18',Lplots$'19',Lplots$'20',Lplots$'21',Lplots$'22',Lplots$'23',Lplots$'24',Lplots$'25',Lplots$'26',Lplots$'27',Lplots$'28',Lplots$'29',common.legend=TRUE,labels = c(16:29),ncol = 5, nrow = 3))
dev.off()





