library(reshape2)
library(tidyverse)
library(ggplot2)
# Get working directory and save into vector
path <- getwd()
# Read pattern from working directory and save file names into a vector
frq.names <- dir(path, pattern ="*.frq$")
freq.names<- dir(path, pattern ="*.frqx")
map.files <- dir(path, pattern ="*.bim")
fam.files <- dir(path, pattern ="*.fam")
#create a blank vector for barplots and line plots
name <- read.table("NAME.txt")
name
frq.df <- read.table(freq.names, header=T,sep="\t")
map <- read.table(map.files, header=F)[,c(1,4,2)]
names(map) <- c("CHR","POS","SNP")
merge.df <- merge(frq.df,map,by=c("SNP","CHR"))
n <- nrow(read.table(fam.files , header=F))
merge.df2 <- merge.df[,c(1,2,11,3,4,5,6,7)]
merge.df2[6:8] <- merge.df2[6:8]/n
names(merge.df2)[1:3] <- c("SNP","CHR","BP")
ANGUS.A2 <- merge.df2[merge.df2$A2 == "A",] 
names(ANGUS.A2) <- c("SNP","CHR","BP","BRH","AN","HOM_BB","HET","HOM_AA")
BRAHMAN.A2 <- merge.df2[merge.df2$A2 == "B",]
names(BRAHMAN.A2) <- c("SNP","CHR","BP","AN","BRH","HOM_AA","HET","HOM_BB")
DF2 <- rbind(ANGUS.A2,BRAHMAN.A2)
data_cum <- DF2  %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP)) %>%
  mutate(BP_add = lag(cumsum(as.numeric(max_BP)), default = 0)) %>%
  select(CHR, BP_add)

DF2 <- DF2 %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(BP_cum = BP + BP_add)

length(DF2$BP_cum)
length(unique(DF2$BP_cum))

DF3 <- DF2[!duplicated(DF2[,c("BP_cum")]),]
gwas_data <- melt(DF3, id.vars =c("SNP","CHR","BP","AN","BRH","BP_cum","BP_add"))
#gwas_data[is.na(gwas_data)] <- 0
data_cum <- gwas_data %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP)) %>%
  mutate(BP_add = lag(cumsum(as.numeric(max_BP)), default = 0)) %>%
  select(CHR, BP_add)

#gwas_data <- gwas_data %>%
#  inner_join(data_cum, by = "CHR") %>%
#  mutate(BP_cum = BP + BP_add)

axis_set <- gwas_data %>%
  group_by(CHR) %>%
  summarize(center = mean(BP_cum))


gwas_data$variable <- factor(gwas_data$variable, levels = c("HOM_AA","HET","HOM_BB"))
#head(gwas_data)
#tail(gwas_data)
df <- gwas_data
#df = df[order(df[,'Date'],-df[,'Depth']),]
#df = df[!duplicated(df$BP_cum),]
#gwas_data <- gwas_data[unique(gwas_data$BP_cum),]
names(df)[names(df) == "variable"]  <- c("Breed of Origin Genotype")
manhplot <- ggplot(df, aes(x=BP_cum,y = value,color=`Breed of Origin Genotype`
                           ,fill=`Breed of Origin Genotype`)) +
  geom_area(stat="identity", size=0.5) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)  +
  ylab("Breed of Origin Frequency") + xlab("Chromosome")+ ggtitle("Breed of Origin Frequency") + 
  scale_color_manual(values= c("orange","green","blue"))+ 
  scale_fill_manual(values= c("orange","green","blue")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y=element_text(size=12,face="bold"),
        legend.key.size = unit(2, 'cm'),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,face="italic"))
#manhplot
ggsave(manhplot,file=paste(name,"Genotype.BO.jpg",sep="."),height = 16, width=25)

df1_14 <- df[df$CHR %in% 1:14,]

manhplot <- ggplot(df, aes(x=BP,y = value,color=`Breed of Origin Genotype`
                           ,fill=`Breed of Origin Genotype`)) +
  geom_area(stat="identity", size=0.5)  + facet_wrap(as.factor(CHR) ~ ., scales='free', ncol=5) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)  +
  ylab("Breed of Origin Frequency") + xlab("Chromosome")+ ggtitle("Breed of Origin Frequency") + 
  scale_color_manual(values= c("orange","green","blue"))+ 
  scale_fill_manual(values= c("orange","green","blue")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y=element_text(size=12,face="bold"),
        legend.key.size = unit(2, 'cm'),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,face="italic"))
#manhplot
ggsave(manhplot,file=paste(name,"Genotype_1_14.BO.jpg",sep="."),height = 16, width=25)

df15_29 <- df[df$CHR %in% 15:29,]
manhplot <- ggplot(df15_29, aes(x=BP,y = value,color=`Breed of Origin Genotype`
                           ,fill=`Breed of Origin Genotype`)) +
  geom_area(stat="identity", size=0.5)  + facet_wrap(as.factor(CHR) ~ ., scales='free', ncol=5) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)  +
  ylab("Breed of Origin Frequency") + xlab("Chromosome")+ ggtitle("Breed of Origin Frequency") + 
  scale_color_manual(values= c("orange","green","blue"))+ 
  scale_fill_manual(values= c("orange","green","blue")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y=element_text(size=12,face="bold"),
        legend.key.size = unit(2, 'cm'),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,face="italic"))
#manhplot
ggsave(manhplot,file=paste(name,"Genotype_15_29.BO.jpg",sep="."),height = 16, width=25)



frq.df <- read.table(frq.names, header=T)
map <- read.table(map.files, header=F)[,c(1,4,2)]
names(map) <- c("CHR","POS","SNP")
merge.df <- merge(frq.df,map,by=c("SNP","CHR"))
n <- nrow(read.table(fam.files , header=F))
merge.df$MJAF <- (1 - merge.df$MAF)
merge.df2 <- merge.df[,c(1,2,7,3,4,5,8)]
names(merge.df2)[1:3] <- c("SNP","CHR","BP")
ANGUS.A2 <- merge.df2[merge.df2$A2 == "A",] 
names(ANGUS.A2) <- c("SNP","CHR","BP","BRH","AN","Brahman","Angus")
BRAHMAN.A2 <- merge.df2[merge.df2$A2 == "B",]
names(BRAHMAN.A2) <- c("SNP","CHR","BP","AN","BRH","Angus","Brahman")
DF2 <- rbind(ANGUS.A2,BRAHMAN.A2)

data_cum <- DF2  %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP)) %>%
  mutate(BP_add = lag(cumsum(as.numeric(max_BP)), default = 0)) %>%
  select(CHR, BP_add)

DF2 <- DF2 %>%
  inner_join(data_cum, by = "CHR") %>%
  mutate(BP_cum = BP + BP_add)

length(DF2$BP_cum)
length(unique(DF2$BP_cum))


DF3 <- DF2[!duplicated(DF2[,c("BP_cum")]),]
gwas_data <- melt(DF3, id.vars =c("SNP","CHR","BP","AN","BRH","BP_cum","BP_add"))
#gwas_data[is.na(gwas_data)] <- 0
data_cum <- gwas_data %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP)) %>%
  mutate(BP_add = lag(cumsum(as.numeric(max_BP)), default = 0)) %>%
  select(CHR, BP_add)

#gwas_data <- gwas_data %>%
#  inner_join(data_cum, by = "CHR") %>%
#  mutate(BP_cum = BP + BP_add)

axis_set <- gwas_data %>%
  group_by(CHR) %>%
  summarize(center = mean(BP_cum))


gwas_data$variable <- factor(gwas_data$variable, levels = c("Angus","Brahman"))
#head(gwas_data)
#tail(gwas_data)
df <- gwas_data
#df = df[order(df[,'Date'],-df[,'Depth']),]
#df = df[!duplicated(df$BP_cum),]
#gwas_data <- gwas_data[unique(gwas_data$BP_cum),]
names(df)[names(df) == "variable"]  <- c("Breed of Origin Gene Frequencies")
manhplot <- ggplot(df, aes(x=BP_cum,y = value,color=`Breed of Origin Gene Frequencies`
                           ,fill=`Breed of Origin Gene Frequencies`)) +
  geom_area(stat="identity", size=0.5) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)  +
  ylab("Breed of Origin Gene Frequency") + xlab("Chromosome")+ ggtitle("Breed of Origin Gene Frequency") + 
  scale_color_manual(values= c("orange","blue"))+ 
  scale_fill_manual(values= c("orange","blue")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y=element_text(size=12,face="bold"),
        legend.key.size = unit(2, 'cm'),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,face="italic"))
#manhplot
ggsave(manhplot,file=paste(name,"Allele.BO.jpg",sep="."),height = 16, width=25)
df1_14 <- df[df$CHR %in% 1:14,]

manhplot1_14 <- ggplot(df1_14, aes(x=BP,y = value,color=`Breed of Origin Gene Frequencies`
                           ,fill=`Breed of Origin Gene Frequencies`)) +
  geom_area(stat="identity", size=0.5) + facet_wrap(as.factor(CHR) ~ ., scales='free', ncol=5) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)  +
  ylab("Breed of Origin Gene Frequency") + xlab("Chromosome")+ ggtitle("Breed of Origin Gene Frequency") + 
  scale_color_manual(values= c("orange","blue"))+ 
  scale_fill_manual(values= c("orange","blue")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=13,face="bold"),
        axis.text.y=element_text(size=10,face="bold"),
        legend.key.size = unit(2, 'cm'),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,face="italic"))
#manhplot1_14
ggsave(manhplot1_14,file=paste(name,"Allele.1_14.BO.jpg",sep="."),height = 16, width=25)

df15_29 <- df[df$CHR %in% 15:29,]
manhplot15_29 <- ggplot(df15_29, aes(x=BP,y = value,color=`Breed of Origin Gene Frequencies`
                           ,fill=`Breed of Origin Gene Frequencies`)) +
  geom_area(stat="identity", size=0.5) + facet_wrap(as.factor(CHR) ~ ., scales='free', ncol=5) +
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center)  +
  ylab("Breed of Origin Gene Frequency") + xlab("Chromosome")+ ggtitle("Breed of Origin Gene Frequency") + 
  scale_color_manual(values= c("orange","blue"))+ 
  scale_fill_manual(values= c("orange","blue")) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=15,face="bold"),
        axis.text.y=element_text(size=12,face="bold"),
        legend.key.size = unit(2, 'cm'),
        plot.title=element_text(hjust=0.5,face="bold"),
        plot.subtitle=element_text(hjust=0.5,face="italic"))
#manhplot15_29
ggsave(manhplot15_29,file=paste(name,"Allele.15_29.BO.jpg",sep="."),height = 16, width=25)