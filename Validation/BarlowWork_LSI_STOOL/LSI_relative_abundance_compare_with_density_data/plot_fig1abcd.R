library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)
# library(extrafont)
# font_import()
# fonts()
# loadfonts()

##########################
##########################
##########################
##########################

##############################
##############################
##############################
##############################
##############################validation
##############################
##############################
##############################
#mutli disease


library(ggplot2)
theme_set(theme_bw())  
library(pheatmap)
library(RColorBrewer)


data<- read.table("LSI.csv", header = TRUE, sep = "," )
ketodata<-data[data$condition=='Keto',]
sampleid<-ketodata$sampleid
ketodata<-ketodata[ , !colnames(ketodata) %in% c("sampleid","condition")]
ketodata<-t(ketodata)
colnames(ketodata)<-sampleid

ketodata<-data.frame(ketodata)

data_group<-c(rep('k1',dim(ketodata)[1]),rep('k2',dim(ketodata)[1]),rep('k3',dim(ketodata)[1]),rep('k4',dim(ketodata)[1]),rep('k5',dim(ketodata)[1]),rep('k6',dim(ketodata)[1]))
abundance<-c(ketodata$k1,ketodata$k2,ketodata$k3,ketodata$k4,ketodata$k5,ketodata$k6)
data_bacteria<-c(rep(rownames(ketodata),dim(ketodata)[2]))
stacked_data<-data.frame(data_group,abundance,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria,levels=c('Faecalibaculum',
'Lactobacillus',
'Lactococcus',
'Akkermansia',
'Muribaculaceae__2',
'Romboutsia',
'others',
'Parasutterella',
'Enterococcus',
'A2',
'Muribaculaceae__',
'Enterorhabdus',
'Bacteroides')) 
stacked_data$data_group <- factor(stacked_data$data_group) 


colourCount = 13
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
g <- ggplot(stacked_data, aes(x=data_group,y=abundance,fill=data_bacteria))+  geom_bar(stat='identity') +   scale_fill_manual(values = getPalette(colourCount))

ggsave(g, filename='LSI_relative_abundance_keto.pdf', width=25, height=30, units=c("cm"))	


controldata<-data[data$condition=='Control',]
sampleid<-controldata$sampleid
controldata<-controldata[ , !colnames(controldata) %in% c("sampleid","condition")]
controldata<-t(controldata)
colnames(controldata)<-sampleid

controldata<-data.frame(controldata)

data_group<-c(rep('c1',dim(controldata)[1]),rep('c2',dim(controldata)[1]),rep('c3',dim(controldata)[1]),rep('c4',dim(controldata)[1]),rep('c5',dim(controldata)[1]))
abundance<-c(controldata$c1,controldata$c2,controldata$c3,controldata$c4,controldata$c5)
data_bacteria<-c(rep(rownames(controldata),dim(controldata)[2]))
stacked_data<-data.frame(data_group,abundance,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria,levels=c('Faecalibaculum',
'Lactobacillus',
'Lactococcus',
'Akkermansia',
'Muribaculaceae__2',
'Romboutsia',
'others',
'Parasutterella',
'Enterococcus',
'A2',
'Muribaculaceae__',
'Enterorhabdus',
'Bacteroides')) 
stacked_data$data_group <- factor(stacked_data$data_group) 


colourCount = 13
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
g <- ggplot(stacked_data, aes(x=data_group,y=abundance,fill=data_bacteria))+  geom_bar(stat='identity') +   scale_fill_manual(values = getPalette(colourCount))

ggsave(g, filename='LSI_relative_abundance_control.pdf', width=25, height=30, units=c("cm"))	



#########################
#########################
#########################
#########################


data<- read.table("LSI_density.csv", header = TRUE, sep = "," )
ketodata<-data[data$condition=='Keto',]
sampleid<-ketodata$sampleid
ketodata<-ketodata[ , !colnames(ketodata) %in% c("sampleid","condition")]
ketodata<-t(ketodata)
colnames(ketodata)<-sampleid

ketodata<-data.frame(ketodata)

data_group<-c(rep('k1',dim(ketodata)[1]),rep('k2',dim(ketodata)[1]),rep('k3',dim(ketodata)[1]),rep('k4',dim(ketodata)[1]),rep('k5',dim(ketodata)[1]),rep('k6',dim(ketodata)[1]))
abundance<-c(ketodata$k1,ketodata$k2,ketodata$k3,ketodata$k4,ketodata$k5,ketodata$k6)
data_bacteria<-c(rep(rownames(ketodata),dim(ketodata)[2]))
stacked_data<-data.frame(data_group,abundance,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria,levels=c('Faecalibaculum',
'Lactobacillus',
'Lactococcus',
'Akkermansia',
'Muribaculaceae__2',
'Romboutsia',
'others',
'Parasutterella',
'Enterococcus',
'A2',
'Muribaculaceae__',
'Enterorhabdus',
'Bacteroides')) 
stacked_data$data_group <- factor(stacked_data$data_group) 


colourCount = 13
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
g <- ggplot(stacked_data, aes(x=data_group,y=abundance,fill=data_bacteria))+  geom_bar(stat='identity') +   scale_fill_manual(values = getPalette(colourCount))+ylim(0,3.5e10)

ggsave(g, filename='LSI_density_keto.pdf', width=25, height=30, units=c("cm"))	


controldata<-data[data$condition=='Control',]
sampleid<-controldata$sampleid
controldata<-controldata[ , !colnames(controldata) %in% c("sampleid","condition")]
controldata<-t(controldata)
colnames(controldata)<-sampleid

controldata<-data.frame(controldata)

data_group<-c(rep('c1',dim(controldata)[1]),rep('c2',dim(controldata)[1]),rep('c3',dim(controldata)[1]),rep('c4',dim(controldata)[1]),rep('c5',dim(controldata)[1]))
abundance<-c(controldata$c1,controldata$c2,controldata$c3,controldata$c4,controldata$c5)
data_bacteria<-c(rep(rownames(controldata),dim(controldata)[2]))
stacked_data<-data.frame(data_group,abundance,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria,levels=c('Faecalibaculum',
'Lactobacillus',
'Lactococcus',
'Akkermansia',
'Muribaculaceae__2',
'Romboutsia',
'others',
'Parasutterella',
'Enterococcus',
'A2',
'Muribaculaceae__',
'Enterorhabdus',
'Bacteroides')) 
stacked_data$data_group <- factor(stacked_data$data_group) 


colourCount = 13
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
g <- ggplot(stacked_data, aes(x=data_group,y=abundance,fill=data_bacteria))+  geom_bar(stat='identity') +   scale_fill_manual(values = getPalette(colourCount))+ylim(0,3.5e10)

ggsave(g, filename='LSI_density_control.pdf', width=25, height=30, units=c("cm"))	
