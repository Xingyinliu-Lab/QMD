library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)


##############################
##############################
##############################
##############################
##############################validation
##############################
##############################
##############################
#Diverging bars


library(ggplot2)
theme_set(theme_bw())  

data<- read.table("LSI_res_forPlot.csv", header = TRUE, sep = "," ,colClasses=c("character","character","numeric","numeric","numeric","numeric"))
data <- data[order(data$observed_microbiota_loads_foldchange), ]  # sort

data_group<-c(rep('observed_microbiota_loads_foldchange',dim(data)[1]),rep('relative_abundance_foldchange',dim(data)[1]),rep('corrected_by_observed_total_microbiota_loads_foldchange',dim(data)[1]),rep('corrected_by_predicted_total_microbiota_loads_foldchange',dim(data)[1]))
foldchange<-c(data$observed_microbiota_loads_foldchange,data$relative_abundance_foldchange,data$corrected_by_observed_total_microbiota_loads_foldchange,data$corrected_by_predicted_total_microbiota_loads_foldchange)
data_bacteria<-c(data$bacteria,data$bacteria,data$bacteria,data$bacteria)
stacked_data<-data.frame(data_group,foldchange,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria , levels = data$bacteria ) 
stacked_data$data_group <- factor(stacked_data$data_group , levels =c('observed_microbiota_loads_foldchange','corrected_by_observed_total_microbiota_loads_foldchange','corrected_by_predicted_total_microbiota_loads_foldchange','relative_abundance_foldchange') ) 


p1<-ggplot(stacked_data, aes(fill=data_group, y=foldchange, x=data_bacteria)) +  geom_bar(position="dodge", stat="identity", width=.5)+coord_flip()+scale_fill_brewer(palette="Set2")

ggsave(p1, filename='LSI_foldchange_validation.pdf', width=25, height=15, units=c("cm"))	


data<- read.table("STOOL_res_forPlot.csv", header = TRUE, sep = "," ,colClasses=c("character","character","numeric","numeric","numeric","numeric"))
# data$abs_change<-abs(data$observed_microbiota_loads_foldchange)
# data <- data[order(data$abs_change), ]  # sort
# data<-data[c(1:10),]
data <- data[order(data$observed_microbiota_loads_foldchange), ]  # sort

data_group<-c(rep('observed_microbiota_loads_foldchange',dim(data)[1]),rep('relative_abundance_foldchange',dim(data)[1]),rep('corrected_by_observed_total_microbiota_loads_foldchange',dim(data)[1]),rep('corrected_by_predicted_total_microbiota_loads_foldchange',dim(data)[1]))
foldchange<-c(data$observed_microbiota_loads_foldchange,data$relative_abundance_foldchange,data$corrected_by_observed_total_microbiota_loads_foldchange,data$corrected_by_predicted_total_microbiota_loads_foldchange)
data_bacteria<-c(data$bacteria,data$bacteria,data$bacteria,data$bacteria)
stacked_data<-data.frame(data_group,foldchange,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria , levels = data$bacteria ) 
stacked_data$data_group <- factor(stacked_data$data_group , levels =c('observed_microbiota_loads_foldchange','corrected_by_observed_total_microbiota_loads_foldchange','corrected_by_predicted_total_microbiota_loads_foldchange','relative_abundance_foldchange') ) 


p2<-ggplot(stacked_data, aes(fill=data_group, y=foldchange, x=data_bacteria)) +  geom_bar(position="dodge", stat="identity", width=.7)+coord_flip()+scale_fill_brewer(palette="Set2")

ggsave(p2, filename='STOOL_foldchange_validation.pdf', width=25, height=15, units=c("cm"))	



library(ggplot2)
theme_set(theme_bw()) 
data<- read.table("LSI_cost.csv", header = TRUE, sep = "," )
LSI_observed_microbiota_loads_foldchange<--1.28
p1<-ggplot(data, aes(x=Phi, y=Target)) + geom_line(col="#00ba38")+  theme_minimal()+ geom_vline(xintercept = c(LSI_observed_microbiota_loads_foldchange,data$Phi[which(data$Target== min(data$Target), arr.ind = TRUE)]))+xlim(-2,-1)+ylim(1.725,1.775)
p2<-ggplot(data, aes(x=Phi, y=Target)) + geom_line(col="#00ba38")+  theme_minimal()+ geom_vline(xintercept = c(LSI_observed_microbiota_loads_foldchange,data$Phi[which(data$Target== min(data$Target), arr.ind = TRUE)]))
ggsave(p1, filename='LSI_cost_magnifier.pdf', width=25, height=15, units=c("cm"))	
ggsave(p2, filename='LSI_cost.pdf', width=25, height=15, units=c("cm"))	


data<- read.table("STOOL_cost.csv", header = TRUE, sep = "," )
STOOL_observed_microbiota_loads_foldchange<--1.514
p1<-ggplot(data, aes(x=Phi, y=Target)) + geom_line(col="#00ba38")+  theme_minimal()+ geom_vline(xintercept = c(STOOL_observed_microbiota_loads_foldchange,data$Phi[which(data$Target== min(data$Target), arr.ind = TRUE)]))+xlim(-2,-0.25)+ylim(1.675,1.8)
p2<-ggplot(data, aes(x=Phi, y=Target)) + geom_line(col="#00ba38")+  theme_minimal()+ geom_vline(xintercept = c(STOOL_observed_microbiota_loads_foldchange,data$Phi[which(data$Target== min(data$Target), arr.ind = TRUE)]))
ggsave(p1, filename='STOOL_cost_magnifier.pdf', width=25, height=15, units=c("cm"))	
ggsave(p2, filename='STOOL_cost.pdf', width=25, height=15, units=c("cm"))	




################
library(ggplot2)
theme_set(theme_bw())  
library(pheatmap)

data <- read.table("LSI_compare.csv", header=T, sep = "," )
data <- data[order(data$Microbial.density.change), ]  # sort
data_group<-c(rep('Microbial density change',dim(data)[1]),rep('QMD',dim(data)[1]),rep('DR',dim(data)[1]),rep('ANCOM_BC',dim(data)[1]))
densitychange<-c(data$Microbial.density.change,data$QMD,data$DR,data$ANCOM.BC)
data_bacteria<-c(data$shortname,data$taxaName,data$shortname,data$taxaName)
stacked_data<-data.frame(data_group,densitychange,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria , levels = data$shortname ) 
stacked_data$data_group <- factor(stacked_data$data_group , levels =c('Microbial density change','QMD','DR','ANCOM_BC') ) 

p1<-ggplot(stacked_data, aes(x=data_bacteria, y=densitychange, colour=data_group,group=data_group,shape=data_group)) +
 geom_line(size=1)+geom_point(size=2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(p1, filename='LSI_compare.pdf', width=25, height=10, units=c("cm"))	


data <- read.table("STOOL_compare.csv", header=T, sep = "," )
data <- data[order(data$Microbial.density.change), ]  # sort
data_group<-c(rep('Microbial density change',dim(data)[1]),rep('QMD',dim(data)[1]),rep('DR',dim(data)[1]),rep('ANCOM_BC',dim(data)[1]))
densitychange<-c(data$Microbial.density.change,data$QMD,data$DR,data$ANCOM.BC)
data_bacteria<-c(data$shortname,data$taxaName,data$shortname,data$taxaName)
stacked_data<-data.frame(data_group,densitychange,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria , levels = data$shortname ) 
stacked_data$data_group <- factor(stacked_data$data_group , levels =c('Microbial density change','QMD','DR','ANCOM_BC') ) 

p2<-ggplot(stacked_data, aes(x=data_bacteria, y=densitychange, colour=data_group,group=data_group,shape=data_group)) +
 geom_line(size=1)+geom_point(size=2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(p2, filename='STOOL_compare.pdf', width=25, height=10, units=c("cm"))	