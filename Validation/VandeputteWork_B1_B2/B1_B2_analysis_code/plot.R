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

################
library(ggplot2)
theme_set(theme_bw())  
library(pheatmap)

data <- read.table("dc_b1_compare.csv", header=T, sep = "," )
data <- data[order(data$Microbial.density.change), ]  # sort
data_group<-c(rep('Microbial density change',dim(data)[1]),rep('QMD',dim(data)[1]),rep('DR',dim(data)[1]),rep('ANCOM_BC',dim(data)[1]))
densitychange<-c(data$Microbial.density.change,data$QMD,data$DR,data$ANCOM.BC)
data_bacteria<-c(data$shortname,data$taxaName,data$shortname,data$taxaName)
stacked_data<-data.frame(data_group,densitychange,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria , levels = data$shortname ) 
stacked_data$data_group <- factor(stacked_data$data_group , levels =c('Microbial density change','QMD','DR','ANCOM_BC') ) 

p1<-ggplot(stacked_data, aes(x=data_bacteria, y=densitychange, colour=data_group,group=data_group,shape=data_group)) +
 geom_line(size=1)+geom_point(size=2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(p1, filename='dc_b1_compare.pdf', width=25, height=10, units=c("cm"))	


data <- read.table("dc_b2_compare.csv", header=T, sep = "," )
data <- data[order(data$Microbial.density.change), ]  # sort
data_group<-c(rep('Microbial density change',dim(data)[1]),rep('QMD',dim(data)[1]),rep('DR',dim(data)[1]),rep('ANCOM_BC',dim(data)[1]))
densitychange<-c(data$Microbial.density.change,data$QMD,data$DR,data$ANCOM.BC)
data_bacteria<-c(data$shortname,data$taxaName,data$shortname,data$taxaName)
stacked_data<-data.frame(data_group,densitychange,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria , levels = data$shortname ) 
stacked_data$data_group <- factor(stacked_data$data_group , levels =c('Microbial density change','QMD','DR','ANCOM_BC') ) 

p2<-ggplot(stacked_data, aes(x=data_bacteria, y=densitychange, colour=data_group,group=data_group,shape=data_group)) +
 geom_line(size=1)+geom_point(size=2)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


ggsave(p2, filename='dc_b2_compare.pdf', width=25, height=10, units=c("cm"))	