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

,colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric","character","character","character","character","character","character"))

data<- read.table("result_stacked.csv", header = TRUE, sep = "," ,colClasses=c("character","numeric","numeric","numeric","numeric","numeric"))
data <- data[order(data$qmd_logged_loads_diff), ]  # sort
data_group<-c(rep('qmd_logged_loads_diff',dim(data)[1]),rep('logged_relative_abundance_diff',dim(data)[1]))
foldchange<-c(data$qmd_logged_loads_diff,data$logged_relative_abundance_diff)
data_bacteria<-c(data$item,data$item)
stacked_data<-data.frame(data_group,foldchange,data_bacteria)
stacked_data$data_bacteria <- factor(stacked_data$data_bacteria , levels = data$item ) 
stacked_data$data_group <- factor(stacked_data$data_group , levels =c('qmd_logged_loads_diff','logged_relative_abundance_diff') ) 

p1<-ggplot(stacked_data, aes(x=data_bacteria, y=foldchange, colour=data_group,group=data_group,shape=data_group)) + geom_line(size=1)+coord_flip()+geom_point(size=1)

ggsave(p1, filename='IBS_foldchange.pdf', width=20, height=50, units=c("cm"))	

data <- data[order(-data$qmd_logged_loads_diff), ] 
data$mwu_logged_pvalue_reject<-0
data$mwu_logged_pvalue_reject[data$mwu_logged_pvalue<0.05]<-1

data$qmd_diff_pvalue_reject<-0
data$qmd_diff_pvalue_reject[data$qmd_diff_pvalue<0.05]<-1

diffTaxa<-data.frame(data$ancomP,data$mwu_logged_pvalue_reject,data$qmd_diff_pvalue_reject)
rownames(diffTaxa)<-data$item
save_pheatmap <- function(x, filename, width=10, height=50) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename,width = width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
p2<-pheatmap(diffTaxa,cluster_rows=FALSE,cluster_cols=FALSE)
save_pheatmap(p2,'crc_diff.pdf')

