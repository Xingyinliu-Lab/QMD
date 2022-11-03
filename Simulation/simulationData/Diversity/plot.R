####R code for analysis and visulize the simulation and validation result

library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)
library(plyr)
library(ggpubr)

# library(extrafont)
# font_import()
# fonts()
# loadfonts()
data<- read.table("H2029_diversity.csv", header = TRUE, sep = ",")
data$typelabel<-'Diversity Decreased'
data$typelabel[data$type<1]<-'Diversity Increased'
data$typelabel<-as.factor(data$typelabel)

p3<-ggplot(data, aes(x=typelabel, y=qmd_MAE, fill=typelabel)) + geom_boxplot()+labs(y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,1))+ stat_compare_means(method = "t.test")
ggsave(p3, device=cairo_pdf, filename='H2029_MAE_diversity.pdf', width=25, height=15, units=c("cm"))

t.test(data$qmd_MAE ~ data$typelabel)

        Welch Two Sample t-test

data:  data$qmd_MAE by data$typelabel
t = 4.0155, df = 24.99, p-value = 0.0004763
alternative hypothesis: true difference in means between group Diversity Decreased and group Diversity Increased is not equal to 0
95 percent confidence interval:
 0.07568556 0.23507564
sample estimates:
mean in group Diversity Decreased mean in group Diversity Increased 
                        0.7847755                         0.6293949 


data<- read.table("Obesity_diversity.csv", header = TRUE, sep = ",")
data$typelabel<-'Diversity Decreased'
data$typelabel[data$type<1]<-'Diversity Increased'
data$typelabel<-as.factor(data$typelabel)

p3<-ggplot(data, aes(x=typelabel, y=qmd_MAE, fill=typelabel)) + geom_boxplot()+labs(y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,1))+ stat_compare_means(method = "t.test")
ggsave(p3, device=cairo_pdf, filename='Obesity_MAE_diversity.pdf', width=25, height=15, units=c("cm"))

t.test(data$qmd_MAE ~ data$typelabel)
        Welch Two Sample t-test

data:  data$qmd_MAE by data$typelabel
t = 5.269, df = 26.528, p-value = 1.559e-05
alternative hypothesis: true difference in means between group Diversity Decreased and group Diversity Increased is not equal to 0
95 percent confidence interval:
 0.1005733 0.2290344
sample estimates:
mean in group Diversity Decreased mean in group Diversity Increased 
                        0.8524157                         0.6876118 
