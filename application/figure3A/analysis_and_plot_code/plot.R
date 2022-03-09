------------------
library(randomForest)
library(dplyr)



a<- read.table("total_microbial_change.csv", header = TRUE, sep = ",")
a$d <- as.factor(a$d)
head(a)
library(ggplot2)
# Basic violin plot
p <- ggplot(a, aes(x=d, y=total_microbial_abundance_change)) + 
  geom_violin()
p







limitslist<-c("CRC","Obesity", 'UC','IBD',"Crohn Disease")
limitslabels<-c("CRC","Obesity", 'UC','IBD',"Crohn Disease")

p<-ggplot(a, aes(x=d, y=total_microbial_abundance_change,color=d))+
  geom_violin(trim=FALSE,scale='width')+ labs(title="",x="", y = "Total microbial abundance change") + 
   geom_blank()+ scale_x_discrete(limits=limitslist,labels=limitslabels) + theme_minimal()+
  stat_summary(fun=median, geom="point", size=2, color="red")+theme(legend.position="none") 


ggsave(p, filename='total_microbial_abundance_change.pdf', width=20, height=10, units=c("cm"))	

