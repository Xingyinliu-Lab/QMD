####R code for analysis and visulize the simulation and validation result

library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)
library(plyr)
# library(extrafont)
# font_import()
# fonts()
# loadfonts()

summary_func <- function(x, col){
         c(mean = mean(x[[col]], na.rm=TRUE), sd = sd(x[[col]], na.rm=TRUE)) 
    } 
data_summary <- function(data, varname, groupnames){ 
    data_sum<-ddply(data, groupnames, .fun=summary_func, varname)
    data_sum <- rename(data_sum, c("mean" = varname)) 
    return(data_sum) 
    }


##########################
##########################
##########################
##########################


setwd('result_Aanlysis/collectedStatistics/')
data<- read.table("H2029_stat.csv", header = TRUE, sep = ",")

# linear regression without intercept
psi_fit = lm(MicrobiotaDensityFoldChange~qmd_Density_delta-1, data=data)
summary(psi_fit)


plot(data$MicrobiotaDensityFoldChange,data$qmd_Density_delta)
p1 <- ggplot(data, aes(x=MicrobiotaDensityFoldChange, y=qmd_Density_delta)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Real Total microbial loads foldchange',y='Predicted Total microbial loads foldchange',title='Total microbial loads foldchange prediction')+ylim(-1,8)+xlim(-1,8)+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))
ggsave(p1, device=cairo_pdf, filename='H2029_microbial_loads_foldchange_prediction_performance.pdf', width=25, height=25, units=c("cm"))


mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))

mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
# p2<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + 
#     scale_fill_brewer(palette="Set2")+
#     geom_boxplot()+labs(x="Changed Taxa Rate", y ="MAE")+ theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
# ggsave(p2, device=cairo_pdf, filename='H2029_MAE_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))

df3 <- data_summary(mae_data, varname="mae_value", groupnames=c("mae_group", "mae_con")) 
df3$mae_con=as.factor(df3$mae_con) 
head(df3)

p2<- ggplot(df3, aes(x=mae_con, y=mae_value, group=mae_group, color=mae_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = mae_group),size=4)+ geom_errorbar(aes(ymin=mae_value-sd, ymax=mae_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="MAE")
ggsave(p2, device=cairo_pdf, filename='H2029_MAE_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))




describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))

describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))



mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))
mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
# p3<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + geom_boxplot()+labs(x="Group Size", y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))


df3 <- data_summary(mae_data, varname="mae_value", groupnames=c("mae_group", "mae_con")) 
df3$mae_con=as.factor(df3$mae_con) 
head(df3)

p3<- ggplot(df3, aes(x=mae_con, y=mae_value, group=mae_group, color=mae_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = mae_group),size=4)+ geom_errorbar(aes(ymin=mae_value-sd, ymax=mae_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group Size", y ="MAE")

ggsave(p3,device=cairo_pdf, filename='H2029_MAE_on_groupsize.pdf', width=25, height=15, units=c("cm"))


describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))






fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_data<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_data<-data.frame(fnr_group,fnr_data)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))

# p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group))+
#   geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FNR") +scale_fill_brewer(palette="Set2")  +
# theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)

p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group,fill=fnr_group))+
  geom_violin()+ geom_blank()+labs(x="", y ="FNR")+
  scale_color_viridis(discrete=TRUE)+scale_fill_viridis(discrete=TRUE,alpha=0.6)+
  theme_ipsum()+ theme(legend.position="none")+ylim(0,1)+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+ stat_summary(fun=median, geom="point", size=2, color="black")

ggsave(p6, device=cairo_pdf, filename='H2029_FNR.pdf', width=20, height=20, units=c("cm"))

t.test(data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)

        Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$qmd_diff_qvalue_FNR
t = -10.052, df = 768.97, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.12345876 -0.08311725
sample estimates:
mean of x mean of y 
0.0952439 0.1985319 

t.test(data$qmd_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR)

        Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$ANCOM_BC_diff_qval_FNR
t = -82.618, df = 928.21, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6875379 -0.6556318
sample estimates:
mean of x mean of y 
0.0952439 0.7668288 

t.test(data$qmd_diff_pvalue_FNR,data$ancomP_FNR)
        Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$ancomP_FNR
t = -35.245, df = 830.64, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.3470908 -0.3104706
sample estimates:
mean of x mean of y 
0.0952439 0.4240246 
t.test(data$qmd_diff_pvalue_FNR,data$mwu_diff_pvalue_FNR)

        Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$mwu_diff_pvalue_FNR
t = -9.4931, df = 994.34, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.08185351 -0.05381004
sample estimates:
mean of x mean of y 
0.0952439 0.1630757 


fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_data<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_data<-data.frame(fpr_group,fpr_data)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group))+
#   geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FPR") +scale_fill_brewer(palette="Set2") +
# theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)

p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group,fill=fpr_group))+
  geom_violin()+ geom_blank()+labs(x="", y ="FPR")+
  scale_color_viridis(discrete=TRUE)+scale_fill_viridis(discrete=TRUE,alpha=0.6)+
  theme_ipsum()+ theme(legend.position="none")+ylim(0,1)+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+ stat_summary(fun=median, geom="point", size=2, color="black")


ggsave(p7,device=cairo_pdf,  filename='H2029_FPR.pdf', width=20, height=20, units=c("cm"))

t.test(data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
        Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$qmd_diff_qvalue_FPR
t = 4.3939, df = 996.26, p-value = 1.233e-05
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.01810867 0.04733773
sample estimates:
 mean of x  mean of y 
0.08644368 0.05372048 

t.test(data$qmd_diff_pvalue_FPR,data$ancomP_FPR)

       Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$ancomP_FPR
t = -10.565, df = 722.28, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1540486 -0.1057663
sample estimates:
 mean of x  mean of y 
0.08644368 0.21635111 


t.test(data$qmd_diff_pvalue_FPR,data$mwu_diff_pvalue_FPR)

       Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$mwu_diff_pvalue_FPR
t = -38.59, df = 628.28, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6387824 -0.5769181
sample estimates:
 mean of x  mean of y 
0.08644368 0.69429389 


t.test(data$qmd_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR)

        Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$ANCOM_BC_diff_qval_FPR
t = 6.3709, df = 729.84, p-value = 3.327e-10
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.02643705 0.04998777
sample estimates:
 mean of x  mean of y 
0.08644368 0.04823127 




fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p2<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fpr_data, varname="fpr_value", groupnames=c("fpr_group", "fpr_con")) 
df3$fpr_con=as.factor(df3$fpr_con) 
p2<- ggplot(df3, aes(x=fpr_con, y=fpr_value, group=fpr_group, color=fpr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fpr_group),size=4)+ geom_errorbar(aes(ymin=fpr_value-sd, ymax=fpr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="FPR") 

ggsave(p2,device=cairo_pdf, filename='H2029_FPR_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))



describeBy(fpr_data$fpr_value,list(fpr_data$fpr_con,fpr_data$fpr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fpr_data$fpr_value,fpr_data$fpr_group,mat=TRUE,digits=2,quant=c(.25,.75))




fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize,data$groupsize)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p3<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Group size", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fpr_data, varname="fpr_value", groupnames=c("fpr_group", "fpr_con")) 
df3$fpr_con=as.factor(df3$fpr_con) 
p3<- ggplot(df3, aes(x=fpr_con, y=fpr_value, group=fpr_group, color=fpr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fpr_group),size=4)+ geom_errorbar(aes(ymin=fpr_value-sd, ymax=fpr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group size", y ="FPR")

ggsave(p3,device=cairo_pdf, filename='H2029_FPR_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(fpr_data$fpr_value,list(fpr_data$fpr_con,fpr_data$fpr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fpr_data$fpr_value,fpr_data$fpr_group,mat=TRUE,digits=2,quant=c(.25,.75))




fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_value<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fnr_data<-data.frame(fnr_group,fnr_value,fnr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fnr_data$fnr_con <- cut(fnr_data$fnr_con, breaks = group, labels =labels)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))
# p2<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fnr_data, varname="fnr_value", groupnames=c("fnr_group", "fnr_con")) 
df3$fnr_con=as.factor(df3$fnr_con) 
p2<- ggplot(df3, aes(x=fnr_con, y=fnr_value, group=fnr_group, color=fnr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fnr_group),size=4)+ geom_errorbar(aes(ymin=fnr_value-sd, ymax=fnr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="FNR")

ggsave(p2,device=cairo_pdf, filename='H2029_FNR_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))




describeBy(fnr_data$fnr_value,list(fnr_data$fnr_con,fnr_data$fnr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fnr_data$fnr_value,fnr_data$fnr_group,mat=TRUE,digits=2,quant=c(.25,.75))

fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_value<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize,data$groupsize)
fnr_data<-data.frame(fnr_group,fnr_value,fnr_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
fnr_data$fnr_con <- cut(fnr_data$fnr_con, breaks = group, labels =labels)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))
# p3<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Group size", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fnr_data, varname="fnr_value", groupnames=c("fnr_group", "fnr_con")) 
df3$fnr_con=as.factor(df3$fnr_con) 
p3<- ggplot(df3, aes(x=fnr_con, y=fnr_value, group=fnr_group, color=fnr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fnr_group),size=4)+ geom_errorbar(aes(ymin=fnr_value-sd, ymax=fnr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group size", y ="FNR")

ggsave(p3,device=cairo_pdf, filename='H2029_FNR_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(fnr_data$fnr_value,list(fnr_data$fnr_con,fnr_data$fnr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fnr_data$fnr_value,fnr_data$fnr_group,mat=TRUE,digits=2,quant=c(.25,.75))


library(scatterplot3d)
cubedraw <- function(res3d,  cex = 2, text. = FALSE)
  {
    ## Purpose: Draw nice cube with corners
    cub <- rbind(c(0,0,4),c(0,0,-4), c(100,0,-4), c(100,1,-4),c(100,1,4),c(0,1,4),  c(100,0,4), c(0,1,-4)) # <- "inner": fore- & back-ground
    ## visibile corners + lines:
    res3d$points3d(cub[c(1:6,1,7,3,7,5) ,], cex = cex, type = 'b', lty = 1)
    ## hidden corner + lines
    res3d$points3d(cub[c(2,8,4,8,6),     ], cex = cex, type = 'b', lty = 3)
    if(text.)## debug
        text(res3d$xyz.convert(cub), labels=1:nrow(cub), col='tomato', cex=2)
  }

simu_con<-data.frame(data$groupsize,data$changedTaxaRate,data$changfold)


pdf('H2029_simu_con.pdf',width = 15, height=15)
rr <- scatterplot3d(simu_con, box = FALSE, angle = 45,xlim = c(-10, 110), ylim = c(-0.1,1.1), zlim = c(-5, 5), color="steelblue",pch = 16)
cubedraw(rr)
dev.off()

p1 <- ggplot(data, aes(x=groupsize, y=changedTaxaRate)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Group size',y='changedTaxaRate',title='')+ylim(0,1)+xlim(0,100)
ggsave(p1, device=cairo_pdf, filename='H2029_groupsize_changedtaxarate.pdf', width=25, height=10, units=c("cm"))

p1 <- ggplot(data, aes(x=groupsize, y=changfold)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Group size',y='changfold',title='')+ylim(-4,4)+xlim(0,100)
ggsave(p1, device=cairo_pdf, filename='H2029_groupsize_changfold.pdf', width=25, height=10, units=c("cm"))

p1 <- ggplot(data, aes(x=changedTaxaRate, y=changfold)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='changedTaxaRate',y='changfold',title='')+ylim(-4,4)+xlim(0,1)
ggsave(p1, device=cairo_pdf, filename='H2029_changfold_changedtaxarate.pdf', width=25, height=10, units=c("cm"))


##########################
##########################
##########################
##########################


data<- read.table("Obesity_stat.csv", header = TRUE, sep = ",")

# linear regression without intercept
psi_fit = lm(MicrobiotaDensityFoldChange~qmd_Density_delta-1, data=data)
summary(psi_fit)

plot(data$MicrobiotaDensityFoldChange,data$qmd_Density_delta)
p1 <- ggplot(data, aes(x=MicrobiotaDensityFoldChange, y=qmd_Density_delta)) + 
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Real Total microbial loads foldchange',y='Predicted Total microbial loads foldchange',title='Total microbial loads foldchange prediction')+ylim(-1,8)+xlim(-1,8)+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))

ggsave(p1, device=cairo_pdf, filename='Obesity_microbial_loads_foldchange_prediction_performance.pdf', width=25, height=25, units=c("cm"))


mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))

mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
# p2<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + 
#     scale_fill_brewer(palette="Set2")+
#     geom_boxplot()+labs(x="Changed Taxa Rate", y ="MAE")+ theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
df3 <- data_summary(mae_data, varname="mae_value", groupnames=c("mae_group", "mae_con")) 
df3$mae_con=as.factor(df3$mae_con) 
head(df3)
p2<- ggplot(df3, aes(x=mae_con, y=mae_value, group=mae_group, color=mae_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = mae_group),size=4)+ geom_errorbar(aes(ymin=mae_value-sd, ymax=mae_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="MAE")

ggsave(p2, device=cairo_pdf, filename='Obesity_MAE_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))


describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))

mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))
mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
# p3<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + geom_boxplot()+labs(x="Group Size", y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
df3 <- data_summary(mae_data, varname="mae_value", groupnames=c("mae_group", "mae_con")) 
df3$mae_con=as.factor(df3$mae_con) 
head(df3)
p3<- ggplot(df3, aes(x=mae_con, y=mae_value, group=mae_group, color=mae_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = mae_group),size=4)+ geom_errorbar(aes(ymin=mae_value-sd, ymax=mae_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group Size", y ="MAE") 

ggsave(p3,device=cairo_pdf, filename='Obesity_MAE_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))




fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_data<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_data<-data.frame(fnr_group,fnr_data)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))

# p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group))+
#   geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FNR") +scale_fill_brewer(palette="Set2")  +
# theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)

p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group,fill=fnr_group))+
  geom_violin()+ geom_blank()+labs(x="", y ="FNR")+
  scale_color_viridis(discrete=TRUE)+scale_fill_viridis(discrete=TRUE,alpha=0.6)+
  theme_ipsum()+ theme(legend.position="none")+ylim(0,1)+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+ stat_summary(fun=median, geom="point", size=2, color="black")

ggsave(p6, device=cairo_pdf, filename='Obesity_FNR.pdf', width=20, height=20, units=c("cm"))

t.test(data$qmd_diff_pvalue_FNR,data$ancomP_FNR)

        Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$ancomP_FNR
t = -39.092, df = 909.08, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.3792266 -0.3429697
sample estimates:
mean of x mean of y 
0.1069659 0.4680641 

t.test(data$qmd_diff_pvalue_FNR,data$mwu_diff_pvalue_FNR)
   Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$mwu_diff_pvalue_FNR
t = -9.1159, df = 993.18, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.08407467 -0.05428946
sample estimates:
mean of x mean of y 
0.1069659 0.1761480 

t.test(data$qmd_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR)
    Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$ANCOM_BC_diff_qval_FNR
t = -100.22, df = 987.55, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.7565196 -0.7274620
sample estimates:
mean of x mean of y 
0.1069659 0.8489568 


t.test(data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
     Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$qmd_diff_qvalue_FNR
t = -9.7925, df = 791.04, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1288074 -0.0857900
sample estimates:
mean of x mean of y 
0.1069659 0.2142646 


fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_data<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_data<-data.frame(fpr_group,fpr_data)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group))+
#   geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FPR") +scale_fill_brewer(palette="Set2") +
# theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)
p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group,fill=fpr_group))+
  geom_violin()+ geom_blank()+labs(x="", y ="FPR")+
  scale_color_viridis(discrete=TRUE)+scale_fill_viridis(discrete=TRUE,alpha=0.6)+
  theme_ipsum()+ theme(legend.position="none")+ylim(0,1)+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+ stat_summary(fun=median, geom="point", size=2, color="black")


ggsave(p7,device=cairo_pdf,  filename='Obesity_FPR.pdf', width=20, height=20, units=c("cm"))


t.test(data$qmd_diff_pvalue_FPR,data$ancomP_FPR)
     Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$ancomP_FPR
t = -7.408, df = 796.26, p-value = 3.274e-13
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.09432386 -0.05480750
sample estimates:
 mean of x  mean of y 
0.08362022 0.15818590 

t.test(data$qmd_diff_pvalue_FPR,data$mwu_diff_pvalue_FPR)
   Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$mwu_diff_pvalue_FPR
t = -40.918, df = 622.14, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6418272 -0.5830419
sample estimates:
 mean of x  mean of y 
0.08362022 0.69605481 

t.test(data$qmd_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR)
     Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$ANCOM_BC_diff_qval_FPR
t = 7.9789, df = 711.2, p-value = 5.914e-15
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.03358565 0.05550841
sample estimates:
 mean of x  mean of y 
0.08362022 0.03907319 

t.test(data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
     Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$qmd_diff_qvalue_FPR
t = 5.0991, df = 988.99, p-value = 4.089e-07
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.02160689 0.04864173
sample estimates:
 mean of x  mean of y 
0.08362022 0.04849591 



fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p2<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fpr_data, varname="fpr_value", groupnames=c("fpr_group", "fpr_con")) 
df3$fpr_con=as.factor(df3$fpr_con) 
p2<- ggplot(df3, aes(x=fpr_con, y=fpr_value, group=fpr_group, color=fpr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fpr_group),size=4)+ geom_errorbar(aes(ymin=fpr_value-sd, ymax=fpr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="FPR") 

ggsave(p2,device=cairo_pdf, filename='Obesity_FPR_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))


describeBy(fpr_data$fpr_value,list(fpr_data$fpr_con,fpr_data$fpr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fpr_data$fpr_value,fpr_data$fpr_group,mat=TRUE,digits=2,quant=c(.25,.75))


fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize,data$groupsize)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p3<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Group size", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fpr_data, varname="fpr_value", groupnames=c("fpr_group", "fpr_con")) 
df3$fpr_con=as.factor(df3$fpr_con) 
p3<- ggplot(df3, aes(x=fpr_con, y=fpr_value, group=fpr_group, color=fpr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fpr_group),size=4)+ geom_errorbar(aes(ymin=fpr_value-sd, ymax=fpr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group size", y ="FPR")

ggsave(p3,device=cairo_pdf, filename='Obesity_FPR_on_groupsize.pdf', width=25, height=15, units=c("cm"))


describeBy(fpr_data$fpr_value,list(fpr_data$fpr_con,fpr_data$fpr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fpr_data$fpr_value,fpr_data$fpr_group,mat=TRUE,digits=2,quant=c(.25,.75))



fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_value<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fnr_data<-data.frame(fnr_group,fnr_value,fnr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fnr_data$fnr_con <- cut(fnr_data$fnr_con, breaks = group, labels =labels)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))
# p2<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fnr_data, varname="fnr_value", groupnames=c("fnr_group", "fnr_con")) 
df3$fnr_con=as.factor(df3$fnr_con) 
p2<- ggplot(df3, aes(x=fnr_con, y=fnr_value, group=fnr_group, color=fnr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fnr_group),size=4)+ geom_errorbar(aes(ymin=fnr_value-sd, ymax=fnr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="FNR")

ggsave(p2,device=cairo_pdf, filename='Obesity_FNR_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))

describeBy(fnr_data$fnr_value,list(fnr_data$fnr_con,fnr_data$fnr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fnr_data$fnr_value,fnr_data$fnr_group,mat=TRUE,digits=2,quant=c(.25,.75))


fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_value<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize,data$groupsize)
fnr_data<-data.frame(fnr_group,fnr_value,fnr_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
fnr_data$fnr_con <- cut(fnr_data$fnr_con, breaks = group, labels =labels)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))
# p3<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Group size", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fnr_data, varname="fnr_value", groupnames=c("fnr_group", "fnr_con")) 
df3$fnr_con=as.factor(df3$fnr_con) 
p3<- ggplot(df3, aes(x=fnr_con, y=fnr_value, group=fnr_group, color=fnr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fnr_group),size=4)+ geom_errorbar(aes(ymin=fnr_value-sd, ymax=fnr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group size", y ="FNR")

ggsave(p3,device=cairo_pdf, filename='Obesity_FNR_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(fnr_data$fnr_value,list(fnr_data$fnr_con,fnr_data$fnr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fnr_data$fnr_value,fnr_data$fnr_group,mat=TRUE,digits=2,quant=c(.25,.75))




library(scatterplot3d)
cubedraw <- function(res3d,  cex = 2, text. = FALSE)
  {
    ## Purpose: Draw nice cube with corners
    cub <- rbind(c(0,0,4),c(0,0,-4), c(100,0,-4), c(100,1,-4),c(100,1,4),c(0,1,4),  c(100,0,4), c(0,1,-4)) # <- "inner": fore- & back-ground
    ## visibile corners + lines:
    res3d$points3d(cub[c(1:6,1,7,3,7,5) ,], cex = cex, type = 'b', lty = 1)
    ## hidden corner + lines
    res3d$points3d(cub[c(2,8,4,8,6),     ], cex = cex, type = 'b', lty = 3)
    if(text.)## debug
        text(res3d$xyz.convert(cub), labels=1:nrow(cub), col='tomato', cex=2)
  }

simu_con<-data.frame(data$groupsize,data$changedTaxaRate,data$changfold)


pdf('Obesity_simu_con.pdf',width = 15, height=15)
rr <- scatterplot3d(simu_con, box = FALSE, angle = 45,xlim = c(-10, 110), ylim = c(-0.1,1.1), zlim = c(-5, 5), color="steelblue",pch = 16)
cubedraw(rr)
dev.off()


p1 <- ggplot(data, aes(x=groupsize, y=changedTaxaRate)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Group size',y='changedTaxaRate',title='')+ylim(0,1)+xlim(0,100)
ggsave(p1, device=cairo_pdf, filename='Obesity_groupsize_changedtaxarate.pdf', width=25, height=10, units=c("cm"))

p1 <- ggplot(data, aes(x=groupsize, y=changfold)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Group size',y='changfold',title='')+ylim(-4,4)+xlim(0,100)
ggsave(p1, device=cairo_pdf, filename='Obesity_groupsize_changfold.pdf', width=25, height=10, units=c("cm"))

p1 <- ggplot(data, aes(x=changedTaxaRate, y=changfold)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='changedTaxaRate',y='changfold',title='')+ylim(-4,4)+xlim(0,1)
ggsave(p1, device=cairo_pdf, filename='Obesity_changfold_changedtaxarate.pdf', width=25, height=10, units=c("cm"))

##########################
##########################
##########################
##########################
data<- read.table("GP_stat.csv", header = TRUE, sep = ",")
# linear regression without intercept
psi_fit = lm(MicrobiotaDensityFoldChange~qmd_Density_delta-1, data=data)
summary(psi_fit)

plot(data$MicrobiotaDensityFoldChange,data$qmd_Density_delta)
p1 <- ggplot(data, aes(x=MicrobiotaDensityFoldChange, y=qmd_Density_delta)) + 
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Real Total microbial loads foldchange',y='Predicted Total microbial loads foldchange',title='Total microbial loads foldchange prediction')+ylim(-1,8)+xlim(-1,8)+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))


ggsave(p1, device=cairo_pdf, filename='GP_microbial_loads_foldchange_prediction_performance.pdf', width=25, height=25, units=c("cm"))




mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))

mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
# p2<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + 
#     scale_fill_brewer(palette="Set2")+
#     geom_boxplot()+labs(x="Changed Taxa Rate", y ="MAE")+ theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
df3 <- data_summary(mae_data, varname="mae_value", groupnames=c("mae_group", "mae_con")) 
df3$mae_con=as.factor(df3$mae_con) 
head(df3)
p2<- ggplot(df3, aes(x=mae_con, y=mae_value, group=mae_group, color=mae_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = mae_group),size=4)+ geom_errorbar(aes(ymin=mae_value-sd, ymax=mae_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="MAE")

ggsave(p2, device=cairo_pdf, filename='GP_MAE_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))



describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))

mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))
mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
# p3<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + geom_boxplot()+labs(x="Group Size", y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
df3 <- data_summary(mae_data, varname="mae_value", groupnames=c("mae_group", "mae_con")) 
df3$mae_con=as.factor(df3$mae_con) 
head(df3)
p3<- ggplot(df3, aes(x=mae_con, y=mae_value, group=mae_group, color=mae_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = mae_group),size=4)+ geom_errorbar(aes(ymin=mae_value-sd, ymax=mae_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group Size", y ="MAE") 


ggsave(p3,device=cairo_pdf, filename='GP_MAE_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))



fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_data<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_data<-data.frame(fnr_group,fnr_data)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))

# p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group))+
#   geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FNR") +scale_fill_brewer(palette="Set2")  +
# theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)
p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group,fill=fnr_group))+
  geom_violin(scale='width')+ geom_blank()+labs(x="", y ="FNR")+
  scale_color_viridis(discrete=TRUE)+scale_fill_viridis(discrete=TRUE,alpha=0.6)+
  theme_ipsum()+ theme(legend.position="none")+ylim(0,1)+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+ stat_summary(fun=median, geom="point", size=2, color="black")

ggsave(p6, device=cairo_pdf, filename='GP_FNR.pdf', width=20, height=20, units=c("cm"))

t.test(data$qmd_diff_pvalue_FNR,data$ancomP_FNR)
    Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$ancomP_FNR
t = -0.38881, df = 993.38, p-value = 0.6975
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.01969680  0.01318227
sample estimates:
 mean of x  mean of y 
0.03410062 0.03735789 


t.test(data$qmd_diff_pvalue_FNR,data$mwu_diff_pvalue_FNR)
      Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$mwu_diff_pvalue_FNR
t = 5.7565, df = 500.42, p-value = 1.499e-08
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.02169792 0.04418342
sample estimates:
  mean of x   mean of y 
0.034100624 0.001159953 

t.test(data$qmd_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR)

     Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$ANCOM_BC_diff_qval_FNR
t = -7.8397, df = 990.42, p-value = 1.163e-14
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.08298109 -0.04975558
sample estimates:
 mean of x  mean of y 
0.03410062 0.10046896 

t.test(data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
     Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FNR and data$qmd_diff_qvalue_FNR
t = -0.25574, df = 994.57, p-value = 0.7982
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.0184886  0.0142253
sample estimates:
 mean of x  mean of y 
0.03410062 0.03623227 



fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_data<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_data<-data.frame(fpr_group,fpr_data)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group))+
#   geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FPR") +scale_fill_brewer(palette="Set2") +
# theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)
p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group,fill=fpr_group))+
  geom_violin(scale='width')+ 
  geom_blank()+labs(x="", y ="FPR")+
  scale_color_viridis(discrete=TRUE)+scale_fill_viridis(discrete=TRUE,alpha=0.6)+
  theme_ipsum()+ theme(legend.position="none")+ylim(0,1)+theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+ stat_summary(fun=median, geom="point", size=2, color="black")

ggsave(p7,device=cairo_pdf,  filename='GP_FPR.pdf', width=20, height=20, units=c("cm"))


t.test(data$qmd_diff_pvalue_FPR,data$ancomP_FPR)
     Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$ancomP_FPR
t = -3.2132, df = 797.42, p-value = 0.001365
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.10326098 -0.02494205
sample estimates:
mean of x mean of y 
0.2239684 0.2880700 

t.test(data$qmd_diff_pvalue_FPR,data$mwu_diff_pvalue_FPR)
   Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$mwu_diff_pvalue_FPR
t = -77.918, df = 499, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.7955779 -0.7564431
sample estimates:
mean of x mean of y 
0.2239684 0.9999789 

t.test(data$qmd_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR)
      Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$ANCOM_BC_diff_qval_FPR
t = 7.303, df = 969.23, p-value = 5.866e-13
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.08267825 0.14343904
sample estimates:
mean of x mean of y 
0.2239684 0.1109098 

t.test(data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)

      Welch Two Sample t-test

data:  data$qmd_diff_pvalue_FPR and data$qmd_diff_qvalue_FPR
t = 2.5749, df = 996.55, p-value = 0.01017
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.00879747 0.06516008
sample estimates:
mean of x mean of y 
0.2239684 0.1869897 




fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p2<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fpr_data, varname="fpr_value", groupnames=c("fpr_group", "fpr_con")) 
df3$fpr_con=as.factor(df3$fpr_con) 
p2<- ggplot(df3, aes(x=fpr_con, y=fpr_value, group=fpr_group, color=fpr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fpr_group),size=4)+ geom_errorbar(aes(ymin=fpr_value-sd, ymax=fpr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="FPR") 

ggsave(p2,device=cairo_pdf, filename='GP_FPR_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))

describeBy(fpr_data$fpr_value,list(fpr_data$fpr_con,fpr_data$fpr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fpr_data$fpr_value,fpr_data$fpr_group,mat=TRUE,digits=2,quant=c(.25,.75))


fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize,data$groupsize)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
# p3<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Group size", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fpr_data, varname="fpr_value", groupnames=c("fpr_group", "fpr_con")) 
df3$fpr_con=as.factor(df3$fpr_con) 
p3<- ggplot(df3, aes(x=fpr_con, y=fpr_value, group=fpr_group, color=fpr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fpr_group),size=4)+ geom_errorbar(aes(ymin=fpr_value-sd, ymax=fpr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group size", y ="FPR")

ggsave(p3,device=cairo_pdf, filename='GP_FPR_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(fpr_data$fpr_value,list(fpr_data$fpr_con,fpr_data$fpr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fpr_data$fpr_value,fpr_data$fpr_group,mat=TRUE,digits=2,quant=c(.25,.75))


fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_value<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fnr_data<-data.frame(fnr_group,fnr_value,fnr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fnr_data$fnr_con <- cut(fnr_data$fnr_con, breaks = group, labels =labels)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))
# p2<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fnr_data, varname="fnr_value", groupnames=c("fnr_group", "fnr_con")) 
df3$fnr_con=as.factor(df3$fnr_con) 
p2<- ggplot(df3, aes(x=fnr_con, y=fnr_value, group=fnr_group, color=fnr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fnr_group),size=4)+ geom_errorbar(aes(ymin=fnr_value-sd, ymax=fnr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Changed Taxa Rate", y ="FNR")

ggsave(p2,device=cairo_pdf, filename='GP_FNR_on_ChangedTaxaRate.pdf', width=25, height=15, units=c("cm"))
describeBy(fnr_data$fnr_value,list(fnr_data$fnr_con,fnr_data$fnr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fnr_data$fnr_value,fnr_data$fnr_group,mat=TRUE,digits=2,quant=c(.25,.75))


fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_value<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_con<-c(data$groupsize,data$groupsize,data$groupsize,data$groupsize,data$groupsize)
fnr_data<-data.frame(fnr_group,fnr_value,fnr_con)
group <- seq(0,100,by=10)
labels<-c('0~10','10~20','20~30','30~40','40~50','50~60','60~70','70~80','80~90','90~100')
fnr_data$fnr_con <- cut(fnr_data$fnr_con, breaks = group, labels =labels)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))
# p3<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Group size", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
df3 <- data_summary(fnr_data, varname="fnr_value", groupnames=c("fnr_group", "fnr_con")) 
df3$fnr_con=as.factor(df3$fnr_con) 
p3<- ggplot(df3, aes(x=fnr_con, y=fnr_value, group=fnr_group, color=fnr_group)) + geom_line(position=position_dodge(0.3)) + geom_point(position=position_dodge(0.3),aes(shape = fnr_group),size=4)+ geom_errorbar(aes(ymin=fnr_value-sd, ymax=fnr_value+sd), width=.6, position=position_dodge(0.3))+theme_ipsum()+scale_fill_brewer(palette="Set2") +theme(axis.line.x=element_line(linetype=1,color="black",size=1))+theme(axis.line.y=element_line(linetype=1,color="black",size=1))+theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+labs(x="Group size", y ="FNR")

ggsave(p3,device=cairo_pdf, filename='GP_FNR_on_groupsize.pdf', width=25, height=15, units=c("cm"))
describeBy(fnr_data$fnr_value,list(fnr_data$fnr_con,fnr_data$fnr_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(fnr_data$fnr_value,fnr_data$fnr_group,mat=TRUE,digits=2,quant=c(.25,.75))





library(scatterplot3d)
cubedraw <- function(res3d,  cex = 2, text. = FALSE)
  {
    ## Purpose: Draw nice cube with corners
    cub <- rbind(c(0,0,4),c(0,0,-4), c(100,0,-4), c(100,1,-4),c(100,1,4),c(0,1,4),  c(100,0,4), c(0,1,-4)) # <- "inner": fore- & back-ground
    ## visibile corners + lines:
    res3d$points3d(cub[c(1:6,1,7,3,7,5) ,], cex = cex, type = 'b', lty = 1)
    ## hidden corner + lines
    res3d$points3d(cub[c(2,8,4,8,6),     ], cex = cex, type = 'b', lty = 3)
    if(text.)## debug
        text(res3d$xyz.convert(cub), labels=1:nrow(cub), col='tomato', cex=2)
  }

simu_con<-data.frame(data$groupsize,data$changedTaxaRate,data$changfold)


pdf('GP_simu_con.pdf',width = 15, height=15)
rr <- scatterplot3d(simu_con, box = FALSE, angle = 45,xlim = c(-10, 110), ylim = c(-0.1,1.1), zlim = c(-5, 5), color="steelblue",pch = 16)
cubedraw(rr)
dev.off()

p1 <- ggplot(data, aes(x=groupsize, y=changedTaxaRate)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Group size',y='changedTaxaRate',title='')+ylim(0,1)+xlim(0,100)
ggsave(p1, device=cairo_pdf, filename='GP_groupsize_changedtaxarate.pdf', width=25, height=10, units=c("cm"))

p1 <- ggplot(data, aes(x=groupsize, y=changfold)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Group size',y='changfold',title='')+ylim(-4,4)+xlim(0,100)
ggsave(p1, device=cairo_pdf, filename='GP_groupsize_changfold.pdf', width=25, height=10, units=c("cm"))

p1 <- ggplot(data, aes(x=changedTaxaRate, y=changfold)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='changedTaxaRate',y='changfold',title='')+ylim(-4,4)+xlim(0,1)
ggsave(p1, device=cairo_pdf, filename='GP_changfold_changedtaxarate.pdf', width=25, height=10, units=c("cm"))
