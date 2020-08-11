####R code for analysis and visulize the simulation and validation result

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


setwd('result_Aanlysis/collectedStatistics/')
data<- read.table("H2029_stat.csv", header = TRUE, sep = ",")

# linear regression without intercept
psi_fit = lm(MicrobiotaDensityFoldChange~qmd_Density_delta-1, data=data)
summary(psi_fit)


plot(data$MicrobiotaDensityFoldChange,data$qmd_Density_delta)
p1 <- ggplot(data, aes(x=MicrobiotaDensityFoldChange, y=qmd_Density_delta)) +
  geom_point( color="#69b3a2",size=0.6) +
  theme_ipsum()+labs(x='Real Total microbial loads foldchange',y='Predicted Total microbial loads foldchange',title='Total microbial loads foldchange prediction')+ylim(-1,8)+xlim(-1,8)
ggsave(p1, device=cairo_pdf, filename='H2029_microbial_loads_foldchange_prediction_performance.pdf', width=25, height=25, units=c("cm"))


mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))

mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
p2<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + 
    scale_fill_brewer(palette="Set2")+
    geom_boxplot()+labs(x="Changed Taxa Rate", y ="MAE")+ theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
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
p3<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + geom_boxplot()+labs(x="Group Size", y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
ggsave(p3,device=cairo_pdf, filename='H2029_MAE_on_groupsize.pdf', width=25, height=15, units=c("cm"))


describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))






fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_data<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_data<-data.frame(fnr_group,fnr_data)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))

p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group))+
  geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FNR") +scale_fill_brewer(palette="Set2")  +
theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)
ggsave(p6, device=cairo_pdf, filename='H2029_FNR.pdf', width=20, height=20, units=c("cm"))


fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_data<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_data<-data.frame(fpr_group,fpr_data)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group))+
  geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FPR") +scale_fill_brewer(palette="Set2") +
theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)

ggsave(p7,device=cairo_pdf,  filename='H2029_FPR.pdf', width=20, height=20, units=c("cm"))





fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
p2<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p3<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Group size", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p2<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p3<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Group size", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
  theme_ipsum()+labs(x='Real Total microbial loads foldchange',y='Predicted Total microbial loads foldchange',title='Total microbial loads foldchange prediction')+ylim(-1,8)+xlim(-1,8)

ggsave(p1, device=cairo_pdf, filename='Obesity_microbial_loads_foldchange_prediction_performance.pdf', width=25, height=25, units=c("cm"))


mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))

mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
p2<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + 
    scale_fill_brewer(palette="Set2")+
    geom_boxplot()+labs(x="Changed Taxa Rate", y ="MAE")+ theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
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
p3<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + geom_boxplot()+labs(x="Group Size", y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
ggsave(p3,device=cairo_pdf, filename='Obesity_MAE_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))




fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_data<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_data<-data.frame(fnr_group,fnr_data)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))

p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group))+
  geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FNR") +scale_fill_brewer(palette="Set2")  +
theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)
ggsave(p6, device=cairo_pdf, filename='Obesity_FNR.pdf', width=20, height=20, units=c("cm"))


fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_data<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_data<-data.frame(fpr_group,fpr_data)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group))+
  geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FPR") +scale_fill_brewer(palette="Set2") +
theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)

ggsave(p7,device=cairo_pdf,  filename='Obesity_FPR.pdf', width=20, height=20, units=c("cm"))





fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
p2<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p3<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Group size", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p2<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p3<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Group size", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
  theme_ipsum()+labs(x='Real Total microbial loads foldchange',y='Predicted Total microbial loads foldchange',title='Total microbial loads foldchange prediction')+ylim(-1,8)+xlim(-1,8)

ggsave(p1, device=cairo_pdf, filename='GP_microbial_loads_foldchange_prediction_performance.pdf', width=25, height=25, units=c("cm"))




mae_group<-c(rep('DR_inferred_foldchange_MAE',dim(data)[1]),rep('logged_relative_abundance_diff_MAE',dim(data)[1]),rep('qmd_logged_Density_diff_MAE',dim(data)[1]),rep('ANCOM_BC_diff_MAE',dim(data)[1]))

mae_value<-c(data$DR_inferred_foldchange_MAE,data$logged_relative_abundance_diff_MAE,data$qmd_logged_Density_diff_MAE,data$ANCOM_BC_diff_MAE)
mae_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
mae_data<-data.frame(mae_group,mae_value,mae_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
mae_data$mae_con <- cut(mae_data$mae_con, breaks = group, labels =labels)
mae_data$mae_group<-factor(mae_data$mae_group,, levels = c('qmd_logged_Density_diff_MAE','logged_relative_abundance_diff_MAE','DR_inferred_foldchange_MAE','ANCOM_BC_diff_MAE'))
p2<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + 
    scale_fill_brewer(palette="Set2")+
    geom_boxplot()+labs(x="Changed Taxa Rate", y ="MAE")+ theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
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
p3<-ggplot(mae_data, aes(x=mae_con, y=mae_value, fill=mae_group)) + geom_boxplot()+labs(x="Group Size", y ="MAE") +scale_fill_brewer(palette="Set2") + theme_ipsum()+scale_y_continuous(limits=c(0,8), breaks=seq(0,8,1))
ggsave(p3,device=cairo_pdf, filename='GP_MAE_on_groupsize.pdf', width=25, height=15, units=c("cm"))

describeBy(mae_data$mae_value,list(mae_data$mae_con,mae_data$mae_group),mat=TRUE,digits=2,quant=c(.25,.75))
describeBy(mae_data$mae_value,mae_data$mae_group,mat=TRUE,digits=2,quant=c(.25,.75))



fnr_group<-c(rep('ancomP_FNR',dim(data)[1]),rep('mwu_diff_pvalue_FNR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FNR',dim(data)[1]),rep('qmd_diff_pvalue_FNR',dim(data)[1]),rep('qmd_diff_qvalue_FNR',dim(data)[1]))
fnr_data<-c(data$ancomP_FNR,data$mwu_diff_pvalue_FNR,data$ANCOM_BC_diff_qval_FNR,data$qmd_diff_pvalue_FNR,data$qmd_diff_qvalue_FNR)
fnr_data<-data.frame(fnr_group,fnr_data)
fnr_data$fnr_group<-factor(fnr_data$fnr_group,levels = c('qmd_diff_pvalue_FNR','qmd_diff_qvalue_FNR','ancomP_FNR','ANCOM_BC_diff_qval_FNR','mwu_diff_pvalue_FNR'))

p6<-ggplot(fnr_data, aes(x=fnr_group, y=fnr_data, color=fnr_group))+
  geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FNR") +scale_fill_brewer(palette="Set2")  +
theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)
ggsave(p6, device=cairo_pdf, filename='GP_FNR.pdf', width=20, height=20, units=c("cm"))


fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_data<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_data<-data.frame(fpr_group,fpr_data)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
p7<-ggplot(fpr_data, aes(x=fpr_group, y=fpr_data, color=fpr_group))+
  geom_boxplot(outlier.colour=NA)+ geom_blank()+labs(x="", y ="FPR") +scale_fill_brewer(palette="Set2") +
theme_ipsum()+ theme(legend.position="none")+geom_jitter(size=0.4, alpha=0.9)+ylim(0,1)

ggsave(p7,device=cairo_pdf,  filename='GP_FPR.pdf', width=20, height=20, units=c("cm"))





fpr_group<-c(rep('ancomP_FPR',dim(data)[1]),rep('mwu_diff_pvalue_FPR',dim(data)[1]),rep('ANCOM_BC_diff_qval_FPR',dim(data)[1]),rep('qmd_diff_pvalue_FPR',dim(data)[1]),rep('qmd_diff_qvalue_FPR',dim(data)[1]))
fpr_value<-c(data$ancomP_FPR,data$mwu_diff_pvalue_FPR,data$ANCOM_BC_diff_qval_FPR,data$qmd_diff_pvalue_FPR,data$qmd_diff_qvalue_FPR)
fpr_con<-c(data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate,data$changedTaxaRate)
fpr_data<-data.frame(fpr_group,fpr_value,fpr_con)
group <- seq(0,1,by=0.1)
labels<-c('0~0.1','0.1~0.2','0.2~0.3','0.3~0.4','0.4~0.5','0.5~0.6','0.6~0.7','0.7~0.8','0.8~0.9','0.9~1')
fpr_data$fpr_con <- cut(fpr_data$fpr_con, breaks = group, labels =labels)
fpr_data$fpr_group<-factor(fpr_data$fpr_group,levels = c('qmd_diff_pvalue_FPR','qmd_diff_qvalue_FPR','ancomP_FPR','ANCOM_BC_diff_qval_FPR','mwu_diff_pvalue_FPR'))
p2<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p3<-ggplot(fpr_data, aes(x=fpr_con, y=fpr_value, fill=fpr_group)) + geom_boxplot()+labs(x="Group size", y ="FPR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p2<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Changed Taxa Rate", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
p3<-ggplot(fnr_data, aes(x=fnr_con, y=fnr_value, fill=fnr_group)) + geom_boxplot()+labs(x="Group size", y ="FNR") +scale_fill_brewer(palette="Set2")  + theme_ipsum()
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
