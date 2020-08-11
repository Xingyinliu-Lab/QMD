
##############################
##############################
##############################
##############################
##############################
#this is the code for QMD implementation and simulation comparision
#the code is slightly different from the final QMD code for for convenience of simulation
#the code is implemented on Python 3.7.7, other python version are not tested
# numpy, scikit-learn, pandas,scipy,ancomP are used in this code
# conda activate python37
# conda install scipy
# ancomP is the python implementation of ANCOM, it can be found in https://github.com/mortonjt/ancomP

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon,mannwhitneyu,ttest_1samp
from ancomP.stats.ancom import ancom
from ancomP.linalg import composition
from sklearn.metrics import confusion_matrix
import copy
import sys
from statsmodels.stats.multitest import fdrcorrection

#parameters used to specify the simulation instances to be analyzed by QMD
#fileplace to store the relative abundance data
fileplace=sys.argv[1]
#permu: permuation loops
permu=int(sys.argv[2])
#lp: lower bound of the confidence interval in the permuation
lp=int(sys.argv[3])
#up: upper bound of the confidence interval in the permuation
up=int(sys.argv[4])
# predix+'_relative_otu_simed.csv' is set as the filename
predix=sys.argv[5]
# control: label of control group
control=sys.argv[6]
# treat: label of treat group
treat=sys.argv[7]

# one can also set the parameters as follows if one wants to run QMD on a instances
# fileplace='D:\\project\\bio\\abs_abundance\\VF\\simulation\\simulationResult\H2029/'
# permu=500
# lp=5
# up=95
# minimum_taxa_detection_num=3
# control='Control'
# treat='Treat'
# predix='H2029_1'


#the optimization process 
#L1 regularization
def cal_cost(x,mod,detectionV):
    tmpmod = np.tile(mod, (x.shape[0], 1))
    tmpdetection = np.tile(detectionV, (x.shape[0], 1))
    delta = x
    tmpmod = np.abs(tmpmod + delta)
    tmpmod = tmpmod * tmpdetection
    j1 = tmpmod.sum(axis=1)
    j = j1 / tmpmod.shape[1]
    return j

# readinto the data preprocessing output
genuslist = pd.read_csv(fileplace +predix+'_' + control + '_' + treat + '_taxa_Into_Model.csv', header=0, index_col=0)
genuslist['D'] = 1 / 2 * (genuslist['controlDetectionRate'] + genuslist['treatDetectionRate'])
dimensions_genus = len(genuslist)
detectionV = np.array(genuslist['D'])
mod = np.asarray(genuslist['logged_mean_diff'])
mod05 = np.asarray(genuslist['logged_mean_diff05'])
mod95 = np.asarray(genuslist['logged_mean_diff95'])
# gdiff=np.load(fileplace +predix+'_'+control+'_'+treat+'_genusDiffIntoModel.npy')
# Traversal the possible change by step=0.01
arr2 = np.array(range(2000)).reshape(2000, 1)
arr2 = arr2 / 100 - 10
res = cal_cost(arr2,mod,detectionV)
minres=np.min(res)
pdelta=arr2[res==minres]
absdelta = np.abs(pdelta)
# get the quantified Psi
pdelta=pdelta[np.argmin(absdelta)][0]


changRecord=pd.read_csv(fileplace+predix+'_changed_record_simed.csv',index_col=False,header=0)
count=len(changRecord)
changRecord.loc[count,'item']='qml_loads_delta'
changRecord.loc[count,'value']=pdelta
changRecord.loc[count,'type']='qml_loads_delta'
count=count+1

def xydiff_confidence_interval(x,y):
    x=x[~np.isnan(x)]
    y = y[~np.isnan(y)]
    x_bootstrap = []
    for i in range(permu):
        np.random.seed(i)
        x_bootstrap.append((np.random.choice(x, size=len(x))))
    x_bootstrap = np.mean(x_bootstrap, axis=1)
    y_bootstrap = []
    for i in range(permu):
        np.random.seed(i)
        y_bootstrap.append((np.random.choice(y, size=len(y))))
    y_bootstrap = np.mean(y_bootstrap, axis=1)
    differences = x_bootstrap - y_bootstrap
    lower_bound = np.percentile(differences, lp)
    upper_bound = np.percentile(differences, up)
    mean_diff=np.mean(differences)
    return mean_diff,lower_bound,upper_bound,differences


genusdata=pd.read_csv(fileplace+predix+'_'+control+'_'+treat+'_taxa_abundance_IntoModel.csv',index_col=0,header=0)

# calculate DA result for QMDD, QMDD with fdr,MWU

taxalist=[]
for taxa in genusdata.columns:
    if taxa=='condition':
        continue
    taxalist.append(taxa)
    controlgdata = genusdata.loc[genusdata['condition'] == control, taxa]
    treatgdata = genusdata.loc[genusdata['condition'] == treat, taxa]
    controlgdata=controlgdata[controlgdata>0]
    treatgdata=treatgdata[treatgdata>0]
    loggedControl=np.log2(controlgdata)
    loggedTreat=np.log2(treatgdata)
    loggedControl=loggedControl[~np.isnan(loggedControl)]
    loggedTreat = loggedTreat[~np.isnan(loggedTreat)]
    if (len(loggedControl)<10) or (len(loggedTreat)<10) :
        mean_diff, _, _, _ = xydiff_confidence_interval(loggedTreat, loggedControl)
    else:
        mean_diff=np.mean(loggedTreat)-np.mean(loggedControl)
    changRecord.loc[count, 'item'] = taxa
    changRecord.loc[count, 'value'] = mean_diff
    changRecord.loc[count, 'type'] = 'logged_relative_abundance_diff'
    count = count + 1
    changRecord.loc[count, 'item'] = taxa
    changRecord.loc[count, 'value'] = pdelta+mean_diff
    changRecord.loc[count, 'type'] = 'qml_logged_loads_diff'
    count = count + 1
    _, qml_pvalue = mannwhitneyu(loggedControl,loggedTreat+pdelta,alternative='two-sided')
    changRecord.loc[count, 'item'] = taxa
    changRecord.loc[count, 'value'] = qml_pvalue
    changRecord.loc[count, 'type'] = 'qml_diff_pvalue'
    count = count + 1
    _, mwu_pvalue = mannwhitneyu(controlgdata,treatgdata,alternative='two-sided')
    changRecord.loc[count, 'item'] = taxa
    changRecord.loc[count, 'value'] = mwu_pvalue
    changRecord.loc[count, 'type'] = 'mwu_diff_pvalue'
    count = count + 1
    _, mwu_logged_pvalue = mannwhitneyu(loggedControl,loggedTreat,alternative='two-sided')
    changRecord.loc[count, 'item'] = taxa
    changRecord.loc[count, 'value'] = mwu_logged_pvalue
    changRecord.loc[count, 'type'] = 'mwu_logged_pvalue'
    count = count + 1


tmpPvalue=changRecord[(changRecord['type']=='qml_diff_pvalue')&(changRecord['value']>=0)]
# fdrcorrection:
# This covers Benjamini/Hochberg for independent or positively correlated and Benjamini/Yekutieli 
# for general or negatively correlated tests. Both are available in the function multipletests, as method=`fdr_bh`, resp. fdr_by.
tmp_taxlist=list(tmpPvalue['item'])
_,qvaluelist = fdrcorrection(tmpPvalue['value'], alpha=0.05, method='indep', is_sorted=False)
for asd in range(len(qvaluelist)):
    changRecord.loc[count, 'item'] = tmp_taxlist[asd]
    changRecord.loc[count, 'value'] = qvaluelist[asd]
    changRecord.loc[count, 'type'] = 'qml_diff_qvalue'
    count = count + 1

# calculate DA result for ANCOM

table=genusdata[taxalist]
table=table.fillna(value=0)
table=composition.multiplicative_replacement(table)
table=pd.DataFrame(table)
table.index=genusdata.index
table.columns=taxalist
labels, uniques = pd.factorize(genusdata['condition'])

grouping = pd.Series(labels,index=genusdata.index)
results = ancom(table, grouping, significance_test='permutative-anova', permutations=permu)
# print(results['reject'])
results['taxaId']=results.index
for j in results.index:
    changRecord.loc[count, 'item'] = results.loc[j,'taxaId']
    if results.loc[j, 'reject']==True:
        tmpvalue=1
    else:
        tmpvalue=0
    changRecord.loc[count, 'value'] = tmpvalue
    changRecord.loc[count, 'type'] = 'ancomP'
    count = count + 1

# integrate the result of ANCOM-BC
ANCOM_BC_res2=pd.read_csv(fileplace+predix+'_ANCOM_BC_res2.csv',index_col=0,header=0)
ANCOM_BC_res1=pd.read_csv(fileplace+predix+'_ANCOM_BC_res1.csv',index_col=0,header=0)
changRecord.loc[count, 'item']=ANCOM_BC_res2.loc[1,'note']
changRecord.loc[count, 'value'] = ANCOM_BC_res2.loc[1,'bias']
changRecord.loc[count, 'type'] = ANCOM_BC_res2.loc[1,'note']
count = count + 1
changRecord.loc[count, 'item']=ANCOM_BC_res2.loc[2,'note']
changRecord.loc[count, 'value'] = ANCOM_BC_res2.loc[2,'bias']
changRecord.loc[count, 'type'] = ANCOM_BC_res2.loc[2,'note']
count = count + 1
for nbv in ANCOM_BC_res1.index:
    changRecord.loc[count, 'item'] = ANCOM_BC_res1.loc[nbv, 'taxon']
    changRecord.loc[count, 'value'] = np.log2(np.exp(ANCOM_BC_res1.loc[nbv, 'mean.difference (Treat - Control)']))
    changRecord.loc[count, 'type'] = 'ANCOM_BC_diff'
    count = count + 1
    changRecord.loc[count, 'item'] = ANCOM_BC_res1.loc[nbv, 'taxon']
    changRecord.loc[count, 'value'] = ANCOM_BC_res1.loc[nbv, 'p.val']
    changRecord.loc[count, 'type'] = 'ANCOM_BC_diff_pval'
    count = count + 1
    changRecord.loc[count, 'item'] = ANCOM_BC_res1.loc[nbv, 'taxon']
    changRecord.loc[count, 'value'] = ANCOM_BC_res1.loc[nbv, 'q.val']
    changRecord.loc[count, 'type'] = 'ANCOM_BC_diff_qval'
    count = count + 1
    changRecord.loc[count, 'item'] = ANCOM_BC_res1.loc[nbv, 'taxon']
    changRecord.loc[count, 'value'] = ANCOM_BC_res1.loc[nbv, 'diff.abn']
    changRecord.loc[count, 'type'] = 'ANCOM_BC_diff_diff_abn'
    count = count + 1


#prepare data for the FPR FNR and MAE calculation
taxalist=pd.DataFrame(taxalist,columns=['taxaId'])
microbiota_loads_foldChange=changRecord.loc[changRecord['type'] == 'TaxaFoldChange',['value','item'] ]
microbiota_loads_foldChange=microbiota_loads_foldChange.rename(columns = {"value": "microbiota_loads_foldChange", "item":"taxaId"})
taxalist=pd.merge(taxalist,microbiota_loads_foldChange,how='left',on='taxaId')

logged_relative_abundance_diff=changRecord.loc[changRecord['type'] == 'logged_relative_abundance_diff',['value','item'] ]
logged_relative_abundance_diff=logged_relative_abundance_diff.rename(columns = {"value": "logged_relative_abundance_diff", "item":"taxaId"})
taxalist=pd.merge(taxalist,logged_relative_abundance_diff,how='left',on='taxaId')

qml_logged_loads_diff=changRecord.loc[changRecord['type'] == 'qml_logged_loads_diff',['value','item'] ]
qml_logged_loads_diff=qml_logged_loads_diff.rename(columns = {"value": "qml_logged_loads_diff", "item":"taxaId"})
taxalist=pd.merge(taxalist,qml_logged_loads_diff,how='left',on='taxaId')

# true
ranksums_on_loads=changRecord.loc[changRecord['type'] == 'ranksums_on_loads',['value','item'] ]
ranksums_on_loads=ranksums_on_loads.rename(columns = {"value": "ranksums_on_loads", "item":"taxaId"})
taxalist=pd.merge(taxalist,ranksums_on_loads,how='left',on='taxaId')

# directly use mwu on relative abundance
mwu_diff_pvalue=changRecord.loc[changRecord['type'] == 'mwu_diff_pvalue',['value','item'] ]
mwu_diff_pvalue=mwu_diff_pvalue.rename(columns = {"value": "mwu_diff_pvalue", "item":"taxaId"})
taxalist=pd.merge(taxalist,mwu_diff_pvalue,how='left',on='taxaId')

# use mwu on logged relative abundance
mwu_logged_pvalue=changRecord.loc[changRecord['type'] == 'mwu_logged_pvalue',['value','item'] ]
mwu_logged_pvalue=mwu_logged_pvalue.rename(columns = {"value": "mwu_logged_pvalue", "item":"taxaId"})
taxalist=pd.merge(taxalist,mwu_logged_pvalue,how='left',on='taxaId')
# use qml
qml_diff_pvalue=changRecord.loc[changRecord['type'] == 'qml_diff_pvalue',['value','item'] ]
qml_diff_pvalue=qml_diff_pvalue.rename(columns = {"value": "qml_diff_pvalue", "item":"taxaId"})
taxalist=pd.merge(taxalist,qml_diff_pvalue,how='left',on='taxaId')
# use ANCOM
ancomP=changRecord.loc[changRecord['type'] == 'ancomP',['value','item'] ]
ancomP=ancomP.rename(columns = {"value": "ancomP", "item":"taxaId"})
taxalist=pd.merge(taxalist,ancomP,how='left',on='taxaId')

qml_diff_qvalue=changRecord.loc[changRecord['type'] == 'qml_diff_qvalue',['value','item'] ]
qml_diff_qvalue=qml_diff_qvalue.rename(columns = {"value": "qml_diff_qvalue", "item":"taxaId"})
taxalist=pd.merge(taxalist,qml_diff_qvalue,how='left',on='taxaId')

ANCOM_BC_diff_qval=changRecord.loc[changRecord['type'] == 'ANCOM_BC_diff_qval',['value','item'] ]
ANCOM_BC_diff_qval=ANCOM_BC_diff_qval.rename(columns = {"value": "ANCOM_BC_diff_qval", "item":"taxaId"})
taxalist=pd.merge(taxalist,ANCOM_BC_diff_qval,how='left',on='taxaId')

ANCOM_BC_diff=changRecord.loc[changRecord['type'] == 'ANCOM_BC_diff',['value','item'] ]
ANCOM_BC_diff=ANCOM_BC_diff.rename(columns = {"value": "ANCOM_BC_diff", "item":"taxaId"})
taxalist=pd.merge(taxalist,ANCOM_BC_diff,how='left',on='taxaId')

values = {'microbiota_loads_foldChange':0}
taxalist=taxalist.fillna(value=values)

# MAE calculation
taxaNumber=len(taxalist)

logged_relative_abundance_diff_MAEdata=taxalist[['microbiota_loads_foldChange','logged_relative_abundance_diff']]
logged_relative_abundance_diff_MAEdata=logged_relative_abundance_diff_MAEdata.dropna()

qml_logged_loads_diff_MAEdata=taxalist[['microbiota_loads_foldChange','qml_logged_loads_diff']]
qml_logged_loads_diff_MAEdata=qml_logged_loads_diff_MAEdata.dropna()

ANCOM_BC_diff_MAEdata=taxalist[['microbiota_loads_foldChange','ANCOM_BC_diff']]
ANCOM_BC_diff_MAEdata=ANCOM_BC_diff_MAEdata.dropna()

logged_relative_abundance_diff_MAE=1/len(logged_relative_abundance_diff_MAEdata)*np.sum(np.abs(logged_relative_abundance_diff_MAEdata['logged_relative_abundance_diff']-logged_relative_abundance_diff_MAEdata['microbiota_loads_foldChange']))
qml_logged_loads_diff_MAE=1/len(qml_logged_loads_diff_MAEdata)*np.sum(np.abs(qml_logged_loads_diff_MAEdata['qml_logged_loads_diff']-qml_logged_loads_diff_MAEdata['microbiota_loads_foldChange']))
ANCOM_BC_diff_MAE=1/len(ANCOM_BC_diff_MAEdata)*np.sum(np.abs(ANCOM_BC_diff_MAEdata['ANCOM_BC_diff']-ANCOM_BC_diff_MAEdata['microbiota_loads_foldChange']))

changRecord.loc[count, 'item'] = 'logged_relative_abundance_diff_MAE'
changRecord.loc[count, 'value'] = logged_relative_abundance_diff_MAE
changRecord.loc[count, 'type'] = 'logged_relative_abundance_diff_MAE'
count = count + 1
changRecord.loc[count, 'item'] = 'qml_logged_loads_diff_MAE'
changRecord.loc[count, 'value'] = qml_logged_loads_diff_MAE
changRecord.loc[count, 'type'] = 'qml_logged_loads_diff_MAE'
count = count + 1
changRecord.loc[count, 'item'] = 'ANCOM_BC_diff_MAE'
changRecord.loc[count, 'value'] = ANCOM_BC_diff_MAE
changRecord.loc[count, 'type'] = 'ANCOM_BC_diff_MAE'
count = count + 1

# FPR FNR calculation
tmptaxalist=copy.deepcopy(taxalist)

trueTreat=np.array(tmptaxalist.loc[:,'ranksums_on_loads'])
trueTreat[trueTreat>0.9]=0.9
trueTreat[trueTreat<0.05]=1
trueTreat[trueTreat<1]=0

mwu_diff_pvalue=np.array(tmptaxalist.loc[:,'mwu_diff_pvalue'])
mwu_diff_pvalue[mwu_diff_pvalue>0.9]=0.9
mwu_diff_pvalue[mwu_diff_pvalue<0.05]=1
mwu_diff_pvalue[mwu_diff_pvalue<1]=0

mwu_logged_pvalue=np.array(tmptaxalist.loc[:,'mwu_logged_pvalue'])
mwu_logged_pvalue[mwu_logged_pvalue>0.9]=0.9
mwu_logged_pvalue[mwu_logged_pvalue<0.05]=1
mwu_logged_pvalue[mwu_logged_pvalue<1]=0


qml_diff_pvalue=np.array(tmptaxalist.loc[:,'qml_diff_pvalue'])
qml_diff_pvalue[qml_diff_pvalue>0.9]=0.9
qml_diff_pvalue[qml_diff_pvalue<0.05]=1
qml_diff_pvalue[qml_diff_pvalue<1]=0

ancomP=np.array(tmptaxalist.loc[:,'ancomP'])


qml_diff_qvalue=np.array(tmptaxalist.loc[:,'qml_diff_qvalue'])
qml_diff_qvalue[qml_diff_qvalue>0.9]=0.9
qml_diff_qvalue[qml_diff_qvalue<0.05]=1
qml_diff_qvalue[qml_diff_qvalue<1]=0

ANCOM_BC_diff_qval=np.array(tmptaxalist.loc[:,'ANCOM_BC_diff_qval'])
ANCOM_BC_diff_qval[ANCOM_BC_diff_qval>0.9]=0.9
ANCOM_BC_diff_qval[ANCOM_BC_diff_qval<0.05]=1
ANCOM_BC_diff_qval[ANCOM_BC_diff_qval<1]=0
# False positive rate (α) = type I error = 1 − specificity = FP / (FP + TN) = 180 / (180 + 1820) = 9%
# False negative rate (β) = type II error = 1 − sensitivity = FN / (TP + FN) = 10 / (20 + 10) = 33%


tn, fp, fn, tp = confusion_matrix(trueTreat[mwu_diff_pvalue>-1].astype(np.int),mwu_diff_pvalue[mwu_diff_pvalue>-1].astype(np.int)).ravel()
FPR=fp/(fp+tn)
FNR=fn/(tp+fn)
changRecord.loc[count, 'item'] = 'mwu_diff_pvalue_FPR'
changRecord.loc[count, 'value'] = FPR
changRecord.loc[count, 'type'] = 'mwu_diff_pvalue_FPR'
count = count + 1
changRecord.loc[count, 'item'] = 'mwu_diff_pvalue_FNR'
changRecord.loc[count, 'value'] = FNR
changRecord.loc[count, 'type'] = 'mwu_diff_pvalue_FNR'
count = count + 1


tn, fp, fn, tp = confusion_matrix(trueTreat[mwu_logged_pvalue>-1].astype(np.int), mwu_logged_pvalue[mwu_logged_pvalue>-1].astype(np.int)).ravel()
FPR=fp/(fp+tn)
FNR=fn/(tp+fn)
changRecord.loc[count, 'item'] = 'mwu_logged_pvalue_FPR'
changRecord.loc[count, 'value'] = FPR
changRecord.loc[count, 'type'] = 'mwu_logged_pvalue_FPR'
count = count + 1
changRecord.loc[count, 'item'] = 'mwu_logged_pvalue_FNR'
changRecord.loc[count, 'value'] = FNR
changRecord.loc[count, 'type'] = 'mwu_logged_pvalue_FNR'
count = count + 1

tn, fp, fn, tp = confusion_matrix(trueTreat[qml_diff_pvalue>-1].astype(np.int), qml_diff_pvalue[qml_diff_pvalue>-1].astype(np.int)).ravel()
FPR=fp/(fp+tn)
FNR=fn/(tp+fn)
changRecord.loc[count, 'item'] = 'qml_diff_pvalue_FPR'
changRecord.loc[count, 'value'] = FPR
changRecord.loc[count, 'type'] = 'qml_diff_pvalue_FPR'
count = count + 1
changRecord.loc[count, 'item'] = 'qml_diff_pvalue_FNR'
changRecord.loc[count, 'value'] = FNR
changRecord.loc[count, 'type'] = 'qml_diff_pvalue_FNR'
count = count + 1

tn, fp, fn, tp = confusion_matrix(trueTreat[ancomP>-1].astype(np.int), ancomP[ancomP>-1].astype(np.int)).ravel()
FPR=fp/(fp+tn)
FNR=fn/(tp+fn)
changRecord.loc[count, 'item'] = 'ancomP_FPR'
changRecord.loc[count, 'value'] = FPR
changRecord.loc[count, 'type'] = 'ancomP_FPR'
count = count + 1
changRecord.loc[count, 'item'] = 'ancomP_FNR'
changRecord.loc[count, 'value'] = FNR
changRecord.loc[count, 'type'] = 'ancomP_FNR'
count = count + 1


tn, fp, fn, tp = confusion_matrix(trueTreat[qml_diff_qvalue>-1].astype(np.int), qml_diff_qvalue[qml_diff_qvalue>-1].astype(np.int)).ravel()
FPR=fp/(fp+tn)
FNR=fn/(tp+fn)
changRecord.loc[count, 'item'] = 'qml_diff_qvalue_FPR'
changRecord.loc[count, 'value'] = FPR
changRecord.loc[count, 'type'] = 'qml_diff_qvalue_FPR'
count = count + 1
changRecord.loc[count, 'item'] = 'qml_diff_qvalue_FNR'
changRecord.loc[count, 'value'] = FNR
changRecord.loc[count, 'type'] = 'qml_diff_qvalue_FNR'
count = count + 1

tn, fp, fn, tp = confusion_matrix(trueTreat[ANCOM_BC_diff_qval>-1].astype(np.int), ANCOM_BC_diff_qval[ANCOM_BC_diff_qval>-1].astype(np.int)).ravel()
FPR=fp/(fp+tn)
FNR=fn/(tp+fn)
changRecord.loc[count, 'item'] = 'ANCOM_BC_diff_qval_FPR'
changRecord.loc[count, 'value'] = FPR
changRecord.loc[count, 'type'] = 'ANCOM_BC_diff_qval_FPR'
count = count + 1
changRecord.loc[count, 'item'] = 'ANCOM_BC_diff_qval_FNR'
changRecord.loc[count, 'value'] = FNR
changRecord.loc[count, 'type'] = 'ANCOM_BC_diff_qval_FNR'
count = count + 1
# FNR
# FPR

# two result are generated 
# one is in stacked table ***_analysis_stacked.csv
# the other is in unstacked table  ***_analysis.csv
# the DR related result are not include in these two file. Use 'batchrun_DR.py' to integrate the DR result
# 
changRecord.to_csv(fileplace+predix+'_analysis_ANCOM_ANCOM-BC_QMD.csv',index=False)
taxalist.to_csv(fileplace+predix+'_analysis_ANCOM_ANCOM-BC_QMD_stacked.csv',index=False)

