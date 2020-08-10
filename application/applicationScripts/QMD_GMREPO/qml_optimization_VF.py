
# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon,mannwhitneyu,ttest_1samp,ranksums
from sklearn.utils import resample
from ancomP.stats.ancom import ancom
from ancomP.linalg import composition
from sklearn.metrics import confusion_matrix
import copy
import sys


fileplace=sys.argv[1]
permu=int(sys.argv[2])
lp=int(sys.argv[3])
up=int(sys.argv[4])
predix=sys.argv[5]
control=sys.argv[6]
treat=sys.argv[7]


## like l1 regularization
def cal_cost(x,mod,detectionV):
    tmpmod = np.tile(mod, (x.shape[0], 1))
    tmpdetection = np.tile(detectionV, (x.shape[0], 1))
    delta = x
    tmpmod = np.abs(tmpmod + delta)
    tmpmod = tmpmod * tmpdetection
    j1 = tmpmod.sum(axis=1)
    j = j1 / tmpmod.shape[1]
    return j


genuslist = pd.read_csv(fileplace +predix+'_' + control + '_' + treat + '_taxa_IntoModel.csv', header=0, index_col=0)
genuslist['D'] = 1 / 2 * (genuslist['controlDetectionRate'] + genuslist['treatDetectionRate'])
dimensions_genus = len(genuslist)
detectionV = np.array(genuslist['D'])
mod = np.asarray(genuslist['logged_mean_diff'])
mod05 = np.asarray(genuslist['logged_mean_diff05'])
mod95 = np.asarray(genuslist['logged_mean_diff95'])

arr2 = np.array(range(2000)).reshape(2000, 1)
arr2 = arr2 / 100 - 10
res = cal_cost(arr2,mod,detectionV)
minres=np.min(res)
pdelta=arr2[res==minres]
absdelta = np.abs(pdelta)
pdelta=pdelta[np.argmin(absdelta)][0]
print(pdelta)

changRecord=pd.DataFrame(columns=['item','value','type'])
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
        x_bootstrap.append((resample(x)))
    x_bootstrap = np.mean(x_bootstrap, axis=1)
    y_bootstrap = []
    for i in range(permu):
        np.random.seed(i)
        y_bootstrap.append((resample(y)))
    y_bootstrap = np.mean(y_bootstrap, axis=1)
    differences = x_bootstrap - y_bootstrap
    lower_bound = np.percentile(differences, lp)
    upper_bound = np.percentile(differences, up)
    mean_diff=np.mean(differences)
    return mean_diff,lower_bound,upper_bound,differences


genusdata=pd.read_csv(fileplace+predix+'_'+control+'_'+treat+'_taxa_abundance_IntoModel.csv',index_col=0,header=0)

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


def stackAnalysisRes(changRecord):
    I = changRecord
    logged_relative_abundance_diff=I.loc[I['type']=='logged_relative_abundance_diff',['item','value']]
    ancomP = I.loc[I['type'] == 'ancomP',['item','value']]
    mwu_diff_pvalue = I.loc[I['type'] == 'mwu_diff_pvalue',['item','value']]
    mwu_logged_pvalue = I.loc[I['type'] == 'mwu_logged_pvalue',['item','value']]
    qml_diff_pvalue = I.loc[I['type'] == 'qml_diff_pvalue',['item','value']]
    qml_logged_loads_diff = I.loc[I['type'] == 'qml_logged_loads_diff',['item','value']]
    logged_relative_abundance_diff=logged_relative_abundance_diff.rename(columns={'value':'logged_relative_abundance_diff'})
    ancomP = ancomP.rename(
        columns={'value': 'ancomP'})
    mwu_diff_pvalue = mwu_diff_pvalue.rename(
        columns={'value': 'mwu_diff_pvalue'})
    mwu_logged_pvalue = mwu_logged_pvalue.rename(
        columns={'value': 'mwu_logged_pvalue'})
    qml_diff_pvalue = qml_diff_pvalue.rename(
        columns={'value': 'qml_diff_pvalue'})
    qml_logged_loads_diff = qml_logged_loads_diff.rename(
        columns={'value': 'qml_logged_loads_diff'})
    res=pd.merge(logged_relative_abundance_diff,ancomP,how='inner',on='item')
    res=pd.merge(res,mwu_diff_pvalue,how='inner',on='item')
    res=pd.merge(res,mwu_logged_pvalue,how='inner',on='item')
    res=pd.merge(res,qml_diff_pvalue,how='inner',on='item')
    res=pd.merge(res,qml_logged_loads_diff,how='inner',on='item')
    res=res.sort_values(by='qml_logged_loads_diff')
    res=res.reset_index(drop=True)
    return res

changRecord.to_csv(fileplace+predix+'_'+control+'_'+treat+'_analysis.csv',index=False)

stackedChangeRecord=stackAnalysisRes(changRecord)
stackedChangeRecord.to_csv(fileplace+predix+'_'+control+'_'+treat+'_analysis_stacked.csv',index=False)

