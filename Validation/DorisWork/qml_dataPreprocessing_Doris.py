# tax detection rate
# 检查在两个组下，tax的检出率，以及是否存在明确的检出率差异

# This uses a simple normal test for proportions.
# It should be the same as running the mean z-test on the data encoded 1 for event and 0 for no event so that the sum corresponds to the count.
#
# In the one and two sample cases with two-sided alternative,
# this test produces the same p-value as proportions_chisquare, since the chisquare is the distribution of the square of a standard normal distribution
import numpy as np

from sklearn.utils import resample
import sys
import pandas as pd
fileplace=sys.argv[1]

permu=int(sys.argv[2])
lp=int(sys.argv[3])
up=int(sys.argv[4])
minimum_taxa_detection_num=int(sys.argv[5])
minimum_taxa_abundance_control=float(sys.argv[6])
predix=sys.argv[7]
control=sys.argv[8]
treat=sys.argv[9]

genusdata=pd.read_csv(fileplace+predix+'.csv',header=0)
genusdata=genusdata.loc[(genusdata['condition']==control) | (genusdata['condition']==treat) ,]
tmpgenuslist=list(genusdata.columns)
genuslist=[]
for g in tmpgenuslist:
    if g!='condition':
        genuslist.append(g)
sumabundance=genusdata[genuslist].sum(axis=1)
genusdata[genuslist]=genusdata[genuslist].div(sumabundance,axis='rows')

genuslist=pd.DataFrame(genuslist)
genuslist.columns=['taxaId']
# genuslist=genuslist[0:20]

controllen=sum(genusdata['condition']==control)
treatlen=sum(genusdata['condition']==treat)
# print(control,controllen)
# print(treat,treatlen)
allcount = np.array([controllen, treatlen])

res_detection_ratio=pd.DataFrame(columns=['taxaId','controlDetectionNum','treatDetectionNum','controlDetectionRate','treatDetectionRate','controlAbundance'])
count=0
for i in genuslist.index:
    gid=genuslist.loc[i,'taxaId']
    # print(i, gid)
    try:
        controlgdata=genusdata.loc[genusdata['condition']==control,gid]
        treatgdata = genusdata.loc[genusdata['condition'] == treat, gid]
        controlobs=sum(controlgdata>0)
        treatobs = sum(treatgdata > 0)
        nobs = np.array([controlobs, treatobs])
        res_detection_ratio.loc[count, 'taxaId'] = gid
        res_detection_ratio.loc[count, 'controlDetectionNum'] = controlobs
        res_detection_ratio.loc[count, 'treatDetectionNum'] = treatobs
        res_detection_ratio.loc[count, 'controlDetectionRate'] = controlobs/controllen
        res_detection_ratio.loc[count,'treatDetectionRate']=treatobs/treatlen
        if controlobs > 1:
            res_detection_ratio.loc[count, 'controlAbundance'] = np.median(controlgdata)
        else:
            res_detection_ratio.loc[count, 'controlAbundance'] = 0
        count=count+1
    except:
        pass

res_detection_ratio.to_csv(fileplace+predix+'_'+control+'_'+treat+'_all_taxa_detectionRate_info.csv')


def xdiff_confidence_interval(x):
    x=x[~np.isnan(x)]
    x_bootstrap = []
    for i in range(permu):
        np.random.seed(i)
        x_bootstrap.append((np.random.choice(x, size=len(x))))
    x_bootstrap = np.mean(x_bootstrap, axis=1)
    lower_bound = np.percentile(x_bootstrap, lp)
    upper_bound = np.percentile(x_bootstrap, up)
    mean_x=np.mean(x_bootstrap)
    return mean_x,lower_bound,upper_bound

genuslist=res_detection_ratio.loc[(res_detection_ratio['controlDetectionNum'] > minimum_taxa_detection_num)
                                  & (res_detection_ratio['treatDetectionNum'] > minimum_taxa_detection_num)
                                  &(res_detection_ratio['controlAbundance'] > minimum_taxa_abundance_control)
,]

genuslist=genuslist.reset_index(drop=True)
genusIntoModel=list(genuslist['taxaId'])
genusdata=genusdata[genusIntoModel+['condition']]
sumabundance=genusdata[genusIntoModel].sum(axis=1)
genusdata[genusIntoModel]=genusdata[genusIntoModel].div(sumabundance,axis='rows')

genusdata.to_csv(fileplace+predix+'_'+control+'_'+treat+'_taxa_abundance_IntoModel.csv')

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
gdiff=[]
for i in genuslist.index:
    gid=genuslist.loc[i,'taxaId']
    taxlist = []
    # print(i, gid)
    try:
        controlgdata = genusdata.loc[genusdata['condition'] == control, gid]
        treatgdata = genusdata.loc[genusdata['condition'] == treat, gid]
        controlobs = sum(controlgdata > 0)
        treatobs = sum(treatgdata > 0)
        controlgdata=controlgdata[controlgdata>0]
        treatgdata=treatgdata[treatgdata>0]
        mean_t, _,_=xdiff_confidence_interval(np.log2(treatgdata))
        mean_c, _, _ = xdiff_confidence_interval(np.log2(controlgdata))
        mean_diff, lower_bound, upper_bound,diff = xydiff_confidence_interval(np.log2(treatgdata), np.log2(controlgdata))
        diff=diff.reshape(-1,1)
        genuslist.loc[i, 'treat_logged_mean'] = mean_t
        genuslist.loc[i, 'control_logged_mean'] = mean_c
        genuslist.loc[i, 'logged_mean_diff'] = mean_diff
        genuslist.loc[i, 'logged_mean_diff05'] = lower_bound
        genuslist.loc[i, 'logged_mean_diff95'] = upper_bound
        if len(gdiff)>0:
            gdiff=np.hstack((gdiff,diff))
        else:
            gdiff=diff
    except:
        pass

genuslist['ID']=genuslist.index
genuslist.to_csv(fileplace+predix+'_'+control+'_'+treat+'_taxa_IntoModel.csv')
# np.save(fileplace+predix+'_'+control+'_'+treat+'_taxa_DiffIntoModel.npy',gdiff)


