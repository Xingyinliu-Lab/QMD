# tax detection rate
# 检查在两个组下，tax的检出率，以及是否存在明确的检出率差异

# This uses a simple normal test for proportions.
# It should be the same as running the mean z-test on the data encoded 1 for event and 0 for no event so that the sum corresponds to the count.
#
# In the one and two sample cases with two-sided alternative,
# this test produces the same p-value as proportions_chisquare, since the chisquare is the distribution of the square of a standard normal distribution
import copy

import numpy as np

from sklearn.utils import resample
import sys
import pandas as pd
fileplace=sys.argv[1]

permu=int(sys.argv[2])
lp=int(sys.argv[3])
up=int(sys.argv[4])
minimum_taxa_detection_num=int(sys.argv[5])

predix=sys.argv[6]
control=sys.argv[7]
treat=sys.argv[8]
#


#

per_res=pd.DataFrame(columns={'permu_index','permu_delta'})

for sss in range(permu):
    genusdata = pd.read_csv(fileplace + predix + '.csv', header=0)
    c_list0 = list(genusdata['condition'])
    # genusdata = copy.copy(ori_genusdata)
    np.random.seed(sss)
    c_list=list(genusdata['condition'].sample(frac=1))
    genusdata['condition']=c_list
    genusdata=genusdata.loc[(genusdata['condition']==control) | (genusdata['condition']==treat) ,]
    tmpgenuslist=list(genusdata.columns)
    genuslist=[]
    for g in tmpgenuslist:
        if g!='condition':
            genuslist.append(g)
    genuslist=pd.DataFrame(genuslist)
    genuslist.columns=['taxaId']
    # genuslist=genuslist[0:20]
    controllen=sum(genusdata['condition']==control)
    treatlen=sum(genusdata['condition']==treat)
    # print(control,controllen)
    # print(treat,treatlen)
    allcount = np.array([controllen, treatlen])

    res_detection_ratio=pd.DataFrame(columns=['taxaId','controlDetectionNum','treatDetectionNum','controlDetectionRate','treatDetectionRate'])
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
            count=count+1
        except:
            pass
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
    # 确保在bench group中，作为分子的菌的丰度大于作为分母的菌的丰度
    # 即ratio>1
    genuslist=res_detection_ratio.loc[(res_detection_ratio['controlDetectionNum'] > minimum_taxa_detection_num)
                                      & (res_detection_ratio['treatDetectionNum'] > minimum_taxa_detection_num)
    ,]
    genuslist=genuslist.reset_index(drop=True)

    genusIntoModel=list(genuslist['taxaId'])
    genusdata=genusdata[genusIntoModel+['condition']]
    sumabundance=genusdata[genusIntoModel].sum(axis=1)
    genusdata[genusIntoModel]=genusdata[genusIntoModel].div(sumabundance,axis='rows')

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


    ## like l1 regularization
    def cal_cost(x, mod, detectionV):
        tmpmod = np.tile(mod, (x.shape[0], 1))
        tmpdetection = np.tile(detectionV, (x.shape[0], 1))
        delta = x
        tmpmod = np.abs(tmpmod + delta)
        tmpmod = tmpmod * tmpdetection
        j1 = tmpmod.sum(axis=1)
        j = j1 / tmpmod.shape[1]
        return j

    genuslist['D'] = 1 / 2 * (genuslist['controlDetectionRate'] + genuslist['treatDetectionRate'])
    dimensions_genus = len(genuslist)
    detectionV = np.array(genuslist['D'])
    mod = np.asarray(genuslist['logged_mean_diff'])
    mod05 = np.asarray(genuslist['logged_mean_diff05'])
    mod95 = np.asarray(genuslist['logged_mean_diff95'])
    # gdiff=np.load(fileplace +predix+'_'+control+'_'+treat+'_taxa_DiffIntoModel.npy')
    arr2 = np.array(range(2000)).reshape(2000, 1)
    arr2 = arr2 / 100 - 10
    res = cal_cost(arr2, mod, detectionV)
    minres = np.min(res)
    pdelta = arr2[res == minres]
    absdelta = np.abs(pdelta)
    pdelta = pdelta[np.argmin(absdelta)][0]
    print(sss,pdelta)
    genusdata['condition'] = c_list0
    per_res.loc[sss,'permu_index']=sss
    per_res.loc[sss, 'permu_delta'] = pdelta
per_res.to_csv(fileplace+predix+'_total_abundance_change_permu_res.csv',index=False)