##############################
##############################
##############################
##############################
##############################
#this is the code for data preprocess of QMD
#the code is slightly different from the final QMD code for for convenience of simulation
#the code is implemented on Python 3.7.7, other python version are not tested
# numpy, scikit-learn, pandas are used in this code
# conda create -n python37 python=3.7 numpy pandas
# conda activate python37
# conda install scikit-learn



import numpy as np
import sys
import pandas as pd


#parameters used to specify the simulation instances to be analyzed by QMD
#fileplace to store the relative abundance data
fileplace=sys.argv[1]
#permu: permuation loops
permu=int(sys.argv[2])
#lp: lower bound of the confidence interval in the permuation
lp=int(sys.argv[3])
#up: upper bound of the confidence interval in the permuation
up=int(sys.argv[4])
#minimum_taxa_detection_num: the minimum detection number of taxa
#taxa whose total detection_samples_number is lower than minimum_taxa_detection_num will not go into the next step analysis
minimum_taxa_detection_num=int(sys.argv[5])
# predix+'_relative_otu_simed.csv' is set as the filename
predix=sys.argv[6]
# control: label of control group
control=sys.argv[7]
# treat: label of treat group
treat=sys.argv[8]

# one can also set the parameters as follows if one wants to run QMD data preprocessing on a instances
# fileplace = 'simulationResult/H2029/'
# permu=1000
# lp=5
# up=95
# minimum_taxa_detection_num=10
# predix='vG'
# control='G'
# treat='FG'

#read into the relative abundance data 
genusdata=pd.read_csv(fileplace+predix+'_relative_otu_simed.csv',header=0)
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
controllen=sum(genusdata['condition']==control)
treatlen=sum(genusdata['condition']==treat)
allcount = np.array([controllen, treatlen])
res_detection_ratio=pd.DataFrame(columns=['taxaId','controlDetectionNum','treatDetectionNum','controlDetectionRate','treatDetectionRate'])
# find the detection rate of taxa
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
# store the detection rate info into ***_all_genes_detectionRate_info
res_detection_ratio.to_csv(fileplace+predix+'_'+control+'_'+treat+'_all_taxa_detectionRate_info.csv')

# using permutation get the mean and the confidence bounds
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

# filter out taxa with low detection num
genuslist=res_detection_ratio.loc[(res_detection_ratio['controlDetectionNum'] > minimum_taxa_detection_num)
                                  & (res_detection_ratio['treatDetectionNum'] > minimum_taxa_detection_num)
,]
genuslist=genuslist.reset_index(drop=True)

genusIntoModel=list(genuslist['taxaId'])
genusdata=genusdata[genusIntoModel+['condition']]
sumabundance=genusdata[genusIntoModel].sum(axis=1)

genusdata[genusIntoModel]=genusdata[genusIntoModel].div(sumabundance,axis='rows')

# regenerate new relative abundance data with high-detected taxa
genusdata.to_csv(fileplace+predix+'_'+control+'_'+treat+'_taxa_abundance_IntoModel.csv')

# using permutation get the mean diff and the confidence bounds
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

# store the data preprocessing result for QMD 
genuslist['ID']=genuslist.index
genuslist.to_csv(fileplace+predix+'_'+control+'_'+treat+'_taxa_Into_Model.csv')
# np.save(fileplace+predix+'_'+control+'_'+treat+'_genusDiffIntoModel.npy',gdiff)


