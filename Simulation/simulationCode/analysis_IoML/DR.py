##############################
##############################
##############################
##############################
##############################
#python code implementation for DR https://doi.org/10.1038/s41467-019-10656-5
#the author did not provide a standalone implementation of DR
#we extracted the implementation code from https://github.com/knightlab-analyses/reference-frames/blob/master/ipynb/simulation-benchmark.ipynb
#
# the DR author recommended to install following packages and setup a conda virtual env to run this code
# conda create -n songbird_env python=3.6 numpy=1.15.4 scikit-bio=0.5.5 seaborn pandas=0.23.4 -c conda-forge
# source activate songbird_env
# conda install tensorflow=1.10 tqdm nomkl
# conda install biom-format h5py -c conda-forge
# conda install jupyter notebook
# conda install songbird -c conda-forge


import tensorflow as tf
import numpy as np
import pandas as pd
from skbio.stats.composition import clr, clr_inv
from songbird.multinomial import MultRegression
import sys

# parameters used to specify the simulation instances to be analyzed by DR
# fileplace to store the relative abundance data
fileplace=sys.argv[1]
# predix: predix of the project  
# predix+'_relative_otu_simed.csv' is set as the filename
predix=sys.argv[2]
# control: label of control group
control=sys.argv[3]
# treat: label of treat group
treat=sys.argv[4]

# one can also set the parameters as follows if one wants to run DR on a instances
# fileplace='simulationResult/H2029/'
# control='Control'
# treat='Treat'
# predix='H2029_1'


# this function run_multinomial is the main function of DR
# this function was extracted from https://github.com/knightlab-analyses/reference-frames/blob/master/ipynb/simulation-benchmark.ipynb
def run_multinomial(table, metadata):
    model = MultRegression(
        batch_size=3, learning_rate=1e-3, beta_scale=1)
    Y = table.values
    X = metadata.values
    trainX = X[:-5]
    trainY = Y[:-5]
    testX = X[-5:]
    testY = Y[-5:]

    with tf.Graph().as_default(), tf.Session() as session:

        model(session, trainX, trainY, testX, testY)

        loss, cv, _ = model.fit(epochs=int(300))

        beta_ = clr(
            clr_inv(
                np.hstack((np.zeros((model.p, 1)), model.B))
            )
        )

    res = pd.DataFrame(
        beta_.T, columns=['intercept', 'statistic'],
        index=table.columns
    )
    return res

# read into the realtive abundance 
genusdata=pd.read_csv(fileplace+predix+'_relative_otu_simed.csv',header=0,index_col=False)
metadata=genusdata[['condition']]
genusdata['sample_id']=genusdata.index
genusdata.index='S'+ genusdata['sample_id'].astype('str')
genusdata=genusdata.drop(['condition','sample_id'],axis=1)
genusdata=genusdata.fillna(value=0)


metadata['intercept']=1
metadata['labels']=1
metadata.loc[metadata['condition']==control,'labels']=-1
metadata=metadata[['intercept', 'labels']]

rel_m1 = run_multinomial(genusdata, metadata)
#rel_a1 = ancom(rel_table1+1, metadata1['labels'])[0]
rel_m1['taxaId']=rel_m1.index

# the effect size estimated from DR is e-based log.
# transform the effect size to microbial density 

rel_m1['DR_inferred_foldchange']=np.log2(np.exp(2*rel_m1['statistic']))

rel_m1=rel_m1[['taxaId','DR_inferred_foldchange']]


rel_m1=rel_m1.reset_index(drop=True)
previous_res=pd.read_csv(fileplace+predix+'_analysis_ANCOM_ANCOM-BC_QMD_stacked.csv',header=0,index_col=False)
final_res=pd.merge(previous_res,rel_m1,how='left',on='taxaId')

previous_res2=pd.read_csv(fileplace+predix+'_analysis_ANCOM_ANCOM-BC_QMD.csv',header=0,index_col=False)
count=len(previous_res2)
for rel_m1index in rel_m1.index:
    previous_res2.loc[count,'item']=rel_m1.loc[rel_m1index,'taxaId']
    previous_res2.loc[count, 'value'] = rel_m1.loc[rel_m1index, 'DR_inferred_foldchange']
    previous_res2.loc[count, 'type'] = 'DR_inferred_foldchange'
    count=count+1

# two result are generated 
# one is in stacked table ***_analysis_stacked_vf.csv
# the other is in unstacked table  ***_analysis_vf.csv
previous_res2.to_csv(fileplace+predix+'_analysis_all_methods.csv',index_label=False)
final_res.to_csv(fileplace+predix+'_analysis_stacked_all_methods.csv',index_label=False)
