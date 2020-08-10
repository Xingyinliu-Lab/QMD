## the code is to prepare simulation instances from Obesity series
# pandas, numpy, scipy are needed


import pandas as pd
import numpy as np
from scipy.stats import wilcoxon, ranksums

fileplace='originalDataPool/Obesity/'
minimun_taxa_relative_abundance=0.001
minimun_taxa_relative_abundance_obs=10



for ooo in range(500):
    i=ooo+1
    oriData = pd.read_csv(fileplace + 'Obesity' + '_relative_otu.csv', header=0, index_col=0)
    oriData = oriData.reset_index(drop=True)
     # filter out low detection rate taxa
    taxalist = []
    for taxa in oriData.columns:
        obs = oriData.loc[:, taxa]
        obsRate = sum(obs > minimun_taxa_relative_abundance)
        if obsRate >= minimun_taxa_relative_abundance_obs:
            taxalist.append(taxa)
    oriData = oriData[taxalist]
    
    # treat the taxa with low abundance as not detected
    oriData[oriData < minimun_taxa_relative_abundance] = 0
    # microbiota density = relative abundance *  total microbial density
    # generate microbiota density table
    np.random.seed(i)
    seed2 = int(np.random.uniform(low=2000, high=3000, size=1)[0])
    for j in oriData.index:
        np.random.seed(i*j*seed2)
        coef = np.random.uniform(low=1e8, high=1e9, size=1)[0]
        oriData.loc[j,:]=oriData.loc[j,:]*coef

    # shuffle oriData @ axis=0
    oriData=oriData.sample(frac=1).reset_index(drop=True)
    # generate random groupsize
    np.random.seed(i+20000)
    groupSize=int(np.random.uniform(low=6,high=100,size=1)[0])
    controlDataset = oriData.iloc[:groupSize, ]
    treatDataset = oriData.iloc[groupSize:(2*groupSize), ]
    # generate random Proportion of differentially abundant taxa 
    newcolumns=np.asarray(oriData.columns).reshape(-1)
    np.random.seed(i+30000)
    foldChangeProportion=np.random.uniform(low=0.05,high=0.95,size=1)[0]
    # choose taxa to be changed 
    totalTaxaNumber=len(newcolumns)
    changedTaxaNumber=int(totalTaxaNumber*foldChangeProportion)
    np.random.seed(i+40000)
    changedTaxa=np.random.choice(newcolumns, size=changedTaxaNumber, replace=False, p=None)

    changRecord=pd.DataFrame(columns=['item','value','type'])
    count=0

    sumControlLoads = controlDataset.sum(axis=1)
    sumTreatLoads = treatDataset.sum(axis=1)
    changRecord.loc[count, 'item'] = 'MicrobiotaLoads_ori'
    changRecord.loc[count, 'value'] = np.mean(np.log2(sumTreatLoads)) - np.mean(np.log2(sumControlLoads))
    changRecord.loc[count, 'type'] = 'MicrobiotaLoadsFoldChange_ori'
    count = count + 1

    np.random.seed(i)
    seed2 = int(np.random.uniform(low=300, high=2000, size=1)[0])
    for taxa in changedTaxa:
        np.random.seed(seed2 * count)
        # change the taxa with random foldchange
        taxaFoldChange=float(np.random.uniform(low=-10,high=10,size=1)[0])
        effectsize=np.power(2,taxaFoldChange)
        treatDataset.loc[:,taxa]=treatDataset.loc[:,taxa]*effectsize
        changRecord.loc[count,'item']=taxa
        changRecord.loc[count, 'value'] = taxaFoldChange
        changRecord.loc[count, 'type'] = 'TaxaFoldChange'
        count=count+1

    changRecord.loc[count, 'item'] = 'changedTaxaCount'
    changRecord.loc[count, 'value'] = count
    changRecord.loc[count, 'type'] = 'changedTaxaCount'
    count=count+1

    changRecord.loc[count, 'item'] = 'totalTaxaCount'
    changRecord.loc[count, 'value'] = totalTaxaNumber
    changRecord.loc[count, 'type'] = 'totalTaxaCount'
    count=count+1

    sumControlLoads=controlDataset.sum(axis=1)
    sumTreatLoads=treatDataset.sum(axis=1)
    changRecord.loc[count, 'item'] = 'MicrobiotaLoads'
    changRecord.loc[count, 'value'] = np.mean(np.log2(sumTreatLoads))-np.mean(np.log2(sumControlLoads))
    changRecord.loc[count, 'type'] = 'MicrobiotaLoadsFoldChange'
    count=count+1
    controlDataset.loc[:,'condition']='Control'
    treatDataset.loc[:,'condition']='Treat'

    loadsDataset=pd.concat([controlDataset,treatDataset],axis=0)
    taxalist=list(np.asarray(oriData.columns).reshape(-1))
    for taxa in taxalist:
        controld=loadsDataset.loc[loadsDataset.loc[:,'condition']=='Control',taxa]
        treatd=loadsDataset.loc[loadsDataset.loc[:,'condition']=='Treat',taxa]
        controld=np.asarray(controld).reshape(-1)
        controld=controld[controld>0]
        controld=np.log2(controld)
        treatd=np.asarray(treatd).reshape(-1)
        treatd=treatd[treatd>0]
        treatd=np.log2(treatd)
        _,p=ranksums(controld,treatd)
        changRecord.loc[count, 'item'] = taxa
        changRecord.loc[count, 'value'] = p
        changRecord.loc[count, 'type'] = 'ranksums_on_loads'
        count=count+1
    # absolute abundance
    absolute_abundance_Dataset = pd.concat([controlDataset, treatDataset], axis=0)
    absolute_abundance_Dataset.to_csv(fileplace + 'Obesity_' + str(i) + '_absolute_abundance_simed.csv', index=False)
    # relative abundance
    controlDataset.loc[:,taxalist]=controlDataset.loc[:,taxalist].div(sumControlLoads,axis='rows')
    treatDataset.loc[:,taxalist]=treatDataset.loc[:,taxalist].div(sumTreatLoads,axis='rows')
    simDataset=pd.concat([controlDataset,treatDataset],axis=0)
    simDataset.to_csv(fileplace+'Obesity_'+str(i)+'_relative_otu_simed.csv',index=False)
    # store change record
    changRecord.to_csv(fileplace+'Obesity_'+str(i)+'_changed_record_simed.csv',index=False)
