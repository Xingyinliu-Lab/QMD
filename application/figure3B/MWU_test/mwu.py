import pandas as pd
from mwutest import mannwhitneyu_np
import numpy as np
genusdata=pd.read_csv('PRJNA763023_single.csv',index_col=None,header=0)
taxalist=list(genusdata.columns)
taxalist.remove('group')

sumabundance = genusdata[taxalist].sum(axis=1)
genusdata[taxalist] = genusdata[taxalist].div(
    sumabundance, axis='rows')

for t in taxalist:
    control=genusdata.loc[(genusdata['group']=='yControl(Fudan cohort)')&(genusdata[t]>0),t]
    treat=genusdata.loc[(genusdata['group']=='yCRC(Fudan cohort)')&(genusdata[t]>0),t]
    if len(control)>3 and len(treat)>3:
        _,p=mannwhitneyu_np(np.log2(control),np.log2(treat))
        print(t,p)




