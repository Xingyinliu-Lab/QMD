import numpy as np
from skbio.diversity.alpha import shannon
import pandas as pd
#skbio.diversity.alpha.shannon(counts, base=2)

fileplace_list=['E:\\github\\xingyinliu_lab\\QMD\\QMD\\Simulation\\simulationData\\GPoriData/',
                'E:\\github\\xingyinliu_lab\\QMD\\QMD\\Simulation\\simulationData\\H2029/',
                'E:\\github\\xingyinliu_lab\\QMD\\QMD\\Simulation\\simulationData\\Obesity/']
prefix_list=['gp_','H2029_','Obesity_']
for j in range(3):
    fileplace=fileplace_list[j]
    prefix=prefix_list[j]
    control_label='Control'
    treat_label='Treat'
    res=pd.DataFrame(columns={'i','control_shannon','treat_shannon'})
    res_count=0
    for q in range(1,501):
        try:
            file=prefix+str(q)+'_relative_otu_simed.csv'
            df=pd.read_csv(fileplace+file,header=0,index_col=False)
            df.fillna(value=0,inplace=True)
            c_dict = df['condition'].to_dict()
            taxa_abundance=df.drop(labels='condition',axis=1)
            control_shannon=[]
            treat_shannon=[]
            for i in taxa_abundance.index:
                try:
                    counts=(10000000*np.asarray(list(taxa_abundance.loc[i]))).astype(int)
                    shannon_value=shannon(counts)
                    label=c_dict.get(i)
                    if label==control_label:
                        control_shannon.append(shannon_value)
                    else:
                        treat_shannon.append(shannon_value)
                except:
                    print(prefix)
                    print(i)
            control_shannon=np.nanmean(control_shannon)
            treat_shannon=np.nanmean(treat_shannon)
            res.loc[res_count,'i']=q
            res.loc[res_count, 'control_shannon'] = control_shannon
            res.loc[res_count, 'treat_shannon'] = treat_shannon
            res_count=res_count+1
        except:
            pass
    res.to_csv(prefix+'diversity.csv',index=False)


