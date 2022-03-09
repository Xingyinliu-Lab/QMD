import numpy as np
import pandas as pd
import math
import powerlaw
# from scipy.stats import powerlaw
from statsmodels.distributions.empirical_distribution import ECDF
from sklearn.metrics import mean_squared_error, r2_score
from matplotlib import pyplot as plt
from scipy.stats import kstest

res=pd.DataFrame(columns={'prefix','instance','pvalue'})
count=0
for prefix in ['Obesity/Obesity','H2029/H2029','GPoriData/gp']:
    for i in range(500):
        try:
            j=i+1
            df=pd.read_csv(prefix+'_'+str(j)+'_relative_otu_simed.csv',index_col=None,header=0)
            col=df.columns
            df_list=[]
            for c in col:
                if c !='condition':
                    df_list=df_list+list(df[c])
            df_list=[xs for xs in df_list if not math.isnan(xs) ]
            df_list=np.asarray(df_list)
            np.sort(df_list)
            fit=powerlaw.Fit(df_list)
            alpha=fit.power_law.alpha
            xmin  = fit.power_law.xmin
            _,p=kstest(df_list, "powerlaw", args=(alpha, xmin), N=len(df_list))
            print(prefix,j,p)
            res.loc[count,'prefix']=prefix
            res.loc[count,'instance']=j
            res.loc[count,'pvalue']=p
            count=count+1
        except:
            pass

res.to_csv('kstest.csv')
        # a,loc,scale=powerlaw.fit(df_list)
        # predict_y=powerlaw.cdf(df_list,a,loc=loc,scale=scale)
        # ecff=ECDF(df_list)
        # true_y=ecff(df_list)
        # plt.plot(df_list,predict_y,'ro')
        # plt.show()
        # plt.plot(df_list,true_y,'ro')
        # plt.show()
        # plt.plot(true_y,predict_y,'ro')
        # plt.show()
        # print(prefix,a,loc,scale,mean_squared_error(true_y,predict_y))
        #
