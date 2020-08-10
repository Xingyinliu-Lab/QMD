
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
gdiff=np.load(fileplace +predix+'_'+control+'_'+treat+'_taxa_DiffIntoModel.npy')
arr2 = np.array(range(2000)).reshape(2000, 1)
arr2 = arr2 / 100 - 10
# arr2=np.array([0.78])
res = cal_cost(arr2,mod,detectionV)
plt.plot(arr2, res)
plt.show()
# res = cal_cost2(arr2,gdiff,detectionV)
minres=np.min(res)
pdelta=arr2[res==minres]
absdelta = np.abs(pdelta)
pdelta=pdelta[np.argmin(absdelta)][0]
print(predix,pdelta)
arr2 = arr2.reshape(-1)
res = res.reshape(-1)
res = pd.DataFrame([arr2, res])
res.to_csv(fileplace+predix+'_'+control+'_'+treat+'_cost.csv')