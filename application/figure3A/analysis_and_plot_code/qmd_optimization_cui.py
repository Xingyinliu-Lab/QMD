# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

def qmd_optimization(genuslist, predix):

    def cal_cost(x, mod, detectionV):
        tmpmod = np.tile(mod, (x.shape[0], 1))
        tmpdetection = np.tile(detectionV, (x.shape[0], 1))
        delta = x
        tmpmod = np.abs(tmpmod + delta)
        tmpmod = tmpmod * tmpdetection
        j1 = tmpmod.sum(axis=1)
        j = j1 / tmpmod.shape[1]
        return j

    genuslist = pd.read_csv(genuslist, header=0, index_col=0)
    genuslist['D'] = 1 / 2 * \
        (genuslist['controlDetectionRate'] + genuslist['treatDetectionRate'])
    detectionV = np.array(genuslist['D'])
    mod = np.asarray(genuslist['logged_mean_diff'])
    arr2 = np.array(range(2000)).reshape(2000, 1)
    arr2 = arr2 / 100 - 10
    res = cal_cost(arr2, mod, detectionV)
    minres = np.min(res)
    pdelta = arr2[res == minres]
    absdelta = np.abs(pdelta)
    pdelta = pdelta[np.argmin(absdelta)][0]
    return pdelta
