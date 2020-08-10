# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import fdrcorrection

def qmd_optimization(fileplace, data_preprocessing_place, control_label, treat_label, group_label, permu_loops, fdr, plot, confidence_level):
    fileplace = fileplace
    permu = int(permu_loops)
    lp = 5
    up = 95
    control = control_label
    treat = treat_label
    confidence_level=float(confidence_level)
    def cal_cost(x, mod, detectionV):
        tmpmod = np.tile(mod, (x.shape[0], 1))
        tmpdetection = np.tile(detectionV, (x.shape[0], 1))
        delta = x
        tmpmod = np.abs(tmpmod + delta)
        tmpmod = tmpmod * tmpdetection
        j1 = tmpmod.sum(axis=1)
        j = j1 / tmpmod.shape[1]
        return j

    genuslist = pd.read_csv(data_preprocessing_place +
                            '/taxa_Into_Model.csv', header=0, index_col=0)
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
    # print(pdelta)

    changRecord = pd.DataFrame(columns=['item', 'value', 'type'])
    count = len(changRecord)
    changRecord.loc[count, 'item'] = 'qmd_total_density_change'
    changRecord.loc[count, 'value'] = pdelta
    changRecord.loc[count, 'type'] = 'qmd_total_density_change'
    count = count+1

    def xydiff_confidence_interval(x, y):
        x = x[~np.isnan(x)]
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
        mean_diff = np.mean(differences)
        return mean_diff, lower_bound, upper_bound, differences

    genusdata = pd.read_csv(
        data_preprocessing_place+'/taxa_abundance_Into_Model.csv', index_col=0, header=0)

    taxalist = []
    for taxa in genusdata.columns:
        if taxa == group_label:
            continue
        taxalist.append(taxa)
        controlgdata = genusdata.loc[genusdata[group_label] == control, taxa]
        treatgdata = genusdata.loc[genusdata[group_label] == treat, taxa]
        controlgdata = controlgdata[controlgdata > 0]
        treatgdata = treatgdata[treatgdata > 0]
        loggedControl = np.log2(controlgdata)
        loggedTreat = np.log2(treatgdata)
        loggedControl = loggedControl[~np.isnan(loggedControl)]
        loggedTreat = loggedTreat[~np.isnan(loggedTreat)]
        mean_diff, _, _, _ = xydiff_confidence_interval(
            loggedTreat, loggedControl)
        changRecord.loc[count, 'item'] = taxa
        changRecord.loc[count, 'value'] = mean_diff
        changRecord.loc[count, 'type'] = 'logged_relative_abundance_diff'
        count = count + 1
        changRecord.loc[count, 'item'] = taxa
        changRecord.loc[count, 'value'] = pdelta+mean_diff
        changRecord.loc[count, 'type'] = 'qmd_density_diff'
        count = count + 1
        _, qmd_pvalue = mannwhitneyu(
            loggedControl, loggedTreat+pdelta, alternative='two-sided')
        changRecord.loc[count, 'item'] = taxa
        changRecord.loc[count, 'value'] = qmd_pvalue
        changRecord.loc[count, 'type'] = 'qmd_diff_pvalue'
        count = count + 1

    tmpPvalue = changRecord[(changRecord['type'] == 'qmd_diff_pvalue') & (
        changRecord['value'] >= 0)]
    # fdrcorrection:
    # This covers Benjamini/Hochberg for independent or positively correlated and Benjamini/Yekutieli
    # for general or negatively correlated tests. Both are available in the function multipletests, as method=`fdr_bh`, resp. fdr_by.
    tmp_taxlist = list(tmpPvalue['item'])
    _, qvaluelist = fdrcorrection(tmpPvalue['value'], alpha=0.05, method='indep', is_sorted=False)
    for asd in range(len(qvaluelist)):
        changRecord.loc[count, 'item'] = tmp_taxlist[asd]
        changRecord.loc[count, 'value'] = qvaluelist[asd]
        changRecord.loc[count, 'type'] = 'qmd_diff_qvalue'
        count = count + 1

    def stackAnalysisRes(changRecord):
        I = changRecord
        logged_relative_abundance_diff = I.loc[I['type'] ==
                                               'logged_relative_abundance_diff', ['item', 'value']]
        qmd_diff_pvalue = I.loc[I['type'] ==
                                'qmd_diff_pvalue', ['item', 'value']]
        qmd_density_diff = I.loc[I['type'] ==
                                 'qmd_density_diff', ['item', 'value']]
        qmd_diff_qvalue = I.loc[I['type'] ==
                                'qmd_diff_qvalue', ['item', 'value']]
        logged_relative_abundance_diff = logged_relative_abundance_diff.rename(
            columns={'value': 'logged_relative_abundance_diff'})
        qmd_diff_pvalue = qmd_diff_pvalue.rename(
            columns={'value': 'qmd_diff_pvalue'})
        qmd_density_diff = qmd_density_diff.rename(
            columns={'value': 'qmd_density_diff'})
        qmd_diff_qvalue = qmd_diff_qvalue.rename(
            columns={'value': 'qmd_diff_qvalue'})
        res = pd.merge(logged_relative_abundance_diff,
                       qmd_diff_pvalue, how='inner', on='item')
        res = pd.merge(res, qmd_density_diff, how='inner', on='item')
        res = pd.merge(res, qmd_diff_qvalue, how='inner', on='item')
        res = res.sort_values(by='qmd_density_diff')
        res = res.reset_index(drop=True)
        return res

    stackedChangeRecord = stackAnalysisRes(changRecord)
    stackedChangeRecord['differentially_abundant_taxa'] = 'No'
    if fdr == 'ON':
        stackedChangeRecord.loc[stackedChangeRecord['qmd_diff_qvalue'] < (
            1-confidence_level), 'differentially_abundant_taxa'] = 'Yes'
    if fdr == 'OFF':
        stackedChangeRecord.loc[stackedChangeRecord['qmd_diff_pvalue'] < (
            1-confidence_level), 'differentially_abundant_taxa'] = 'Yes'
        stackedChangeRecord = stackedChangeRecord.drop(
            ['qmd_diff_qvalue'], axis=1)
    for ads in stackedChangeRecord.index:
        changRecord.loc[count, 'item'] = stackedChangeRecord.loc[ads, 'item']
        if stackedChangeRecord.loc[ads, 'differentially_abundant_taxa'] == 'No':
            changRecord.loc[count, 'value'] = 0
        if stackedChangeRecord.loc[ads, 'differentially_abundant_taxa'] == 'Yes':
            changRecord.loc[count, 'value'] = 1
        changRecord.loc[count, 'type'] = 'QMDD'
        count = count + 1

    changRecord.to_csv(fileplace+'/QMD_result.csv', index=False)
    stackedChangeRecord.to_csv(
        fileplace+'/QMD_result_stacked.csv', index=False)

    if plot == 'ON':
        stackedChangeRecord['colors'] = [
            'green' if x < 0 else 'red' for x in stackedChangeRecord['qmd_density_diff']]
        stackedChangeRecord.sort_values('qmd_density_diff', inplace=True)
        stackedChangeRecord.reset_index(inplace=True)
        taxacount=len(stackedChangeRecord)
        taxalength=0
        for adf in stackedChangeRecord.index:
            tmpcharlength=len(stackedChangeRecord.loc[adf,'item'])
            if tmpcharlength>taxalength:
                taxalength=tmpcharlength
        # print(taxacount,taxalength)
        plt.figure(figsize=(np.max([taxalength/30,6]),np.max([taxacount/6,6])))
        plt.hlines(y=stackedChangeRecord.index, xmin=0, xmax=stackedChangeRecord.qmd_density_diff,
                   color=stackedChangeRecord.colors, alpha=0.4, linewidth=5)
        plt.yticks(stackedChangeRecord.index,
                   stackedChangeRecord.item, fontsize=8)
        plt.title('Microbial density changes', fontdict={'size': 12})
        plt.grid(linestyle='--', alpha=0.5)
        # plt.tight_layout()
        plt.savefig(fileplace+'/qmd_microbial_density_changes_plot.pdf',
                    dpi=600, bbox_inches='tight')
