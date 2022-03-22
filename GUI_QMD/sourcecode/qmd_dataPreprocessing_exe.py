# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd


def preprocess_qmd(data_preprocessing_place, datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, os_p,manner):
    # print(manner)
    permu = int(permu_loops)
    lp = 5
    up = 95
    minimum_taxa_detection_num = int(minimum_taxa_detection_num)
    minimum_taxa_abundance_control = float(minimum_taxa_abundance_control)
    control = control_label
    treat = treat_label
    genusdata = pd.read_csv(datafilename, header=0, index_col=False)
    genusdata = genusdata.loc[(genusdata[group_label] == control) | (
        genusdata[group_label] == treat), ]
    tmpgenuslist = list(genusdata.columns)
    genuslist = []
    for g in tmpgenuslist:
        if g != group_label:
            genuslist.append(g)
    sumabundance = genusdata[genuslist].sum(axis=1)
    genusdata[genuslist] = genusdata[genuslist].div(sumabundance, axis='rows')
    genuslist = pd.DataFrame(genuslist)
    genuslist.columns = ['taxaId']
    controllen = sum(genusdata[group_label] == control)
    treatlen = sum(genusdata[group_label] == treat)
    res_detection_ratio = pd.DataFrame(columns=[
                                       'taxaId', 'controlDetectionNum', 'treatDetectionNum', 'controlDetectionRate', 'treatDetectionRate', 'controlAbundance'])
    count = 0
    for i in genuslist.index:
        gid = genuslist.loc[i, 'taxaId']
        # print(i, gid)
        try:
            controlgdata = genusdata.loc[genusdata[group_label]
                                         == control, gid]
            treatgdata = genusdata.loc[genusdata[group_label] == treat, gid]
            controlobs = sum(controlgdata > 0)
            treatobs = sum(treatgdata > 0)
            res_detection_ratio.loc[count, 'taxaId'] = gid
            res_detection_ratio.loc[count, 'controlDetectionNum'] = controlobs
            res_detection_ratio.loc[count, 'treatDetectionNum'] = treatobs
            res_detection_ratio.loc[count,
                                    'controlDetectionRate'] = controlobs/controllen
            res_detection_ratio.loc[count,
                                    'treatDetectionRate'] = treatobs/treatlen
            if controlobs > 1:
                res_detection_ratio.loc[count, 'controlAbundance'] = np.nanmedian(
                    controlgdata)
            else:
                res_detection_ratio.loc[count, 'controlAbundance'] = 0
            count = count+1
        except:
            pass

    res_detection_ratio.to_csv(
        data_preprocessing_place+'/taxa_detectionRate_info.csv')

    def xdiff_confidence_interval(x):
        x = x[~np.isnan(x)]
        x_bootstrap = []
        for i in range(permu):
            np.random.seed(i)
            x_bootstrap.append((np.random.choice(x, size=len(x))))
        x_bootstrap = np.mean(x_bootstrap, axis=1)
        lower_bound = np.percentile(x_bootstrap, lp)
        upper_bound = np.percentile(x_bootstrap, up)
        mean_x = np.mean(x_bootstrap)
        return mean_x, lower_bound, upper_bound




    genuslist = res_detection_ratio.loc[(res_detection_ratio['controlDetectionNum'] > minimum_taxa_detection_num)
                                        & (res_detection_ratio['treatDetectionNum'] > minimum_taxa_detection_num)
                                        & (res_detection_ratio['controlAbundance'] >= minimum_taxa_abundance_control), ]

    genuslist = genuslist.reset_index(drop=True)
    genusIntoModel = list(genuslist['taxaId'])

    filtered_genuslist=res_detection_ratio[~res_detection_ratio['taxaId'].isin(genusIntoModel)]

    if manner=='Drop directly' or len(filtered_genuslist)==0:
        genusdata = genusdata[genusIntoModel+[group_label]]
    if manner=='Summarize as others':
        if len(filtered_genuslist)>0:
            filtered_genuslist=list(filtered_genuslist['taxaId'])
            print(filtered_genuslist)
            genusdata['others'] = genusdata[filtered_genuslist].apply(lambda x: x.sum(), axis=1)
            genusIntoModel=genusIntoModel+['others']
            genusdata = genusdata[genusIntoModel+[group_label]]
            genuslist.loc[count,'taxaId']='others'
            controlgdata = genusdata.loc[genusdata[group_label]== control, 'others']
            treatgdata = genusdata.loc[genusdata[group_label] == treat, 'others']
            controlobs = sum(controlgdata > 0)
            treatobs = sum(treatgdata > 0)
            genuslist.loc[count, 'controlDetectionNum'] = controlobs
            genuslist.loc[count, 'treatDetectionNum'] = treatobs
            genuslist.loc[count,
                                    'controlDetectionRate'] = controlobs/controllen
            genuslist.loc[count,
                                    'treatDetectionRate'] = treatobs/treatlen
            if controlobs > 1:
                genuslist.loc[count, 'controlAbundance'] = np.median(
                    controlgdata)
            else:
                genuslist.loc[count, 'controlAbundance'] = 0

    genusIntoModel = list(genuslist['taxaId'])
    genusdata = genusdata[genusIntoModel+[group_label]]
    sumabundance = genusdata[genusIntoModel].sum(axis=1)
    genusdata[genusIntoModel] = genusdata[genusIntoModel].div(
        sumabundance, axis='rows')
    genusdata.to_csv(data_preprocessing_place+'/abundance_data_intoModel.csv')



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
    for i in genuslist.index:
        gid = genuslist.loc[i, 'taxaId']
        try:
            controlgdata = genusdata.loc[genusdata[group_label]
                                         == control, gid]
            treatgdata = genusdata.loc[genusdata[group_label] == treat, gid]
            controlgdata = controlgdata[controlgdata > 0]
            treatgdata = treatgdata[treatgdata > 0]
            mean_t, _, _ = xdiff_confidence_interval(np.log2(treatgdata))
            mean_c, _, _ = xdiff_confidence_interval(np.log2(controlgdata))
            mean_diff, lower_bound, upper_bound, _ = xydiff_confidence_interval(
                np.log2(treatgdata), np.log2(controlgdata))
            genuslist.loc[i, 'treat_logged_mean'] = mean_t
            genuslist.loc[i, 'control_logged_mean'] = mean_c
            genuslist.loc[i, 'logged_mean_diff'] = mean_diff
            genuslist.loc[i, 'logged_mean_diff05'] = lower_bound
            genuslist.loc[i, 'logged_mean_diff95'] = upper_bound
        except:
            pass
    genuslist['ID'] = genuslist.index
    genuslist = genuslist.reset_index(drop=True)
    genuslist.to_csv(data_preprocessing_place+'/taxa_Into_Model.csv')

#
# preprocess_qmd('E:\\pyqt5\\QMD_v2\\Project\\af1\\datapreprocessing/', 'E:\\pyqt5\\QMD_v2\\demoData\\dc_b2.csv', 5, 'Healthy', 'CD', 'condition',
#                500,0, '','Summarize as others')