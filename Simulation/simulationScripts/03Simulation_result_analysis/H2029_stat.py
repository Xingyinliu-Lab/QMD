# the code is used to make basic statistics of 500 simualtion instances from H2029
import pandas as pd
import numpy as np
import os
fileplace='simulationData/H2029/'
stat_res=pd.DataFrame()
for oo in range(500):
    i=oo+1
    predix = 'H2029_' + str(i)
    if not os.path.exists(fileplace+predix+'_analysis_ANCOM_ANCOM-BC_QMD.csv'):
        continue
    changRecord=pd.read_csv(fileplace+predix+'_analysis_ANCOM_ANCOM-BC_QMD.csv', index_col=False,header=0)
    changRecord=changRecord[changRecord['type']!='ANCOM_BC_diff_diff_abn']
    changRecord['value'] = changRecord['value'].astype(float)

    changedTaxaCount=list(changRecord.loc[changRecord['type']=='changedTaxaCount','value'])[0]
    totalTaxaCount=list(changRecord.loc[changRecord['type']=='totalTaxaCount','value'])[0]
    MicrobiotaLoadsFoldChange=list(changRecord.loc[changRecord['type']=='MicrobiotaLoadsFoldChange','value'])[0]
    MicrobiotaLoadsFoldChange_ori = list(changRecord.loc[changRecord['type'] == 'MicrobiotaLoadsFoldChange_ori', 'value'])[0]
    changfold=np.mean(list(changRecord.loc[changRecord['type'] == 'TaxaFoldChange', 'value']))
    stat_res.loc[oo,'i']=i
    stat_res.loc[oo, 'changedTaxaCount'] = changedTaxaCount#the number of differentially abundant taxa
    stat_res.loc[oo, 'changedTaxaRate'] = changedTaxaCount/totalTaxaCount#the proportion of differentially abundant taxa
    stat_res.loc[oo, 'MicrobiotaDensityFoldChange'] = MicrobiotaLoadsFoldChange#the true/real total microbial density change
    stat_res.loc[oo, 'qmd_Density_delta'] = list(changRecord.loc[changRecord['type'] == 'qml_loads_delta', 'value'])[0]  # total microbial density change quantified by QMD
    stat_res.loc[oo, 'MicrobiotaDensityFoldChange_ori'] = MicrobiotaLoadsFoldChange_ori# the original total microbial density change before the simulated changed were adpoted to the artificial treatment group
    stat_res.loc[oo, 'changfold'] = changfold # average changes of taxa
    stat_res.loc[oo, 'ancomP_FNR'] = list(changRecord.loc[changRecord['type'] == 'ancomP_FNR', 'value'])[0]#FNR of ANCOM
    stat_res.loc[oo, 'ancomP_FPR'] = list(changRecord.loc[changRecord['type'] == 'ancomP_FPR', 'value'])[0]#FPR of ANCOM
    stat_res.loc[oo, 'mwu_diff_pvalue_FNR'] = list(changRecord.loc[changRecord['type'] == 'mwu_diff_pvalue_FNR', 'value'])[0]#FNR of MWU
    stat_res.loc[oo, 'mwu_diff_pvalue_FPR'] = list(changRecord.loc[changRecord['type'] == 'mwu_diff_pvalue_FPR', 'value'])[0]#FPR of MWU
    stat_res.loc[oo, 'mwu_logged_pvalue_FNR'] = list(changRecord.loc[changRecord['type'] == 'mwu_logged_pvalue_FNR', 'value'])[0]#FNR of MWU on logged relative abundance 
    stat_res.loc[oo, 'mwu_logged_pvalue_FPR'] = list(changRecord.loc[changRecord['type'] == 'mwu_logged_pvalue_FPR', 'value'])[0]#FPR of MWU on logged relative abundance 
    stat_res.loc[oo, 'qmd_diff_pvalue_FNR'] = list(changRecord.loc[changRecord['type'] == 'qml_diff_pvalue_FNR', 'value'])[0]#FNR of QMD
    stat_res.loc[oo, 'qmd_diff_pvalue_FPR'] = list(changRecord.loc[changRecord['type'] == 'qml_diff_pvalue_FPR', 'value'])[0]#FPR of QMD
    stat_res.loc[oo, 'qmd_diff_qvalue_FPR'] = list(changRecord.loc[changRecord['type'] == 'qml_diff_qvalue_FPR', 'value'])[0]  #FNR of QMD with fdr
    stat_res.loc[oo, 'qmd_diff_qvalue_FNR'] = list(changRecord.loc[changRecord['type'] == 'qml_diff_qvalue_FNR', 'value'])[0]  #FPR of QMD with fdr
    stat_res.loc[oo, 'ANCOM_BC_diff_qval_FPR'] =  list(changRecord.loc[changRecord['type'] == 'ANCOM_BC_diff_qval_FPR', 'value'])[0]  #FNR of ANCOM-BC with fdr
    stat_res.loc[oo, 'ANCOM_BC_diff_qval_FNR'] = list(changRecord.loc[changRecord['type'] == 'ANCOM_BC_diff_qval_FNR', 'value'])[0]   #FPR of ANCOM-BC with fdr
    stat_res.loc[oo, 'qmd_logged_Density_diff_MAE'] = list(changRecord.loc[changRecord['type'] == 'qml_logged_loads_diff_MAE', 'value'])[0]#MAE by QMD
    stat_res.loc[oo, 'logged_relative_abundance_diff_MAE'] = list(changRecord.loc[changRecord['type'] == 'logged_relative_abundance_diff_MAE', 'value'])[0]#MAE by directly using logged_relative_abundance_diff
    stat_res.loc[oo, 'ANCOM_BC_diff_MAE'] = list(changRecord.loc[changRecord['type'] == 'ANCOM_BC_diff_MAE', 'value'])[0]   #MAE of ANCOM-BC

    genusData=pd.read_csv(fileplace+predix+'_Control_Treat_taxa_abundance_IntoModel.csv', index_col=False,header=0)
    stat_res.loc[oo, 'groupsize'] =len(genusData)/2# groupsize into analysis
    ana_stacked=pd.read_csv(fileplace+predix+'_analysis_stacked_all_methods.csv', index_col=0,header=0)
    DR_inferred_foldchange=ana_stacked[['microbiota_loads_foldChange','DR_inferred_foldchange']]
    DR_inferred_foldchange=DR_inferred_foldchange.dropna()

    DR_inferred_foldchange_MAE=1/len(DR_inferred_foldchange)*np.sum(np.abs(DR_inferred_foldchange['DR_inferred_foldchange']-DR_inferred_foldchange['microbiota_loads_foldChange']))
    stat_res.loc[oo, 'DR_inferred_foldchange_MAE'] =DR_inferred_foldchange_MAE# MAE of DR

stat_res.to_csv('H2029_stat.csv',index=False)
