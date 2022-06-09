# the code is used to make basic statistics of 500 simualtion instances from GP
import pandas as pd
import numpy as np
import os
fileplace='T:\\sd\\GPoriData/'
stat_res=pd.DataFrame()
control='Control'
treat='Treat'

for oo in range(500):
    i=oo+1
    predix = 'gp_' + str(i)
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

    stat_res.loc[oo, 'MicrobiotaDensityFoldChange'] = MicrobiotaLoadsFoldChange#the true/real total microbial density change
    stat_res.loc[oo, 'qmd_Density_delta'] = list(changRecord.loc[changRecord['type'] == 'qml_loads_delta', 'value'])[0]  # total microbial density change quantified by QMD

    permu_df=pd.read_csv(fileplace+predix+'_'+control+'_'+treat+'_total_abundance_change_permu_res.csv', index_col=False,header=0)
    permu_list=np.array(permu_df['permu_delta'])
    stat_res.loc[oo, 'permutation_quantitle'] = sum(permu_list<stat_res.loc[oo, 'qmd_Density_delta'])/len(permu_list)
    stat_res.loc[oo, 'permutation_pvalue'] = 1
    stat_res.loc[oo, 'permutation_test'] = False

    if stat_res.loc[oo, 'qmd_Density_delta']>=np.mean(permu_list):
        stat_res.loc[oo, 'permutation_pvalue'] = 1-stat_res.loc[oo, 'permutation_quantitle']
    else:
        stat_res.loc[oo, 'permutation_pvalue'] = stat_res.loc[oo, 'permutation_quantitle']

    if stat_res.loc[oo, 'permutation_pvalue']<=0.05:
        stat_res.loc[oo, 'permutation_test'] = True

stat_res.to_csv('GP_Stat_on_qmd_permutation_on_total_abundance_change_stat.csv',index=False)


fileplace='X:\\sd\\H2029/'
stat_res=pd.DataFrame()
control='Control'
treat='Treat'

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

    stat_res.loc[oo, 'MicrobiotaDensityFoldChange'] = MicrobiotaLoadsFoldChange#the true/real total microbial density change
    stat_res.loc[oo, 'qmd_Density_delta'] = list(changRecord.loc[changRecord['type'] == 'qml_loads_delta', 'value'])[0]  # total microbial density change quantified by QMD

    permu_df=pd.read_csv(fileplace+predix+'_'+control+'_'+treat+'_total_abundance_change_permu_res.csv', index_col=False,header=0)
    permu_list=np.array(permu_df['permu_delta'])
    stat_res.loc[oo, 'permutation_quantitle'] = sum(permu_list<stat_res.loc[oo, 'qmd_Density_delta'])/len(permu_list)
    stat_res.loc[oo, 'permutation_pvalue'] = 1
    stat_res.loc[oo, 'permutation_test'] = False

    if stat_res.loc[oo, 'qmd_Density_delta']>=np.mean(permu_list):
        stat_res.loc[oo, 'permutation_pvalue'] = 1-stat_res.loc[oo, 'permutation_quantitle']
    else:
        stat_res.loc[oo, 'permutation_pvalue'] = stat_res.loc[oo, 'permutation_quantitle']

    if stat_res.loc[oo, 'permutation_pvalue']<=0.05:
        stat_res.loc[oo, 'permutation_test'] = True

stat_res.to_csv('H2029_Stat_on_qmd_permutation_on_total_abundance_change_stat.csv',index=False)



fileplace='X:\\sd\\Obesity/'
stat_res=pd.DataFrame()
control='Control'
treat='Treat'

for oo in range(500):
    i=oo+1
    predix = 'Obesity_' + str(i)
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

    stat_res.loc[oo, 'MicrobiotaDensityFoldChange'] = MicrobiotaLoadsFoldChange#the true/real total microbial density change
    stat_res.loc[oo, 'qmd_Density_delta'] = list(changRecord.loc[changRecord['type'] == 'qml_loads_delta', 'value'])[0]  # total microbial density change quantified by QMD

    permu_df=pd.read_csv(fileplace+predix+'_'+control+'_'+treat+'_total_abundance_change_permu_res.csv', index_col=False,header=0)
    permu_list=np.array(permu_df['permu_delta'])
    stat_res.loc[oo, 'permutation_quantitle'] = sum(permu_list<stat_res.loc[oo, 'qmd_Density_delta'])/len(permu_list)
    stat_res.loc[oo, 'permutation_pvalue'] = 1
    stat_res.loc[oo, 'permutation_test'] = False

    if stat_res.loc[oo, 'qmd_Density_delta']>=np.mean(permu_list):
        stat_res.loc[oo, 'permutation_pvalue'] = 1-stat_res.loc[oo, 'permutation_quantitle']
    else:
        stat_res.loc[oo, 'permutation_pvalue'] = stat_res.loc[oo, 'permutation_quantitle']

    if stat_res.loc[oo, 'permutation_pvalue']<=0.05:
        stat_res.loc[oo, 'permutation_test'] = True

stat_res.to_csv('Obesity_Stat_on_qmd_permutation_on_total_abundance_change_stat.csv',index=False)
