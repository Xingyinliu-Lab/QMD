import  numpy as np
import pandas as pd
import matplotlib.pyplot as plt
predix = 'genus_AbundanceData'

fileplace='../Data/gutMicrobiomeAbundanceData/'


def getRes2(control,treat):

    I = pd.read_csv(fileplace + predix + '_' + control + '_' + treat + '_analysis.csv', index_col=False, header=0)

    logged_relative_abundance_diff=I.loc[I['type']=='logged_relative_abundance_diff',['item','value']]
    ancomP = I.loc[I['type'] == 'ancomP',['item','value']]
    mwu_diff_pvalue = I.loc[I['type'] == 'mwu_diff_pvalue',['item','value']]
    mwu_logged_pvalue = I.loc[I['type'] == 'mwu_logged_pvalue',['item','value']]
    qml_diff_pvalue = I.loc[I['type'] == 'qml_diff_pvalue',['item','value']]
    qml_logged_loads_diff = I.loc[I['type'] == 'qml_logged_loads_diff',['item','value']]


    logged_relative_abundance_diff=logged_relative_abundance_diff.rename(columns={'value':'logged_relative_abundance_diff'})
    ancomP = ancomP.rename(
        columns={'value': 'ancomP'})
    mwu_diff_pvalue = mwu_diff_pvalue.rename(
        columns={'value': 'mwu_diff_pvalue'})
    mwu_logged_pvalue = mwu_logged_pvalue.rename(
        columns={'value': 'mwu_logged_pvalue'})
    qml_diff_pvalue = qml_diff_pvalue.rename(
        columns={'value': 'qml_diff_pvalue'})
    qml_logged_loads_diff = qml_logged_loads_diff.rename(
        columns={'value': 'qml_logged_loads_diff'})

    res=pd.merge(logged_relative_abundance_diff,ancomP,how='inner',on='item')
    res=pd.merge(res,mwu_diff_pvalue,how='inner',on='item')
    res=pd.merge(res,mwu_logged_pvalue,how='inner',on='item')
    res=pd.merge(res,qml_diff_pvalue,how='inner',on='item')
    res=pd.merge(res,qml_logged_loads_diff,how='inner',on='item')
    res = res.sort_values(by='qml_logged_loads_diff')
    res=res.reset_index(drop=True)
    return res

control='H2029'
treat='IBS'
res=getRes2(control,treat)
taxid=pd.read_csv(fileplace+'genus_taxaid_list.csv',index_col=0,header=0,converters={'taxa':str})
taxid['item']='taxa'+taxid['taxa']

res=pd.merge(res,taxid,how='inner',on='item')
res['taxaName']=res['genus']

res.to_csv(fileplace + predix + '_' + control + '_' + treat + '_analysis_stacked.csv',index=False)


control='H2029'
treat='hypertension'
res=getRes2(control,treat)
taxid=pd.read_csv(fileplace+'genus_taxaid_list.csv',index_col=0,header=0,converters={'taxa':str})
taxid['item']='taxa'+taxid['taxa']

res=pd.merge(res,taxid,how='inner',on='item')
res['taxaName']=res['genus']

res.to_csv(fileplace + predix + '_' + control + '_' + treat + '_analysis_stacked.csv',index=False)

