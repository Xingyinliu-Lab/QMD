##############################
##############################
##############################
##############################
##############################
# this is the batch run script for QMD analysis on GP simulation instances


import os
import time
import numpy as np
from multiprocessing import Pool
import os
import pandas as pd

# processors_num: threads of multi threading 
# this para should be modified before start this script
processors_num=20

fileplace='simulationData/GPoriData/'
permu=500
lp=5
up=95
minimum_taxa_detection_num=5
control='Control'
treat='Treat'


def cal_dataPreprocessing(cl):
    for oo in cl:
        i=oo+1
        print(i)
        predix = 'gp_'+str(i)
        if os.path.exists(fileplace + predix + '_' + control + '_' + treat + '_taxa_Into_Model.npy'):
            continue
        cmdStr='python qmd_dataPreprocessing.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(predix)+' '+control+' '+treat
        os.system(cmdStr)

def cal_optimization(cl):
    for oo in cl:
        i=oo+1
        print(i)
        predix = 'gp_'+str(i)
        if os.path.exists(fileplace+predix+'_analysis_ANCOM_ANCOM-BC_QMD.csv'):
            continue
        cmdStr='python qmd_ANCOM.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
        os.system(cmdStr)



if __name__ == '__main__':
    processor = processors_num
    cl = np.array_split(np.asarray(range(500)), processor, axis=0)
    reslist = []
    p = Pool(processor)
    for i in range(processor):
        reslist.append(p.apply_async(cal_dataPreprocessing, args=(cl[i],)))
    p.close()
    p.join()
    print('-------------------------------\n-------------------------------\n-------------------------------\n-------------------------------\n')
    print('DATA preprocessing done')
    print('Start QML analysis')
    time.sleep(5)
    reslist = []
    p = Pool(processor)
    for i in range(processor):
        reslist.append(p.apply_async(cal_optimization, args=(cl[i],)))
    p.close()
    p.join()
