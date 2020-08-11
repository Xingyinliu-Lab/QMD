##############################
##############################
##############################
##############################
##############################
# this is the batch run script for DR analysis on all simulation instances
import os
import time
import numpy as np
from multiprocessing import Pool
import os
import pandas as pd

control='Control'
treat='Treat'
# processors_num: threads of multi threading 
# this para should be modified before start this script
processors_num=20


def cal_H2029(cl):
    fileplace='simulationData/H2029/'
    for oo in cl:
        i=oo+1
        print(i)
        predix = 'H2029_'+str(i)
        if os.path.exists(fileplace+predix+'_analysis_stacked_all_methods.csv'):
            continue
        cmdStr='python DR.py'+' '+fileplace+' '+str(predix)+' '+control+' '+treat
        os.system(cmdStr)

def cal_Obesity(cl):
    fileplace='simulationData/Obesity/'
    for oo in cl:
        i=oo+1
        print(i)
        predix = 'Obesity_'+str(i)
        if os.path.exists(fileplace+predix+'_analysis_stacked_all_methods.csv'):
            continue
        cmdStr='python DR.py'+' '+fileplace+' '+str(predix)+' '+control+' '+treat
        os.system(cmdStr)


def cal_gp(cl):
    fileplace='simulationData/GPoriData/'
    for oo in cl:
        i=oo+1
        print(i)
        predix = 'gp_'+str(i)
        if os.path.exists(fileplace+predix+'_analysis_stacked_all_methods.csv'):
            continue
        cmdStr='python DR.py'+' '+fileplace+' '+str(predix)+' '+control+' '+treat
        os.system(cmdStr)



if __name__ == '__main__':
    
    processor = processors_num
    cl = np.array_split(np.asarray(range(500)), processor, axis=0)
    reslist = []
    p = Pool(processor)
    for i in range(processor):
        reslist.append(p.apply_async(cal_H2029, args=(cl[i],)))
    p.close()
    p.join()
    print('-------------------------------\n-------------------------------\n-------------------------------\n-------------------------------\n')
    print('H2029 done')
    time.sleep(5)
    reslist = []
    p = Pool(processor)
    for i in range(processor):
        reslist.append(p.apply_async(cal_Obesity, args=(cl[i],)))
    p.close()
    p.join()
    print('-------------------------------\n-------------------------------\n-------------------------------\n-------------------------------\n')
    print('Obesity done')
    time.sleep(5)
    reslist = []
    p = Pool(processor)
    for i in range(processor):
        reslist.append(p.apply_async(cal_gp, args=(cl[i],)))
    p.close()
    p.join()
