import os
import time
import numpy as np
from multiprocessing import Pool
import os

fileplace='BarlowWork_LSI_STOOL/Data/'
permu=500
lp=5
up=95
minimum_taxa_detection_num=3
control='Control'
treat='Keto'
import pandas as pd


predix = 'LSI'
cmdStr='python qmd_dataPreprocessingJacob.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python qmd_optimizationJacob.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)



fileplace='JacobWork/Data/'
permu=500
lp=5
up=95
minimum_taxa_detection_num=3
control='Control'
treat='Keto'
import pandas as pd

predix = 'STOOL'
cmdStr='python qmd_dataPreprocessingJacob.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python qmd_optimizationJacob.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)




