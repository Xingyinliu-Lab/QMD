
import os

fileplace='VandeputteWork_B1_B2/Data/'
permu=500
lp=5
up=95
minimum_taxa_detection_num=2
minimum_taxa_abundance_control=0.0001
# processers=5
predix = 'dc_b1'

control='Healthy'
treat='CD'

cmdStr='python qmd_dataPreprocessing_Doris.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '
os.system(cmdStr)

cmdStr='python qmd_optimization_Doris.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

fileplace='DorisWork/Data/'
permu=500
lp=5
up=95
minimum_taxa_detection_num=2
minimum_taxa_abundance_control=0.0001
# processers=5
predix = 'dc_b2'

control='Healthy'
treat='CD'

cmdStr='python qmd_dataPreprocessing_Doris.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '
os.system(cmdStr)

cmdStr='python qmd_optimization_Doris.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)
