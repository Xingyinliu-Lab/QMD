
import os

fileplace='../Data/gutMicrobiomeAbundanceData/'
permu=500
lp=5
up=95
processers=10
minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='crohn'
#
print(treat)
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='crc'
#
print(treat)
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='UC'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='IBS'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='H50'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='DiabetesMellitusType2'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=30
minimum_taxa_abundance_control=0
control='H2029'
treat='Kidney'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='liver'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)


###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='autoimmuneD'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

###############

minimum_taxa_detection_num=15
minimum_taxa_abundance_control=0
control='H2029'
treat='CysticFibrosis'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='lung'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='Obesity'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='CardiovascularDiseases'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

###############

minimum_taxa_detection_num=30
minimum_taxa_abundance_control=0
control='H2029'
treat='GDM'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

###############

minimum_taxa_detection_num=50
minimum_taxa_abundance_control=0
control='H2029'
treat='hypertension'
print(treat)
#
predix = 'genus_AbundanceData'
cmdStr='python ../QMD_GMREPO/qmd_dataPreprocessing_VF_multi_thread.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(minimum_taxa_detection_num)+' '+str(minimum_taxa_abundance_control)+' '+str(predix)+' '+control+' '+treat+' '+str(processers)
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/plot_cost_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)

cmdStr='python ../QMD_GMREPO/qmd_optimization_VF.py'+' '+fileplace+' '+str(permu)+' '+str(lp)+' '+str(up)+' '+str(predix)+' '+control+' '+treat
os.system(cmdStr)
