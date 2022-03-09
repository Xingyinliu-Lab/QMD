from qmd_dataPreprocessing_cui import preprocess_qmd
from qmd_optimization_cui import qmd_optimization
import sys
control_label='Health'
group_label='group'
permu_loops=500
manner=1
minimum_taxa_abundance_control=0
minimum_taxa_detection_num=5


predix='PRJNA763023_paired'
datafilename=predix+'.csv'
treat_label='yCRC(huadong cohort)'
control_label='yControl(huadong cohort)'
preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA763023_paired'
datafilename=predix+'.csv'
treat_label='oCRC(huadong cohort)'
control_label='oControl(huadong cohort)'
preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA763023_single'
datafilename=predix+'.csv'
treat_label='oCRC(Fudan cohort)'
control_label='oControl(Fudan cohort)'
preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA763023_single'
datafilename=predix+'.csv'
treat_label='yCRC(Fudan cohort)'
control_label='yControl(Fudan cohort)'
preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

predix='PRJNA433269'
datafilename=predix+'.csv'
treat_label='obesity'
control_label='normal'
preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJEB6070'
datafilename=predix+'.csv'
treat_label='CRC'
control_label='Normal'
preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

predix='PRJNA397219'
datafilename=predix+'.csv'
treat_label='Colorectal Neoplasms'
control_label='Health'
preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)





control_label='Health'

predix='PRJNA268708'
datafilename=predix+'.csv'
treat_label='Irritable Bowel Syndrome'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

predix='PRJEB11419'
datafilename=predix+'.csv'
treat_label='Irritable Bowel Syndrome'
control_label='healthy'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

###################
control_label='Health'
predix='PRJEB1220'
datafilename=predix+'.csv'
treat_label='Colitis, Ulcerative'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJEB1220'
datafilename=predix+'.csv'
treat_label='Crohn Disease'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJEB4336'
datafilename=predix+'.csv'
treat_label='Obesity'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJEB7949'
datafilename=predix+'.csv'
treat_label='Inflammatory Bowel Diseases'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

predix='PRJNA232056'
datafilename=predix+'.csv'
treat_label='Crohn Disease'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA282013'
datafilename=predix+'.csv'
treat_label='Autistic Disorder'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

predix='PRJNA289586'
datafilename=predix+'.csv'
treat_label='Diabetes Mellitus, Type 1'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA353587'
datafilename=predix+'.csv'
treat_label='Clostridium difficile infection'
control_label='healthy control for CDI'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

predix='PRJNA355023'
datafilename=predix+'.csv'
treat_label='Autistic Disorder'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA368966'
datafilename=predix+'.csv'
treat_label='Inflammatory Bowel Diseases'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA368966'
datafilename=predix+'.csv'
treat_label='Colitis, Ulcerative'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)

predix='PRJNA379979'
datafilename=predix+'.csv'
treat_label='Clostridium difficile colitis'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)



predix='PRJNA385949'
datafilename=predix+'.csv'
treat_label='Inflammatory Bowel Diseases'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA388210'
datafilename=predix+'.csv'
treat_label='Colitis, Ulcerative'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA389280'
datafilename=predix+'.csv'
treat_label='Crohn Disease'
control_label='nonIBD'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)



predix='PRJNA389280'
datafilename=predix+'.csv'
treat_label='Colitis, Ulcerative'
control_label='nonIBD'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)





predix='PRJNA397219'
datafilename=predix+'.csv'
treat_label='Colorectal Neoplasms'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)



predix='PRJNA401981'
datafilename=predix+'.csv'
treat_label='Obesity'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)




predix='PRJNA422434'
datafilename=predix+'.csv'
treat_label='Diabetes Mellitus, Type 2'
control_label='Health'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)



predix='PRJNA428898'
datafilename=predix+'.csv'
treat_label='Crohn Disease'
control_label='Healthy'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)



predix='PRJNA438164'
datafilename=predix+'.csv'
treat_label='ulcerative colitis'
control_label='Healthy'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA445640'
datafilename=predix+'.csv'
treat_label='CRC'
control_label='Healthy'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)


predix='PRJNA448494'
datafilename=predix+'.csv'
treat_label='T2D'
control_label='Healthy'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)





predix='PRJNA453829'
datafilename=predix+'.csv'
treat_label='CRC'
control_label='Healthy'

preprocess_qmd(predix,datafilename, minimum_taxa_detection_num, control_label, treat_label, group_label, permu_loops, minimum_taxa_abundance_control, manner)
pdelta=qmd_optimization(predix+'_taxa_Into_Model.csv', predix)
print(predix,treat_label,pdelta)