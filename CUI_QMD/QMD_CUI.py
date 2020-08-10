import qmd_dataPreprocessing_cui
import qmd_optimization_cui
import sys

fileplace=sys.argv[1]# place to store analysis result
data_preprocessing_place=sys.argv[2]# place to store datapreprocessing result
datafilename=sys.argv[3]# filename of to-be-analyzed datafile
minimum_taxa_detection_num=int(sys.argv[4])#minimum taxa detection number in each group
minimum_taxa_abundance_control=float(sys.argv[5])#minimum taxa abundance median in control group
control_label=sys.argv[6]# label of control group
treat_label=sys.argv[7]# label of treatment group
group_label=sys.argv[8]# column name of the group indicator in the datafile
permu_loops=int(sys.argv[9])#permuation loops
fdr=sys.argv[10]# whether turn FDR adjustment on
plot=sys.argv[11]# whether plot the taxa changes bar chart
confidence_level=float(sys.argv[12])#confidence level to identify statistically significant differentilly abundant taxa


qmd_dataPreprocessing_cui.preprocess_qmd(data_preprocessing_place, datafilename, minimum_taxa_detection_num,
                                         control_label, treat_label,
                                         group_label, permu_loops, minimum_taxa_abundance_control)

qmd_optimization_cui.qmd_optimization(fileplace, data_preprocessing_place, control_label,
                                      treat_label, group_label, permu_loops, fdr, plot, confidence_level)

