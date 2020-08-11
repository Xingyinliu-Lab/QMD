Character User Interface QMD@XYL_Lab
=================================

## Installation recommendations
It is recommended to install CUI-QMD@XYL_Lab as follows
```python
conda create -n qmd python=3.7 numpy pandas
source activate qmd
conda install scipy
conda install statsmodels
conda install matplotlib
```
## Instructions for use
### Arguments
The QMD_CUI.py is the entrance of the software. It takes following arguments.
> 1. fileplace. The place to store analysis result.
> 2. data_preprocessing_place. The place to store datapreprocessing result
> 3. datafilename. The path of the csv file to be analyzed. Both the datafile name and its relative or absolute path should be provided.
> 4. minimum_taxa_detection_num. Minimum taxa detection samples in each group. The value should be an integer. If a taxa is only detected in a few samples less than the set value in either the control group or the treatment group, it will be filtered out. It is suggested the value is set larger than 3. If the group sample size is quite small, e.g. 5 or 6, this value can be set 2.
> 5. minimum_taxa_abundance_control. Minimum taxa median abundance in control group.The value should be in [0,1). If the median relative abundance of taxa is less than the set value, the taxa will be filtered out. This is another way to control the prevalence of taxa. One can specify this value to 0 if one do not want to filter out low abundance taxa.
> 6. control_label. Label of control group. Specify the control group label. For the demo project in this folder, it should be *Healthy*.
> 7. treat_label. Label of treatment group. Specify the treatment group label. For the demo project in this folder, it should be *CD*.
> 8. group_label. Column name of group indicator. Specify the column name of group indicator. This should be consistent with the datafile. For the demo project in this folder, it should be *condition*.
> 9. permu_loops. The software uses permutation test to make the quantification. Permutation loops affects the results’ stability. A larger permutation loop number also cost a longer time to conduct the analysis. It is suggested permutation loops is at least 500. The value should be an integer.
> 10. fdr. Whether turn FDR adjustment on. Benjamini/Hochberg FDR p value adjustment are provided in the software. One can choose turning it on or off. It is suggested if the sample size in each group are larger than 50, the FDR adjust is on. Valid values include "ON" and "OFF".
> 11. plot. Whether plot the taxa changes bar chart. One can choose whether to plot the microbial density changes diverging bars chart. Valid values include "ON" and "OFF".
> 12. confidence_level. The confidence interval for differentially abundant taxa identification. The value should be in (0.7,1). 
> 

### Outputs

QMD@XYL_Lab will output the analysis result in the specified folder. It includes

> 1. The *qmd_microbial_density_changes_plot.pdf*, the generated diverging bars chart when plot is on. 
> 2. The *QMD_result_stacked.csv*. This file contains a table list the relative abundance changes (in the column *logged_relative_abundance_diff*), the p value of QMDD(in the column *qmd_diff_pvalue*), the q value(FDR adjusted pvalue) of QMDD(in the column *qmd_diff_qvalue*), the identified differentially abundant taxa(in the column *differentially_abundant_taxa*) between group for every taxa. The taxa with YES in the *differentially_abundant_taxa* column are identified as statistically significant changed taxa. If FDR adjust is off, q value will be not provided. 
> 3. The *QMD_result.csv*, unstacked version of *QMD_result_stacked.csv*. The quantified are also recorded in this file.


## Examples
A demo project can be found in this folder. To analysis the demo project, users can using following command.
### Using absolute path in linux platform
```python
python QMD_CUI.py /mnt/d/project/CUI_QMD/demoProject /mnt/d/project/CUI_QMD/demoProject /mnt/d/project/CUI_QMD/demoProject/dc_b2.csv 2 0 Healthy CD condition 500 ON ON 0.95
```
### Using relative path in linux platform
```python
python QMD_CUI.py ./demoProject ./demoProject ./demoProject/dc_b2.csv 2 0 Healthy CD condition 500 ON ON 0.95
```
### Using absolute path in windows platform
```python
python .\QMD_CUI.py D:\project\CUI_QMD\demoProject D:\project\CUI_QMD\demoProject D:\project\CUI_QMD\demoProject\dc_b2.csv 2 0 Healthy CD condition 500 ON ON 0.95
```
### Using relative path in windows platform
```python
python .\QMD_CUI.py ./demoProject ./demoProject ./demoProject/dc_b2.csv 2 0 Healthy CD condition 500 ON ON 0.95
```

