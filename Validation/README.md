QMD application on human gut microbiome on common diseases compared with healthy population
=================
## Folder structure
```
Validation # LSI STOOL B1 B2 Validation code and results
  ├─BarlowWork_LSI_STOOL # LSI STOOL validation
  │  ├─Data # input and output data in LSI STOOL validation
  │  │      LSI.csv # relative abundance data
  │  │      LSI_absolute_abundance.csv # absoulte abundance data 
  │  │      LSI_ANCOM_BC_res1.csv # ANCOM-BC bias estimation result
  │  │      LSI_ANCOM_BC_res2.csv # ANCOM-BC microbial density estimation and differentially abundant taxa (DA) identification result
  │  │      LSI_compare.csv # ANCOM-BC, QMD, DR and other methods performance comparision on microbial density estimation
  │  │      LSI_compare.xlsx # excel version of LSI_compare.csv
  │  │      LSI_Control_Keto_all_taxa_detectionRate_info.csv # taxa detection rate info for LSI
  │  │      LSI_Control_Keto_analysis.csv # analysis reslut of LSI validation
  │  │      LSI_Control_Keto_analysis_stacked.csv # stacked analysis reslut of LSI validation
  │  │      LSI_Control_Keto_taxa_abundance_IntoModel.csv # relative abundance of remained taxa after removing low detection rate taxa
  │  │      LSI_Control_Keto_taxa_IntoModel.csv # Remained taxa after removing low detection rate taxa
  │  │      LSI_cost.csv #  the optimization process records
  │  │      LSI_cost.xlsx # xlsx version of the optimization process records
  │  │      LSI_DR_analysis.csv # DR analysis result
  │  │      LSI_loads.csv # microbial density of LSI
  │  │      LSI_res_forPlot.csv # plot data for LSI
  |  |      ...# STOOL validation input and output result
  │  │      
  │  ├─figures # Generate figure 2 and extended figure 1.
  │  │      LSI_compare.pdf # figure 2 c
  │  │      LSI_cost.pdf # figure 2 b
  │  │      LSI_cost_magnifier.pdf # figure 2 b
  │  │      LSI_foldchange_validation.pdf # figure 2 a
  │  │      STOOL_compare.pdf # extended figure 1 c
  │  │      STOOL_cost.pdf # extended figure 1 b
  |  |      STOOL_cost_magnifier.pdf # extended figure 1 b
  │  │      STOOL_foldchange_validation.pdf # extended figure 1 a
  │  │      
  │  ├─LSI_relative_abundance_compare_with_density_data # Generate figure 1.a
  │  │      LSI.csv # LSI relative abundance data
  │  │      LSI_density_control.pdf # figure 1 a
  │  │      LSI_density_keto.pdf # figure 1 b 
  │  │      LSI_density.csv # LSI microbial density table
  │  │      LSI_relative_abundance_control.pdf # figure 1 c
  │  │      LSI_relative_abundance_keto.pdf # figure 1 d
  │  │      plot_.R # scripts for figure 1 a, b, c, d plotting
  │  │      
  │  └─LSI_STOOL_analysis_code # validation scripts
  │          ANCOM_BC.R # ANCOM_BC implementation 
  │          DR_Jacob.py # DR implementation
  │          plot.R # chart plotting scripts for figure 2 and extended figure 1.
  │          qmd_dataPreprocessingJacob.py # datapreprocessing scripts for LSI and STOOL
  │          QMD_Jacob.py # QMD implementation entrance
  │          qmd_optimizationJacob.py # QMD scripts for LSI and STOOL
  │          
  └─VandeputteWork_B1_B2 # B1 B2 validation
      ├─B1_B2_analysis_code # validation scripts
      │      ANCOM_BC.R # ANCOM_BC implementation 
      │      DR_Doris.py # DR implementation
      │      plot.R # chart plotting scripts for extended figure 2.
      │      QMD_Doris.py # datapreprocessing scripts for B1 and B2
      │      qmd_dataPreprocessing_Doris.py # QMD implementation entrance
      │      qmd_optimization_Doris.py # QMD scripts for B1 and B2.
      │      
      ├─Data # input and output data in B1 B2 validation
      │      dc_b1.csv # relative abundance data
      │      dc_b1_absolute_abundance.csv # absoulte abundance data
      │      dc_b1_ANCOM_BC_res1.csv # ANCOM-BC bias estimation result
      │      dc_b1_ANCOM_BC_res2.csv # ANCOM-BC microbial density estimation and differentially abundant taxa (DA) identification result
      │      dc_b1_compare.csv # ANCOM-BC, QMD, DR and other methods performance comparision on microbial density estimation
      │      dc_b1_compare.xlsx # excel version of B1_compare.csv
      │      dc_b1_DR_analysis.csv # DR analysis result
      │      dc_b1_Healthy_CD_all_taxa_detectionRate_info.csv # taxa detection rate info for B1
      │      dc_b1_Healthy_CD_analysis.csv # analysis reslut of B1 validation
      │      dc_b1_Healthy_CD_analysis_stacked.csv # stacked analysis reslut of B1 validation
      │      dc_b1_Healthy_CD_taxa_abundance_IntoModel.csv # relative abundance of remained taxa after removing low detection rate taxa
      │      dc_b1_Healthy_CD_taxa_IntoModel.csv # Remained taxa after removing low detection rate taxa
      │      dc_b1_loads.csv # microbial density of B1
      |      ...
      │      
      └─figures # Generate extended figure 2.
              dc_b1_compare.pdf # extended figure 2 a.
              dc_b2_compare.pdf # extended figure 2 a.

```
## Software dependencies

R version 3.6.1 (2019-07-05) was used to conduct the chart plotting. Following librarys are needed.
```
library(forcats)
library(ggplot2)
library(dplyr)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(psych)
library(pheatmap)
library(nloptr)
library(dplyr)
```
Note 1: R code implementation of [ANCOM-BC](https://doi.org/10.1038/s41467-020-17041-7). The code is a copy from [the author's github](https://github.com/FrederickHuangLin/ANCOM-BC/blob/master/scripts/ancom_bc_v1.0.R)

Python 3.7.7 was used to conduct the QMD analysis processes. Following packages are needed.
```
numpy
pandas
matplotlib
scipy
sklearn
ancomP
multiprocessing
tensorflow
skbio
songbird
```
Note 2: ancomP is the python implementation of ANCOM, it can be found in https://github.com/mortonjt/ancomP
Note 3: Python code implementation for [DR](https://doi.org/10.1038/s41467-019-10656-5). The author did not provide a standalone implementation of DR. We extracted the implementation code from [the simulatin pynote from the author](https://github.com/knightlab-analyses/reference-frames/blob/master/ipynb/simulation-benchmark.ipynb). The author recommended to install following packages and setup a conda virtual env to run the DR code
```
conda create -n songbird_env python=3.6 numpy=1.15.4 scikit-bio=0.5.5 seaborn pandas=0.23.4 -c conda-forge
source activate songbird_env
conda install tensorflow=1.10 tqdm nomkl
conda install biom-format h5py -c conda-forge
conda install jupyter notebook
conda install songbird -c conda-forge
```



## Reproduce the analysis results

To reproduce the analysis results of LSI and STOOL validation, please follows the steps below:
> 1. Set the workspace of R to the fileplace storing abundance data of all analyzed cases, ANCOM_BC.R Line 427.
> 2. Run ANCOM_BC.R.
> 3. Go to conda songbird_env.
> 4. Modify DR_Jacob.py line 84 the fileplace. Make sure it points to the fileplace storing abundance data of all analyzed cases. 
> 5. Run DR_Jacob.py.
> 6. Deactivate songbird_env.
> 7. Modify QMD_Jacob.py line 7 the fileplace. Make sure it points to the fileplace storing abundance data of all analyzed cases. 
> 8. Run QMD_Jacob.py.
> 9. Set the workspace of R LSI_STOOL_analysis_code/plot.R to the fileplace storing abundance data of all analyzed cases.
> 10. Generate figure 2 and extended figure 1

To reproduce figure 1 a, b, c, d:
>1. Set the workspace of R LSI_relative_abundance_compare_with_density_data/plot_fig1abcd.R to the fileplace LSI_relative_abundance_compare_with_density_data.
>2. Run plot_fig1abcd.R.

To reproduce the analysis results of B1 and B2 validation, please follows the steps below:
> 1. Set the workspace of R to the fileplace storing abundance data of all analyzed cases, ANCOM_BC.R Line 427.
> 2. Run ANCOM_BC.R.
> 3. Go to conda songbird_env.
> 4. Modify DR_Doris.py line 84 the fileplace. Make sure it points to the fileplace storing abundance data of all analyzed cases. 
> 5. Run DR_Doris.py.
> 6. Deactivate songbird_env.
> 7. Modify QMD_Doris.py line 7 the fileplace. Make sure it points to the fileplace storing abundance data of all analyzed cases. 
> 8. Run QMD_Doris.py.
> 9. Set the workspace of R B1_B2_analysis_code/plot.R to the fileplace storing abundance data of all analyzed cases.
> 10. Generate figure 2 and extended figure 1


