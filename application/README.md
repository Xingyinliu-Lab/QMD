QMD application on human gut microbiome on common diseases compared with healthy population
=================
## Folder structure
```
application # QMD applications on human gut microbiome in common diseases compared with the healthy population
│  README.md # this file
│  
├─applicationScripts # place to store analysis and chart-plotting scripts
│  │  batchrun_VGMrepo.py # analysis entrance
│  │  plot_GMrepo.R # plot charts shown in ./Data/figures
│  │  stack_GMrepo_result.py # stack the analysis result
│  │  
│  └─QMD_GMREPO # a slightly modified version of QMD in the application analysis
│          plot_cost_VF.py # plot the optimization process 
│          qmd_dataPreprocessing_VF.py # single-thread data preprocessing
│          qmd_dataPreprocessing_VF_multi_thread.py # multi-threads data preprocessing
│          qmd_optimization_VF.py # conduct the QMD
│          
└─Data # QMD application analysis results
    ├─figures #generate figure 4.
    │      hypertension_diff.pdf #generate figure 4.b
    │      hypertension_foldchange.pdf #generate figure 4.b
    │      IBS_diff.pdf #generate figure 4.a
    │      IBS_foldchange.pdf #generate figure 4.a
    │      
    ├─gutMicrobiomeAbundanceData # input and output of QMD application
    │      alltaxaID.csv # ncbi taxa id list
    │      genus_AbundanceData.csv # genus level relative abundance data of all analyzed cases
    │      genus_AbundanceData.xlsx #xlsx version of genus level abundance data
    │      genus_AbundanceData_H2029_hypertension_all_taxa_detectionRate_info.csv # the taxa detection rate info in hypertension and healthy population QMD analysis.
    │      genus_AbundanceData_H2029_hypertension_analysis.csv # all analysis results in hypertension and healthy population QMD analysis.
    │      genus_AbundanceData_H2029_hypertension_analysis_stacked.csv # stacked analysis results in hypertension and healthy population QMD analysis.
    │      genus_AbundanceData_H2029_hypertension_cost.csv # the optimization process records
    │      genus_AbundanceData_H2029_hypertension_taxa_abundance_IntoModel.csv # relative abundance of remained taxa after removing low detection rate taxa
    │      genus_AbundanceData_H2029_hypertension_taxa_IntoModel.csv # Remained taxa after removing low detection rate taxa
    |      ... # data of other diseases
    │      genus_taxaid_list.csv # taxa name
    │      
    └─metadata # metadata for data pipelined into analysis
            autoimmuneD.csv # metadata for auto immune diseases
            CardiovascularDiseases.csv # metadata for cardiovascular diseases
            ...# metadata for other conditions
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
```
Python 3.7.7 was used to conduct the QMD analysis processes. Following packages are needed.
```
numpy
pandas
matplotlib
scipy
sklearn
ancomP
multiprocessing
```
Note that ancomP is the python implementation of ANCOM, it can be found in https://github.com/mortonjt/ancomP

## Reproduce the analysis results
To reproduce the analysis results, please follows the steps below:
> 1. Make sure all dependencies are installed. A conda environment are recommended.
> 2. Modify batchrun_VGMrepo.py line 4 the fileplace. Make sure it points to the fileplace storing genus level relative abundance data of all analyzed cases. The processers in Line 10 should also be modified to suit your computer configuration.
> 3. Run batchrun_VGMrepo.py
> 4. Modify stack_GMrepo_result.py line 6 the fileplace. Make sure it points to the fileplace storing genus level relative abundance data of all analyzed cases.
> 5. Run stack_GMrepo_result.py
> 6. Set the workspace of R to the fileplace storing genus level relative abundance data of all analyzed cases.
> 7. Run plot_GMrepo.R.

