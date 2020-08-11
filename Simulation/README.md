QMD benchmark simulation 
=================
## Folder structure

```
Simulation # benchmark simulation scripts and results
 ├─result_Aanlysis # simulation result statistics comparision
 │  ├─collectedStatistics # H2029, GP, Obesity simulation statistics
 │  │      GP_stat.csv # GP simulation statistics
 │  │      H2029_stat.csv # H2029 simulation statistics
 │  │      Obesity_stat.csv # Obesity simulation statistics
 │  │      
 │  └─figures # generate figure 3, extended figure 3, extended figure 4, extended figure 6 from collectedStatistics.
 │          GP_changfold_changedtaxarate.pdf # extended figure 6 l
 │          GP_FNR.pdf # extended figure 4 d
 │          GP_FNR_on_ChangedTaxaRate.pdf # extended figure 4 e
 │          GP_FNR_on_groupsize.pdf # extended figure 4 f
 │          GP_FPR.pdf # extended figure 4 g
 │          GP_FPR_on_ChangedTaxaRate.pdf # extended figure 4 h
 │          GP_FPR_on_groupsize.pdf # extended figure 4 i
 │          GP_groupsize_changedtaxarate.pdf #  # extended figure 6 k
 │          GP_groupsize_changfold.pdf # extended figure 6 j
 │          GP_MAE_on_ChangedTaxaRate.pdf # extended figure 4 b
 │          GP_MAE_on_groupsize.pdf # extended figure 4 c
 │          GP_microbial_loads_foldchange_prediction_performance.pdf # extended figure 4 a
 │          GP_simu_con.pdf # extended figure 6 i
 │          H2029_changfold_changedtaxarate.pdf # extended figure 6 d
 │          H2029_FNR.pdf # figure 3 d
 │          H2029_FNR_on_ChangedTaxaRate.pdf # figure 3 e
 │          H2029_FNR_on_groupsize.pdf # figure 3 f
 │          H2029_FPR.pdf # figure 3 g
 │          H2029_FPR_on_ChangedTaxaRate.pdf  # figure 3 h
 │          H2029_FPR_on_groupsize.pdf # figure 3 i
 │          H2029_groupsize_changedtaxarate.pdf # extended figure 6 c
 │          H2029_groupsize_changfold.pdf # extended figure 6 b
 │          H2029_MAE_on_ChangedTaxaRate.pdf  # figure 3 b
 │          H2029_MAE_on_groupsize.pdf # figure 3 c
 │          H2029_microbial_loads_foldchange_prediction_performance.pdf # figure 3 a
 │          H2029_simu_con.pdf # extended figure 6 a
 │          Obesity_changfold_changedtaxarate.pdf  # extended figure 6 h
 │          Obesity_FNR.pdf # extended figure 3 d
 │          Obesity_FNR_on_ChangedTaxaRate.pdf # extended figure 3 e
 │          Obesity_FNR_on_groupsize.pdf # extended figure 3 f
 │          Obesity_FPR.pdf # extended figure 3 g
 │          Obesity_FPR_on_ChangedTaxaRate.pdf # extended figure 3 h
 │          Obesity_FPR_on_groupsize.pdf # extended figure 3 i
 │          Obesity_groupsize_changedtaxarate.pdf # extended figure 6 g
 │          Obesity_groupsize_changfold.pdf # extended figure 6 f
 │          Obesity_MAE_on_ChangedTaxaRate.pdf # extended figure 3 b
 │          Obesity_MAE_on_groupsize.pdf # extended figure 3 c
 │          Obesity_microbial_loads_foldchange_prediction_performance.pdf # extended figure 3 a
 │          Obesity_simu_con.pdf # extended figure 6 e
 │          
 ├─simulationData # generated simulation input and output data
 │  ├─GPoriData # input data(relative abundance, absolute abundance, simulation change record) and output data for GP series
 │  │      gp_100_absolute_abundance_simed.csv # absolute abundance table of the 100th instance
 │  │      gp_100_analysis_all_methods.csv # final simulation results including ANCOM, ANCOM-BC, DR, QMD, QDM_qvalue, MWU, RAC of the 100th instance
 │  │      gp_100_analysis_ANCOM_ANCOM-BC_QMD.csv #simulation results including ANCOM, ANCOM-BC, QMD, QDM_qvalue of the 100th instance
 │  │      gp_100_analysis_ANCOM_ANCOM-BC_QMD_stacked.csv #stacked simulation results including ANCOM, ANCOM-BC, QMD, QDM_qvalue of the 100th instance
 │  │      gp_100_analysis_stacked_all_methods.csv # stacked final simulation results including ANCOM, ANCOM-BC, DR, QMD, QDM_qvalue, MWU, RAC of the 100th instance
 │  │      gp_100_ANCOM_BC_res1.csv # ANCOM-BC estimated bias of the 100th instance
 │  │      gp_100_ANCOM_BC_res2.csv # ANCOM-BC microbial density estimation and differentially abundant taxa (DA) identification result of the 100th instance
 │  │      gp_100_changed_record_simed.csv # simulation change record (the true change in data generation processes) of the 100th instance
 │  │      gp_100_Control_Treat_all_taxa_detectionRate_info.csv # taxa detection rate info for the 100th instance of the 100th instance
 │  │      gp_100_Control_Treat_taxa_abundance_IntoModel.csv # relative abundance of remained taxa after removing low detection rate taxa of the 100th instance
 │  │      gp_100_Control_Treat_taxa_Into_Model.csv # Remained taxa after removing low detection rate taxa of the 100th instance
 │  │      gp_100_relative_otu_simed.csv # relative abundance table of the 100th instance
 │  │      ... # other GP simulation instances data (totally 500)
 │  │      
 │  ├─H2029 # input data(relative abundance, absolute abundance, simulation change record) and output data for H2029 series
 │  │      H2029_100_absolute_abundance_simed.csv # absolute abundance table of the 100th instance
 │  │      H2029_100_analysis_all_methods.csv # final simulation results including ANCOM, ANCOM-BC, DR, QMD, QDM_qvalue, MWU, RAC of the 100th instance
 │  │      H2029_100_analysis_ANCOM_ANCOM-BC_QMD.csv #simulation results including ANCOM, ANCOM-BC, QMD, QDM_qvalue of the 100th instance
 │  │      H2029_100_analysis_ANCOM_ANCOM-BC_QMD_stacked.csv #stacked simulation results including ANCOM, ANCOM-BC, QMD, QDM_qvalue of the 100th instance
 │  │      H2029_100_analysis_stacked_all_methods.csv # stacked final simulation results including ANCOM, ANCOM-BC, DR, QMD, QDM_qvalue, MWU, RAC of the 100th instance
 │  │      H2029_100_ANCOM_BC_res1.csv # ANCOM-BC estimated bias of the 100th instance
 │  │      H2029_100_ANCOM_BC_res2.csv # ANCOM-BC microbial density estimation and differentially abundant taxa (DA) identification result of the 100th instance
 │  │      H2029_100_changed_record_simed.csv # simulation change record (the true change in data generation processes) of the 100th instance
 │  │      H2029_100_Control_Treat_all_taxa_detectionRate_info.csv # taxa detection rate info for the 100th instance of the 100th instance
 │  │      H2029_100_Control_Treat_taxa_abundance_IntoModel.csv # relative abundance of remained taxa after removing low detection rate taxa of the 100th instance
 │  │      H2029_100_Control_Treat_taxa_Into_Model.csv # Remained taxa after removing low detection rate taxa of the 100th instance
 │  │      H2029_100_relative_otu_simed.csv # relative abundance table of the 100th instance
 │  │      ... # other H2029 simulation instances data (totally 500)
 │  │      
 │  ├─Obesity # input data(relative abundance, absolute abundance, simulation change record) and output data for Obesity series
 │  │      Obesity_100_absolute_abundance_simed.csv # absolute abundance table of the 100th instance
 │  │      Obesity_100_analysis_all_methods.csv # final simulation results including ANCOM, ANCOM-BC, DR, QMD, QDM_qvalue, MWU, RAC of the 100th instance
 │  │      Obesity_100_analysis_ANCOM_ANCOM-BC_QMD.csv #simulation results including ANCOM, ANCOM-BC, QMD, QDM_qvalue of the 100th instance
 │  │      Obesity_100_analysis_ANCOM_ANCOM-BC_QMD_stacked.csv #stacked simulation results including ANCOM, ANCOM-BC, QMD, QDM_qvalue of the 100th instance
 │  │      Obesity_100_analysis_stacked_all_methods.csv # stacked final simulation results including ANCOM, ANCOM-BC, DR, QMD, QDM_qvalue, MWU, RAC of the 100th instance
 │  │      Obesity_100_ANCOM_BC_res1.csv # ANCOM-BC estimated bias of the 100th instance
 │  │      Obesity_100_ANCOM_BC_res2.csv # ANCOM-BC microbial density estimation and differentially abundant taxa (DA) identification result of the 100th instance
 │  │      Obesity_100_changed_record_simed.csv # simulation change record (the true change in data generation processes) of the 100th instance
 │  │      Obesity_100_Control_Treat_all_taxa_detectionRate_info.csv # taxa detection rate info for the 100th instance of the 100th instance
 │  │      Obesity_100_Control_Treat_taxa_abundance_IntoModel.csv # relative abundance of remained taxa after removing low detection rate taxa of the 100th instance
 │  │      Obesity_100_Control_Treat_taxa_Into_Model.csv # Remained taxa after removing low detection rate taxa of the 100th instance
 │  │      Obesity_100_relative_otu_simed.csv # relative abundance table of the 100th instance
 │  │      ... # other Obesity simulation instances data (totally 500)
 │  │      
 │  └─originalDataPool # data pool to generate relative abundance and absolute abundance data
 │      ├─GP # data pool to generate relative abundance and absolute abundance data for GP series
 │      │      gp_100_relative_otu.csv # the 100th GP instance data pool
 │      │      gp_101_relative_otu.csv # the 101th GP instance data pool
 │      │      ... # other GP instance data pools (totally 500)
 │      │      
 │      ├─H2029 # data pool to generate relative abundance and absolute abundance data for H2029 series
 │      │      genus_taxaid_list.csv # taxa name
 │      │      H2029_metadata.csv # H2029 metadata, including the related NCBI project etc.
 │      │      H2029_relative_abundance.csv # H2029 relative abundance including runid
 │      │      H2029_relative_otu.csv # H2029 relative abundance
 │      │      
 │      └─Obesity # data pool to generate relative abundance and absolute abundance data for Obesity series
 │              genus_taxaid_list.csv # taxa name
 │              Obesity_metadata.csv # Obesity metadata, including the related NCBI project etc.
 │              Obesity_relative_abundance.csv # Obesity relative abundance including runid
 │              Obesity_relative_otu.csv # Obesity relative abundance
 │              
 └─simulationScripts # place to store the scripts in simulation
     ├─01GenerateData # scripts generate all input data
     │      01generate_GlobalPatternsData.R # Following McMurdie, Paul J's work generating GP data pool
     │      02RandomGlobalPatternsData.py # From GP data pool, generate simualted relative abundance and absolute abundance data for GP series
     │      03RandomH2029Data.py # From H2029 data pool, generate simualted relative abundance and absolute abundance data for H2029 series
     │      04RandomObesityData.py # From Obesity data pool, generate simualted relative abundance and absolute abundance data for Obesity series
     │      
     ├─02ANCOM_ANCOM-BC_QMD_DR # conducting simualtion 
     │      ANCOM_BC.R # ANCOM-BC analysis scripts
     │      batchrun_DR.py # DR analysis scripts entrance
     │      batchrun_GP.py # QMD analysis entrance for GP series
     │      batchrun_H2029.py # QMD analysis entrance for H2029 series
     │      batchrun_Obesity.py # QMD analysis entrance for Obesity series
     │      DR.py  # DR analysis code
     │      qmd_ANCOM.py # QMD and ANCOM analysis code
     │      qmd_dataPreprocessing.py # QMD data preprocessing code
     │      
     └─03Simulation_result_analysis # analysis result and plot charts
             GlobalPatterns_stat.py # collect statistics from GP data series simulation analysis result
             H2029_stat.py # collect statistics from H2029 data series simulation analysis result
             Obesity_stat.py # collect statistics from Obesity data series simulation analysis result
             plot_and_analysis_onSimulationResult.R # plot statistics chart
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
library(cluster)
library(doParallel)
library(edgeR)
library(DESeq)
library(DESeq2)
library(foreach)
library(grid)
library(scales)
library(metagenomeSeq)
library(phyloseq)
library(plyr)
library(reshape2)
library(ROCR)
```
Note 1: R code implementation of [ANCOM-BC](https://doi.org/10.1038/s41467-020-17041-7). The code is a copy from [the author's github](https://github.com/FrederickHuangLin/ANCOM-BC/blob/master/scripts/ancom_bc_v1.0.R)

Note 2: R code implementation of GP data pool generation was copied from Paul J. McMurdie's work [Waste Not, Want Not: Why Rarefying Microbiome Data Is Inadmissible](https://doi.org/10.1371/journal.pcbi.1003531.s001).

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
Note 3: ancomP is the python implementation of ANCOM, it can be found in https://github.com/mortonjt/ancomP
Note 4: Python code implementation for [DR](https://doi.org/10.1038/s41467-019-10656-5). The author did not provide a standalone implementation of DR. We extracted the implementation code from [the simulatin pynote from the author](https://github.com/knightlab-analyses/reference-frames/blob/master/ipynb/simulation-benchmark.ipynb). The author recommended to install following packages and setup a conda virtual env to run the DR code
```
conda create -n songbird_env python=3.6 numpy=1.15.4 scikit-bio=0.5.5 seaborn pandas=0.23.4 -c conda-forge
source activate songbird_env
conda install tensorflow=1.10 tqdm nomkl
conda install biom-format h5py -c conda-forge
conda install jupyter notebook
conda install songbird -c conda-forge
```


## Reproduce the analysis results

To regenerate data pool for simulation, please follows the steps below:
Note 5: the generated data pool might be a litte different from the repo. This is caused by the random seed set.
Note 6: there's a little difference in generating GP simualtion data and generating H2029 and Obesity data. The GP series following McMurdie's work and we generated 500 data instances as the data pool. Each simulation dataset was resampled and generated from one data instance. For the H2029 and Obesity data, there's only one dataset as the data pool for each simulation series and  each simulation dataset was resampled from the same data pool.
> 1. Run 01generate_GlobalPatternsData.R, this will output 500 GP data pool in originalDataPool/GP
> 2. Modify 02RandomGlobalPatternsData.py line 8 to the fileplace storing GP data pool
> 3. Run 02RandomGlobalPatternsData.py, this will output 500 GP simulated absolute abundance table, relative abundance table and the changed records in simulationData/GPoriData
> 4. Modify 03RandomH2029Data.py line 8 to the fileplace storing H2029 data pool
> 5. Run 03RandomH2029Data.py, this will output 500 H2029 simulated absolute abundance table, relative abundance table and the changed records in simulationData/H2029
> 6. Modify 04RandomObesityData.py line 8 to the fileplace storing Obesity data pool
> 7. Run 04RandomObesityData.py, this will output 500 Obesity simulated absolute abundance table, relative abundance table and the changed records in simulationData/Obesity



To reproduce the analysis results of simulation, please follows the steps below:
> 1. Set the workspace of R to the fileplace storing abundance data of all analyzed cases, ANCOM_BC.R Line 427,469,511.
> 2. Run ANCOM_BC.R.
> 3. Modify batch_GP.py line 20 the fileplace to the fileplace storing abundance data of GP series.  
> 4. Run batch_GP.py.
> 5. Modify batch_Obesity.py line 20 the fileplace to the fileplace storing abundance data of Obesity series. 
> 6. Run batch_Obesity.py.
> 7. Modify batch_H2029.py line 20 the fileplace to the fileplace storing abundance data of H2029 series. 
> 8. Run batch_H2029.py.
> 9. Go to conda songbird_env.
> 10. Modify batchrun_DR.py line 22, 33, 45 the fileplace to the fileplace storing abundance data of H2029, Obesity, GP series respectively. 
> 11. Run batchrun_DR.py.
> 12. Deactivate songbird_env.


To reproduce simulation result analysis and chart plotting:
> 1. Modify GlobalPatterns_stat.py line 5 the fileplace to the fileplace storing simulation results of GP series.  
> 2. Run GlobalPatterns_stat.py.
> 3. Modify H2029_stat.py line 5 the fileplace to the fileplace storing simulation results of H2029 series.  
> 4. Run H2029_stat.py.
> 5. Modify Obesity_stat.py line 5 the fileplace to the fileplace storing simulation results of Obesity series.  
> 6. Run Obesity_stat.py.
> 7. Modify plot_and_analysis_onSimulationResult.R line 21, set the R workspace to the fileplace storing simulation results statistics.  
> 8. Run plot_and_analysis_onSimulationResult.R.

