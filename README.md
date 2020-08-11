# <center><font color=#8A2BE2F>QMD@XYL_Lab</font></center>

Development of DNA sequencing techniques have enabled our insights into composition and dynamics of complex microbial communities from human gut to soils and oceans. In practice, the microbiome data can be divided into four categories.

> 1. The *Microbial Load* refers to the number of a taxa in a sample.
> 2. The *Microbial Density* refers to the microbial loads per unit sample mass (e.g. volume or weight).
> 3. The *Absolute Abundance* refers to the microbial feature raw counts of samples generated from the sequencing platforms.
> 4. The *Relative Abundance* are compositional data and summed to a constant. 

Microbial density is the most unbiased data presenting the real microbiome world. Quantification of Microbial Density Changes (QMD) is aiming for estimation of the microbial density changes between two experimental groups. Based on the result from QMD, statistical tests for differentially abundant taxa (QMDD) are constructed and p value are provided to nail the statistically significant differentially abundant taxa.

## <font color=#8A2BE2F>Assumption of QMD and QMDD</font>

> The microbial density changes of a considerable part of taxa are relatively small.

## <font color=#8A2BE2F>QMD@XYL_Lab</font>
QMD@XYL_Lab helps researchers to conduct QMD and QMDD in differentially abundant taxa identification and different abundance estimations between experimental groups. This software was implemented by python3. QMD@XYL_Lab is open source software. One can find source code of QMD@XYL_Lab at [github](https://github.com/Xingyinliu-Lab/QMD). Both GUI and CUI version of the software are provided. 

## <font color=#8A2BE2F>Folder Structure of this repository</font>
```
QMD # root of the repository
│
├─CUI_QMD # Character User Interface QMD@XYL_Lab
│  │  
│  └─demoProject # a demo project in CUI_QMD
│          
├─GUI_QMD # QMD@XYL_Lab with Graphical User Interface
│  │  
│  ├─release # zipped release version of the GUI_QMD software
│  │      
│  ├─release_v1.0.20200810 # v1.0.20200810 GUI_QMD
│  │  │  
│  │  ├─demoData # demo data for GUI_QMD
│  │  │      
│  │  ├─Project # place to store the project analysis result
│  │  │  └─demo_project_b2 # place to store the demoproject analysis result
│  │  │      │  
│  │  │      └─datapreprocessing # place to store the demoproject preprocessing result
│  │  │              
│  │  └─source # source data file containing manual.html etc. to support GUI_QMD
│  │      │  
│  │      └─fonts # source fonts data
│  │              
│  └─sourcecode # source code for GUI_QMD
│      │  
│      └─source # source data file containing manual.html etc. to support GUI_QMD
│          │  
│          └─fonts # source fonts data
│  
├─Validation # LSI STOOL B1 B2 Validation code and results
│   ├─BarlowWork_LSI_STOOL # LSI STOOL validation
│   │  ├─LSI_STOOL_analysis_code # validation scripts
│   │  │      
│   │  ├─Data # input and output data in LSI STOOL validation
│   │  │      
│   │  ├─figures # Generate figure 2 and extended figure 1.
│   │  │      
│   │  └─LSI_relative_abundance_compare_with_density_data # Generate figure 1.a
│   │          
│   └─VandeputteWork_B1_B2 # B1 B2 validation
│       ├─B1_B2_analysis_code # validation scripts
│       │      
│       ├─Data # input and output data in B1 B2 validation
│       │      
│       └─figures # Generate extended figure 2.
│
│                  
├─Simulation # benchmark simulation scripts and results
│  ├─result_Aanlysis # simulation result statistics comparision
│  │  ├─collectedStatistics # H2029, GP, Obesity simulation statistics
│  │  │      
│  │  └─figures # generate figure 3, extended figure 3, extended figure 4, extended figure 6.
│  │          
│  ├─simulationData # generated simulation input and output data
│  │  ├─GPoriData # input data(relative abundance, absolute abundance, simulation change record) and output data for GP series
│  │  │      
│  │  ├─H2029 # input data(relative abundance, absolute abundance, simulation change record) and output data for H2029 series
│  │  │      
│  │  ├─Obesity # input data(relative abundance, absolute abundance, simulation change record) and output data for Obesity series
│  │  │      
│  │  └─originalDataPool # data pool to generate relative abundance and absolute abundance data
│  │      ├─GP # data pool to generate relative abundance and absolute abundance data for GP series
│  │      │      
│  │      ├─H2029 # data pool to generate relative abundance and absolute abundance data for H2029 series
│  │      │      
│  │      └─Obesity # data pool to generate relative abundance and absolute abundance data for Obesity series
│  │              
│  └─simulationScripts # place to store the scripts in simulation
│      ├─01GenerateData # scripts generate all input data
│      │      
│      ├─02ANCOM_ANCOM-BC_QMD_DR # conducting simualtion 
│      │      
│      └─03Simulation_result_analysis # analysis result and plot charts
│              
└─application # QMD applications on human gut microbiome in common diseases compared with the healthy population
   │  
   ├─applicationScripts # place to store analysis and chart-plotting scripts
   │  │  
   │  └─QMD_GMREPO # a slightly modified version of QMD in the application analysis
   │          
   └─Data # QMD application analysis results
       ├─figures #generate figure 4.
       │      
       ├─gutMicrobiomeAbundanceData # input and output of QMD application
       │      
       └─metadata # metadata for data pipelined into analysis
```
## <font color=#8A2BE2F>Getting started</font>
Both character user interface and graphical user interface version QMD@XYL_Lab are provided in this repo. Users familiar with python lanugage can start with the CUI version. This version provides more flexible functions. For MS windows users, one can also try the graphical user interface version. 
If one wants to reproduce the simulations, validations or applications results, one should install the corresponding dependencies and modify the file place set in the scripts. 

## <font color=#8A2BE2F>Citing QMD@XYL_Lab</font>
If you use QMD@XYL_Lab for any published research, please include the following citation:
```
QMD: Quantification of microbial density changes and its application on microbiome differential abundance analysis
```
