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
│          
├─GUI_QMD # QMD@XYL_Lab with Graphical User Interface
│  
├─Validation # LSI STOOL B1 B2 Validation code and results
│                  
├─Simulation # benchmark simulation scripts and results
│              
└─application # QMD applications on human gut microbiome in common diseases compared with the healthy population
```
## <font color=#8A2BE2F>Getting started</font>
Both character user interface and graphical user interface version QMD@XYL_Lab are provided in this repo. Users familiar with python lanugage can start with the CUI version. This version provides more flexible functions. For MS windows users, one can also try the graphical user interface version. 
If one wants to reproduce the simulations, validations or applications results, one should install the corresponding dependencies and modify the file place set in the scripts. 

## <font color=#8A2BE2F>Citing QMD@XYL_Lab</font>
If you use QMD@XYL_Lab for any published research, please include the following citation:
```
QMD: Quantification of microbial density changes and its application on microbiome differential abundance analysis
```
