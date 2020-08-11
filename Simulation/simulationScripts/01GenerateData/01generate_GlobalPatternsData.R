## the R code generate original abundance data from Global Patterns
## this code was copied from https://doi.org/10.1371/journal.pcbi.1003531.s001

# The required package list:
reqpkg = c("cluster", "doParallel", "edgeR", "DESeq", "DESeq2", "foreach", "ggplot2", 
    "grid", "scales", "metagenomeSeq", "phyloseq", "plyr", "reshape2", "ROCR")
# Load all required packages and show version
for (i in reqpkg) {
    print(i)
    print(packageVersion(i))
    library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}


data("GlobalPatterns")
# load(savedatafile)

# The delimiter in the command parameter string
comdelim = "_"
minobs=1L
# Define the different biological source templates to use sampletypes =
# c('Ocean', 'Soil')
# sampletypes = levels(get_variable(GlobalPatterns, "SampleType"))
sampletypes =c('Feces')
Ncores=7
# Define the ceiling in the number of OTUs to consider in the template
# multinomial nOTUs = 353L
nOTUs = 500L

# Define the number of samples in each class of a simulated experiment J =
# c(3, 10)
J = c(200)

# The different values of effect size to apply foldeffect = c(1.5, 5, 20)
foldeffect = c(0)

## The number of reads per sample, ns ns = c(2000, 5E4)
ns = c(1E6)

# Vector of the replicate numbers to repeat for each comb of simulation
# parameters (n, etc) reps=1:2
reps =1

# Define the simulation parameters combinations
simparams = apply(expand.grid(ns, sampletypes, reps, foldeffect, J), 1, paste0, 
    collapse = comdelim)
# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simparamslabels = c("nreads", "SampleType", "Replicate", "EffectSize", "nsamples")

# Define date-stamp for file names
datestamp = gsub(":", "_", gsub("[[:space:]]+", "-", date()), fixed = TRUE)
print(datestamp)


sampsums = sample_sums(GlobalPatterns)
keepsamples = sample_data(GlobalPatterns)$SampleType %in% sampletypes

template = prune_samples(keepsamples, GlobalPatterns)
# Make a list of source templates
templatelist = lapply(sampletypes, function(i, tempall, minobs, nOTUs) {
    cat(i, "\n")
    whtemp = (get_variable(tempall, "SampleType") %in% i)
    templatei = prune_samples(whtemp, tempall)
    samobs = apply(otu_table(templatei), 1, function(x, m) sum(x > m), m = minobs)
    otudf = data.frame(prev = samobs, sums = taxa_sums(templatei))
    otudf = otudf[order(-otudf$prev, -otudf$sums), ]
    # Trim all but the first nOTUs
    return(prune_taxa(rownames(otudf)[1:nOTUs], templatei))
}, template, minobs, nOTUs)


####
####Simulation Function
####Input
####postfix, template, J, n (number or reads)
####Output
####phyloseq object, incorporating simulation results and inputs
###

microbesim = function(postfix = "sim", template, J, n = 10000) {
    # Generate `J` simulated microbiomes with `n` total reads each (all the
    # same, or n has length equal to the value of `J`), with subsamples drawn
    # from `template`.  `postfix` is a dummy idenitifer added to help
    # distinguish simulated samples in downstream code.
    require("phyloseq")
    # call the proporitions vector `pi`, similar to nomenclature from DMN
    pi = taxa_sums(template)
    # n must be a scalar (recycled as the number of reads for every simulation)
    # or it can be vector of length equal to J, the number of samples being
    # simulated.
    if (length(J) != 1) {
        stop("Length of J should be 1.")
    }
    if (length(n) != 1 & length(n) != J) {
        stop("n should be length 1, or length J.")
    }
    # Actually create the simulated abundance table
    simat = mapply(function(i, x, sample.size) {
        if (FALSE) 
            {
                print(i)
            }  # i is a dummy iterator
        phyloseq:::rarefaction_subsample(x, sample.size)
    }, i = 1:J, sample.size = n, MoreArgs = list(x = pi), SIMPLIFY = TRUE)
    simat = t(simat)
    # Add the OTU names to the OTU (column) indices
    colnames(simat) <- names(pi)
    # Add new simulated sample_names to the row (sample) indices
    rownames(simat) <- paste(i, "::", 1:nrow(simat), postfix, sep = "")
    # Put simulated abundances together with metadata as a phyloseq object
    OTU = otu_table(simat, taxa_are_rows = FALSE)
    # Define data.frame that will become sample_data
    SDF = data.frame(sample = sample_names(OTU), TableNumber = i, type = "simulated")
    SDF$postfix <- postfix
    rownames(SDF) <- sample_names(OTU)
    SD = sample_data(SDF)
    # Return a phyloseq object
    return(phyloseq(OTU, SD))
}


# rescale the sum of reads in the original raw(-ish) template data to the
# expected library size being requested here
sumsim = function(n, sumtemplate, J) {
    # `n` - expected size target `sumtemplate` - the template vector of library
    # sizes observed in template `J` - The number of sample sizes to return
    scaledSums = round(n * (sumtemplate/median(sumtemplate)))
    return(sample(scaledSums, size = J, replace = TRUE))
}
########################
########################
########################
########################
########################
###########single thread

i = simparams[1]
params = strsplit(i, comdelim)[[1]]
names(params) <- simparamslabels
# Initialize
n = sim = sim1 = sim2 = n1 = n2 = NULL
# cat(i, '\n')
n = as.numeric(params["nreads"])
sampletypei = params["SampleType"]
# The number of samples to use for each class in this simulation
Ji = as.integer(params["nsamples"])
templatei = templatelist[[1]]


for (i in 1:5){
    tryAgain = TRUE
    infiniteloopcounter = 1
    while (tryAgain & infiniteloopcounter < 5) {
        n1 = sumsim(n, sampsums, Ji)
        sim1 = microbesim(paste0(sampletypei, ";grp1"), templatei, Ji, n1)
        if (is.null(sim1) | is.null(n1)  | inherits(sim1,  "try-error") ) {
            tryAgain = TRUE
            infiniteloopcounter = infiniteloopcounter + 1
        } else {
            tryAgain = FALSE
        }
    }
    if (infiniteloopcounter >= 5) {
        stop("Consistent error found during simulation. Need to investigate cause.")
    }
    # Merge the two simulated datasets together into one phyloseq object and add
    # back tree.
    sim = merge_phyloseq(sim1, tax_table(GlobalPatterns), phy_tree(GlobalPatterns))
    otu_data<-sim@otu_table@.Data
    otu_sum<-apply(otu_data, 1, sum)
    otu_sum<-matrix(rep(otu_sum, nOTUs),ncol=nOTUs)
    relative_otu<-otu_data/otu_sum
    filename= paste('gp_',i,'_relative_otu.csv',sep='')
    write.csv(relative_otu,filename)
}
# Rarely a simulation has a weird value and fails.  Catch these with `try`,
# and repeat the simulation call if error (it will be a new seed)
########################
########################
########################
########################
########################
###########multi thread

cl <- makeCluster(Ncores)
registerDoParallel(cl)

# generate 500 relative abundance data table from Global patterns
foreach(oo = 1:500, .packages = c("phyloseq")) %dopar% {
        
    i = simparams[1]
    params = strsplit(i, comdelim)[[1]]
    names(params) <- simparamslabels
    # Initialize
    n = sim = sim1 = sim2 = n1 = n2 = NULL
    # cat(i, '\n')
    n = as.numeric(params["nreads"])
    sampletypei = params["SampleType"]
    # The number of samples to use for each class in this simulation
    Ji = as.integer(params["nsamples"])
    templatei = templatelist[[1]]

    tryAgain = TRUE
    infiniteloopcounter = 1
    while (tryAgain & infiniteloopcounter < 5) {
        n1 = sumsim(n, sampsums, Ji)
        sim1 = microbesim(paste0(sampletypei, ";grp1"), templatei, Ji, n1)
        if (is.null(sim1) | is.null(n1)  | inherits(sim1,  "try-error") ) {
            tryAgain = TRUE
            infiniteloopcounter = infiniteloopcounter + 1
        } else {
            tryAgain = FALSE
        }
    }
    if (infiniteloopcounter >= 5) {
        stop("Consistent error found during simulation. Need to investigate cause.")
    }
    # Merge the two simulated datasets together into one phyloseq object and add
    # back tree.
    sim = merge_phyloseq(sim1, tax_table(GlobalPatterns), phy_tree(GlobalPatterns))
    otu_data<-sim@otu_table@.Data
    otu_sum<-apply(otu_data, 1, sum)
    otu_sum<-matrix(rep(otu_sum, nOTUs),ncol=nOTUs)
    relative_otu<-otu_data/otu_sum
    filename= paste('gp_',oo,'_relative_otu.csv',sep='')
    write.csv(relative_otu,filename)

}