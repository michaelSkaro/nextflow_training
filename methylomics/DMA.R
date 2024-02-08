#!/usr/bin/env Rscript
#####################################
# Date:February 7, 2024
# Author: Michael Skaro
# Purpose: Use the skeleton script to add methylation analsyis from illuminaEPIC array data to a nextflow pipeline
# 
# Usage: R script to be invoked and interacted with from the terminal.
# Parameters: 
# Outputs: 

# install the optparse package if not already installed
if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse", repos = "http://cran.us.r-project.org")
}

# install the libraries from the cran repo and the bioconductor repo to conduct the differential expression analysis in a nexflow pipeline

# install the tidyverse package if not already installed
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse", repos = "http://cran.us.r-project.org")
}

# install data.table package if not already installed
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table", repos = "http://cran.us.r-project.org")
}

# install the bioconductor and biocmanager packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
}



# install complexHeatmap package if not already installed
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

# install circlize package if not already installed
if (!requireNamespace("circlize", quietly = TRUE)) {
  BiocManager::install("circlize")
}

# install the grid package if not already installed
if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid", repos = "http://cran.us.r-project.org")
}

# install missMethyl package if not already installed
if (!requireNamespace("missMethyl", quietly = TRUE)) {
  BiocManager::install("missMethyl")
}

# install the minfi package if not already installed
if (!requireNamespace("minfi", quietly = TRUE)) {
  BiocManager::install("minfi")
}

# install the DMRcate package if not already installed
if (!requireNamespace("DMRcate", quietly = TRUE)) {
  BiocManager::install("DMRcate")
}

# install the illuminaHumanMethylationEPICanno.ilm10b2.hg19 package if not already installed
if (!requireNamespace("illuminaHumanMethylationEPICanno.ilm10b2.hg19", quietly = TRUE)) {
  BiocManager::install("illuminaHumanMethylationEPICanno.ilm10b2.hg19")
}

# install the minfiData package if not already installed
if (!requireNamespace("minfiData", quietly = TRUE)) {
  BiocManager::install("minfiData")
}

# install the illuminahumanmethlation450kanno.ilmn12.hg19 package if not already installed
if (!requireNamespace("illuminaHumanmethlation450kanno.ilmn12.hg19", quietly = TRUE)) {
  BiocManager::install("illuminaHumanmethlation450kanno.ilmn12.hg19")
}

# install the iiuminaHumanMethylation450kmanifest package if not already installed
if (!requireNamespace("iiuminaHumanMethylation450kmanifest", quietly = TRUE)) {
  BiocManager::install("iiuminaHumanMethylation450kmanifest")
}

# install the limma package if not already installed
if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}



# Call the libraries
library("optparse")
library("tidyverse")
library("data.table")
library("ComplexHeatmap")
library("circlize")
library("grid")
library("missMethyl")
library("minfi")
library("DMRcate")
library("illuminaHumanMethylationEPICanno.ilm10b2.hg19")
library("minfiData")
library("illuminaHumanmethlation450kanno.ilmn12.hg19")
library("iiuminaHumanMethylation450kmanifest")
library("limma")

# Export the loaded libraries in the session info into a formatted text file
sessionInfo()$otherPkgs %>% as.data.frame() %>% write.table("session_info.txt", sep="\t")

# create an opt list for the aruguements
option_list = list(
    make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-c", "--coldata"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-v", "--variable"), type="character", default="output/", 
              help="output file directory [default= %default]", metavar="character"),
    make_option(c("-a", "--annotation"), type="character", default="output/", 
              help="output file directory [default= %default]", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="output/", 
              help="output file directory [default= %default]", metavar="character"),
    make_option(c("-E", "--experiment"), type="character", default="experiment", 
              help="experiment name [default= %default]", metavar="character"))

# parse the arguments using thr arg parser. 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$out)) {
  stop("Please specify an output file directory with -o or --out");
}
# load the coldata from the argument --coldata
coldata <- data.table::fread(opt$coldata)


anno450k <- getAnnotation(illuminaHumanMethylation450kanno.ilmn12.hg19)

# ead in the epic annotation merged from the annotation arguement
annoEPIC <- getAnnotation(illuminaHumanMethylationEPICanno.ilm10b2.hg19)
# checkjokeroo for this file

# assign the idatWD from the argument --directory
idatWD <- opt$directory

# copy the coldata object to the output directory
write.csv(coldata, file=paste0(opt$out, "coldata.csv"))



prepare_targets <- function(idatWD){
    # list the files in the idatWD
    idatFiles <- list.files(idatWD, pattern="idat", full.names=TRUE)
    # check if the idatFiles are empty
    if (length(idatFiles) == 0) {
        stop("No idat files found in the directory")
    }
    
    # use the metharray package to read the idat files
    targets <- metharray::read.metharray.exp(targets = idatFiles)
    # return the targets
    return(targets)

}

prepare_RG_set <- function(targets){
    Arguements: targets
    returns: swan normalized RGset as an mSet 

    RGset <- read.metharray.exp(targets = targets, force = TRUE)
    # annotate the set 
    sampleNames(RGset) <- targets$Sample_Name
    ann_epic <- getAnnotation(illuminaHumanMethylationEPICanno.ilm10b2.hg19)
    annotation(RGset) <- ann_epic
    RGset@annotation <- c(array = "IlluminaMethylationEPICv2", annotation = "20a.hg38")
    # filter low quality probes
    detP <- detectionP(RGset)
    keep <- detP > 0.01
    RGset <- RGset[,keep]
    # Swan normalization
    msetSq <- preprocessSWAN(RGset)
    # save the objects 
    save(RGset, file = paste0(opt$out, "RGset.RData"))
    save(msetSq, file = paste0(opt$out, "msetSq.RData"))
    save(detP, file = paste0(opt$out, "detP.RData"))
    rm(keep, detP)

    return(msetSq)
  
    
}

# Prepare the betas 

prpepareBetas <- function(msetSq, RGset, targets){
    # Arguments: msetSq
    # Returns: betas as a data frame
    var <- data.frame(getBeta(msetSq)) 
    betas <- data.frame(probes = rownames(var), var)
    fwrite(betas, file=paste0(opt$out, "betas.txt", row.names=FALSE, sep="\t"))
    # sort the RGset

    RGset <- RGset[, sort[colnames(RGset)]]
    # sort the targets by the sample_Name
    targets <- targets[order(targets$Sample_Name),]
    # save the RGset and targets
    save(RGset, file=paste0(opt$out, "RGset.RData"))
    save(targets, file=paste0(opt$out, "targets.RData"))

    rm(var, RGset, targets)
    return(betas) 


}

# pepare the beta long value

prepareBetaLong <- function(betas, coldata){
    # Arguments: betas
    # Returns: betas long as a data frame
    betas_long <- betas %>% 
        pivot_longer(cols = -probes, names_to = "Sample_Name", values_to = "Beta") %>%
        left_join(coldata, by = "Sample_Name")
    fwrite(betas_long, file=paste0(opt$out, "betas_long.txt", row.names=FALSE, sep="\t"))
    return(betas_long)
}

# prepare the mVals

prepareMvals <- function(msetSq){
    # Arguments: msetSq
    # Returns: mvals as a data frame
    mvals <- data.frame(getM(msetSq))
    mvals <- data.frame(probes = rownames(mvals), mvals)
    fwrite(mvals, file=paste0(opt$out, "mvals.txt", row.names=FALSE, sep="\t"))
    
    return(mvals)
}

# preparemValslong

prepareMvalsLong <- function(mvals, coldata){
    # Arguments: mvals
    # Returns: mvals long as a data frame
    mvals_long <- mvals %>% 
        pivot_longer(cols = -probes, names_to = "Sample_Name", values_to = "mVal") %>%
        left_join(coldata, by = "Sample_Name")
    fwrite(mvals_long, file=paste0(opt$out, "mvals_long.txt", row.names=FALSE, sep="\t"))
    return(mvals_long)
}

# DMRCate

myAnnotation <- getAnnotation(IlluminaMethylationEPICv2)

# prepare the DMRcate

prepareDMRcate(coldata, mvals, annoEPIC, RGset, msetSq, targets){
    # Arguments: coldata, mvals, annoEPIC, RGset, msetSq, targets
    # Returns: DMRcate object
    # prepare the DMRcate object
    design <- model.matrix(~ 0 + Group, data = coldata)
    fit <- lmfit(mvals, design)
    contMat <- makeContrasts(Group1 - Group2, levels = design)
    fit2 <- contrasts.fit(fit, contMat)
    fit2 <- eBayes(fit2)
    dat <- as.data.frame(summary(decideTests(fit2)))
    dat.sig <- dat[dat$adj.P.Val < 0.05,]

    # create a DMPs object from the toptable fiunctions

    dmp <- toptable(fit2, coef = 1, number = Inf, adjust.method = "BH", p.value = 0.05, genelist = annoEPIC, lfc = 0.1) %>%
      left_join(annoEPIC, by = "probes") %>%
      filter(!is.na(chr)) %>%
      filter(chr != "X" & chr != "Y" & chr != "NA")

    mvals2 <- as.matrix(mvals)
    region_locs <- cpg.annotate(mvals2, datatype = "array", what = "M", analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contMat, 
    arraytype = "EPIC", annotation = annoEPIC, p.adjust.method = "BH", p.value = 0.05, lfc = 0.1, block.size = 100, cores = 1, verbose = TRUE)

    DMRs <- dmrcate(region_locs, lamda = 1000, C=2)
    # get the results ranges
    results <- getResults(DMRs)
    results.ranges extractRanges(DMR, genome = "hg38")

    # save the results
    save(dmp, file=paste0(opt$out, "dmp.RData"))
    save(DMRs, file=paste0(opt$out, "DMRs.RData"))
    save(results, file=paste0(opt$out, "results.RData"))
    save(results.ranges, file=paste0(opt$out, "results.ranges.RData"))
    rm(design, fit, contMat, fit2, dat, dat.sig, dmp, mvals2, region_locs, results, results.ranges)
    return(DMRs)
}



