#!/usr/bin/env Rscript
#####################################
# Date:February 7, 2024
# Author: Michael Skaro
# Purpose: Use the skeleton script to add the differential expression analysis functions
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

# install the DESeq2 package using the biocmanager if not already installed

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}

# install complexHeatmap package if not already installed
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}

# install tximport package if not already installed
if (!requireNamespace("tximport", quietly = TRUE)) {
  BiocManager::install("tximport")
}

# install circlize package if not already installed
if (!requireNamespace("circlize", quietly = TRUE)) {
  BiocManager::install("circlize")
}

# install the grid package if not already installed
if (!requireNamespace("grid", quietly = TRUE)) {
  install.packages("grid", repos = "http://cran.us.r-project.org")
}

# Call the libraries
library("optparse")
library("tidyverse")
library("data.table")
library("DESeq2")
library("ComplexHeatmap")
library("tximport")
library("circlize")
library("grid")

# Export the loaded libraries in the session info into a formatted text file
sessionInfo()$otherPkgs %>% as.data.frame() %>% write.table("session_info.txt", sep="\t")

# create an opt list for the aruguements
option_list = list(
    make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-c", "--coldata"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-g", "--genes"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-v", "--variable"), type="character", default="output/", 
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
# load the genes files from the argument --genes
genes <- data.table::fread(opt$genes)
# create a files object that uses the directory from the arguements, the SampleDir column in the coldata files, and quant.sf as a files list
files <- file.path(opt$directory, coldata$SampleDir, "quant.sf")
# create a tx2gene object that uses the genes file and the gene_id and transcript_id columns
tx2gene <- genes[, c("gene_id", "transcript_id")]
# create a txi object that uses the tximport function and the files and type arguments
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# use the value from the variable argument to design the dds object.
# if the variable == "Compound", then use the Compound column in the coldata file, if the variable == "Time", then use the Time column in the coldata file, if the variable == "Dose", 
# then use the Dose column in the coldata file, if the variable == "Treatment", if the variable == "Disease", then use the Disease column in the coldata file.
if(opt$variable == "Compound") {
  dds <- DESeqDataSetFromTximport(txi, coldata, ~ Compound)
} else if(opt$variable == "Time") {
  dds <- DESeqDataSetFromTximport(txi, coldata, ~ Time)
} else if(opt$variable == "Dose") {
  dds <- DESeqDataSetFromTximport(txi, coldata, ~ Dose)
} else if(opt$variable == "Treatment") {
  dds <- DESeqDataSetFromTximport(txi, coldata, ~ Treatment)
} else if(opt$variable == "Disease") {
  dds <- DESeqDataSetFromTximport(txi, coldata, ~ Disease)
} else {
  stop("Please specify a variable with -v or --variable")
}

# create a filter function for the dds object for low TPM values
tpm_cutoff <- 1
filter.tpm <- function(txi, tpm_cutoff,coldata, variable) {
  cutoff <- tpm_cutoff
  tpm.table <-data.frame(txi$abundance)
  names(tpm.table) <- coldata$SampleID
  tpm.table.genes <- tpm.table %>%
    rownames_to_column(var = gene.names) %>%
    as.data.frame()
  # the number of conidtions
  if(variable == "Compound") {
    ncond <- as.vector(unique(coldata$Compound))
  } else if(variable == "Time") {
    ncond <- as.vector(unique(coldata$Time))
  } else if(variable == "Dose") {
    ncond <- as.vector(unique(coldata$Dose))
  } else if(variable == "Treatment") {
    ncond <- as.vector(unique(coldata$Treatment))
  } else if(variable == "Disease") {
    ncond <- as.vector(unique(coldata$Disease))
  } else {
    stop("Please specify a variable with -v or --variable")
  }
  tpm.table.long <- tpm.table.genes %>%
    pivot_longer(cols = -gene.names, names_to = "SampleID", values_to = "TPM") %>%
    mutate(condition = if_else(grepl(ncond[1], SampleID), ncond[1], ncond[2])) %>%
    group_by(gene.names, conditions) %>%
    dplyr::summarize(median _TPM = median(TPM), .groups = "drop") %>%
    group_by(gene.names) %>%
    dplyr::summarize(median_TPM = as.numeric(median(median_TPM))) %>%
    filter(median_TPM >= cutoff)
  filter.names <- unique(tpm.table.long$gene.names)
  txi.gene.list <- rownames(txi$abundance)
  # filter the txi object
  txi$abundance <- txi$abundance[which(txi.gene.list %in% filter.names),]
  txi$counts <- txi$counts[which(txi.gene.list %in% filter.names),]
  txi$length <- txi$length[which(txi.gene.list %in% filter.names)]

  return(txi)
}


# create a function to run the DESeq2 analysis using the txi and the dds objects
run.deseq2 <- function(txi, variable, outdir) {
    # filter the txi object
    txi.filter <- filter.tpm(txi, tpm_cutoff, coldata, variable)
    dds <- DESeqDataSetFromTximport(txi, coldata, ~ variable)
    # if the variable == "Compound", then use the Compound column in the coldata file, if the variable == "Time", then use the Time column in the coldata file, if the variable == "Dose",
    # then use the Dose column in the coldata file, if the variable == "Treatment", if the variable == "Disease", then use the Disease column in the coldata file.
    if(variable == "Compound") {
      dds <- DESeqDataSetFromTximport(txi.filter, coldata, ~ Compound)
      # relevel the dds object to DMSO as the ref 
      dds$Compound <- relevel(dds$Compound, ref = "DMSO")
    } else if(variable == "Time") {
      dds <- DESeqDataSetFromTximport(txi.filter, coldata, ~ Time)
      # relevel the dds object to 0 as the ref
      dds$Time <- relevel(dds$Time, ref = "0")

    } else if(variable == "Dose") {
      dds <- DESeqDataSetFromTximport(txi.filter, coldata, ~ Dose)

      # relevel the dds object to 0 as the ref
      dds$Dose <- relevel(dds$Dose, ref = "0")

    } else if(variable == "Treatment") {
      dds <- DESeqDataSetFromTximport(ttxi.filterxi, coldata, ~ Treatment)
      # relevel the dds object to naive as the ref
      dds$Treatment <- relevel(dds$Treatment, ref = "naive")
      

    } else if(variable == "Disease") {
      dds <- DESeqDataSetFromTximport(txi.filter, coldata, ~ Disease)

        # relevel the dds object to healthy as the ref
        dds$Disease <- relevel(dds$Disease, ref = "healthy")
    } else {
      stop("Please specify a variable with -v or --variable")
    }

    featureData <- data.frame(gene = rownames(dds))
    colnames(dds) <- coldata$SampleID
    mcols(dds) <- DataFrame(mcols(dds), featureData)
    dds <- DESeq(dds, test="Wald", fitType="parametric")
    res <- results(dds)

    # save the dds and res objects to Rdata files
    save(dds, file = file.path(outdir, paste0(variable, "_dds.Rdata")))
    save(res, file = file.path(outdir, paste0(variable, "_res.Rdata")))
    return(res)  
}

# run the run.deseq2 function
res <- run.deseq2(txi, opt$variable, opt$out)

# create a resOrdered object that uses the res object and the padj column
resOrdered <- res[order(res$padj),]

# assign a level of significance to the resOrdered object where the log2FoldChange is greater than 1 and the padj is less than 0.05 if it is upregulated, and the log2FoldChange is less than -1 and the padj is less than 0.05 if it is downregulated
resOrdered$diffexpressed <- "NotSig"
resOrdered$diffexpressed[resOrdered$padj < 0.05 & resOrdered$log2FoldChange >= 1] <- "Upregulated"
resOrdered$diffexpressed[resOrdered$padj < 0.05 & resOrdered$log2FoldChange <= -1] <- "Downregulated"
# write the resOrdered object to a file in the output directory
write.table(resOrdered, file = file.path(opt$out, paste0(opt$variable, "_resOrdered.txt")), sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)


# Count the number of up and down regulated genes in the resOrdered object and plot them in a ggplot barplot where the red bars are upregulated and the blue bars are downregulated
# create a heatmap
# create a gene set enrichment analsyis

