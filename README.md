### Who: Michael Skaro PhD
### What: RNAseq pipeline project
### When: 2024-01-24
### Why: The purpose of this project is to learn how to use the NextFlow workflow manager. The project will be a RNAseq pipeline that will take in fastq files and output a list of differentially expressed genes from multiple time points and doses on two cell lines.

##### ----------------------------------------------------------------------------------------------------------------------------------------------------------

#### 1. Download the data [X]
#### 2. Quality control with trimming with trimgalore, multiqc []
#### 3. Alignment with STAR, featureCounts []
#### 4. Quantify with salmon, and RSEM []
#### 5. Differential expression with DESeq2 against DMSO for each drug at
####    each time point. []

#### Overall design:RNAseq of MCF7 and LNCaP cells each treated with
#### DMSO, 3 drugs, and the 3 combinations of the 3 drugs, in trplicate 
#### well as a triple combination in LNCaP, at 0, 3, 6, 9, 12, and 24 
#### hours

#### The drugs are: Tamoxifen, Mefloquine, and Withaferin A
#### The cell lines are: MCF7 and LNCaP

#### Data format from methods section of PMID: 32945258 : The TruSeq Stranded mRNA Library Prep Kit (RS-122–2101/RS-122–2102) was then used to prepare the samples for 30 million reads of single end sequencing (100 bp) with the Illumina HiSeq2500. Three replicates were used for the combination experiments and two replicates for the dose experiments.

#### The data is sinlge end, 100 bp reads, This modify the pipeline from a paired end to asingle end stranded approach. 
#### We will put in a trim step to remove the polyA tails, and then proceed with the rest of the pipeline as normal. 


