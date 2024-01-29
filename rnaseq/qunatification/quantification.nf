#!/usr/bin/ env nextflow

/*
* Who: Michael Skaro
* What: Nexflow script for quantification of gene expression data using salmon
* When: January 2024
*/

/* Declare the global parameters */

params {
  // The input directory containing the fastq files
  inputDir = 'data/fastq'

  // The output directory
  outputDir = 'results'

  // The reference genome
  genome = 'data/genome.fa'

  // The GTF file
  gtf = 'data/genes.gtf'

  // The number of threads to use
  threads = 8

  // The memory to use
  memory = '16GB'
}

/* Generate the index for the reference genome using salmon */

process indexGenome {
  publishDir "${params.outputDir}/salmon_index", mode: 'copy'

  input:
  file genome

  output:
  file 'genome.fa'

  script:
  """
  salmon index -t ${genome} -i ${params.outputDir}/salmon_index
  """
}

/* Run fastqc on the input fastq files */

process fastqc {
  publishDir "${params.outputDir}/fastqc", mode: 'copy'

  input:
  file '*.fastq' from inputDir

  output:
  file '*.html' into fastqcDir

  script:
  """
  fastqc -o ${params.outputDir}/fastqc ${input}
  """
}

/* Use multiqc to aggregate the fastqc results */

process multiqc {
  publishDir "${params.outputDir}/multiqc", mode: 'copy'

  input:
  file '*.html' from fastqcDir

  output:
  file 'multiqc_report.html' into multiqcDir

  script:
  """
  multiqc ${fastqcDir} -o ${params.outputDir}/multiqc
  """
}

/* Trim the input fastq files using trimGalore */

process trimGaloreSingleEnd {
  publishDir "${params.outputDir}/trimmed_fastq", mode: 'copy'

  input:
  file '*.fastq' from inputDir

  output:
  file '*.fq' into trimmedDir

  script:
  """
  trim_galore --fastqc --output_dir ${params.outputDir}/trimmed_fastq ${input}
  """
}

/* Create a process for paired-end reads trimGalore */


process trimGalorePairedEnd {
  publishDir "${params.outputDir}/trimmed_fastq", mode: 'copy'

  input:
  file '*.fastq' from inputDir

  output:
  file '*.fq' into trimmedDir

  script:
  """
  trim_galore --fastqc --output_dir ${params.outputDir}/trimmed_fastq --paired ${input}
  """
}

/* Quantification of gene expression using salmon from single end reads */

process salmonQuantSingle {
  publishDir "${params.outputDir}/salmon_quant", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.fq' from trimmedDir

  output:
  file '*.sf' into salmonDir

  script:
  """
  salmon quant -i ${indexGenome} -l A -r ${input} -p ${params.threads} -o ${salmonDir}
  """
}



/* Quantify the gene expression using salmon from the trimmed reads */

process salmonQuantPaired {
  publishDir "${params.outputDir}/salmon_quant", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.fq' from trimmedDir

  output:
  file '*.sf' into salmonDir

  script:
  """
  salmon quant -i ${indexGenome} -l A -1 ${input} -2 ${input} -p ${params.threads} -o ${salmonDir}
  """
}

/* Create a qua

/* Create a multiqc report for the trimmed reads the salmon results */

process multiqcSalmon {
  publishDir "${params.outputDir}/multiqc_salmon", mode: 'copy'

  input:
  file '*.html' from fastqcDir
  file '*.sf' from salmonDir

  output:
  file 'multiqc_report.html' into multiqcSalmonDir

  script:
  """
  multiqc ${params.outputDir}/trimmed_fastq ${salmonDir} -o ${params.outputDir}/multiqc_salmon
  multiqc ${params.outputDir}/multiqc_salmon -o ${params.outputDir}/multiqc_salmon
  """
}

/* Run the workflow */
/* Create a check if the reads are paired-end or single-end */
workflow {
  if (params.pairedEnd) {
    call indexGenome
    call fastqc
    call multiqc
    call trimGalorePairedEnd
    call salmonQuantPaired
    call multiqcSalmon
  } else {
    call indexGenome
    call fastqc
    call multiqc
    call trimGaloreSingleEnd
    call salmonQuantSingle
    call multiqcSalmon
  }
}
