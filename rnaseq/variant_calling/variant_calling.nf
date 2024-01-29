#!/usr/bin/ env nextflow

/*
* Who: Michael Skaro
* What: Variant calling pipeline for RNAseq data following the nfcore/rnaseq pipeline
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

/* First process will be to index the reference genome with star */

process indexGenome {
  publishDir "${params.outputDir}/star_index", mode: 'copy'

  input:
  file genome

  output:
  file 'genome.fa'

  script:
  """
  STAR --runMode genomeGenerate --genomeDir ${params.outputDir}/star_index --genomeFastaFiles ${genome} --sjdbGTFfile ${params.gtf} --runThreadN ${params.threads}
  """
}

/* Second process will be to QC the reads we have downloaded from the SRA */

process fastqc {
  publishDir "${params.outputDir}/fastqc", mode: 'copy'

  input:
  file '*.fastq' from inputDir

  output:
  file '*.html' into fastqcDir

  script:
  """
  fastqc -o ${fastqcDir} ${input}
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
  multiqc ${fastqcDir} -o ${multiqcDir}
  """
}

/* Align the reads to the reference genome, create a sorted bam file, and index it */

process starAlign {
  publishDir "${params.outputDir}/star_align", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.fastq' from inputDir

  output:
  file 'Aligned.sortedByCoord.out.bam' into starAlignDir

  script:
  """
  STAR --genomeDir ${indexGenome}/star_index --readFilesIn ${input} --runThreadN ${params.threads} --outFileNamePrefix ${starAlignDir}/
  samtools sort -@ ${params.threads} -o ${starAlignDir}/Aligned.sortedByCoord.out.bam ${starAlignDir}/Aligned.out.sam
  samtools index -@ ${params.threads} ${starAlignDir}/Aligned.sortedByCoord.out.bam
  """
}

/* picard to mark duplicates */

process picard {
  publishDir "${params.outputDir}/picard", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.fastq' from inputDir

  output:
  file 'Aligned.sortedByCoord.out.bam' into picardDir

  script:
  """
  java -jar picard.jar MarkDuplicates I=${starAlignDir}/Aligned.sortedByCoord.out.bam O=${picardDir}/output.bam M=${picardDir}/output.metrics.txt
  """
}



/* Pass the BAM files to RSEM to estimate the expression levels */
/* might have screwed up this part, but I think it's close to being correct */
process rsem {
  publishDir "${params.outputDir}/rsem", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.fastq' from inputDir

  output:
  file '*.genes.results' into rsemDir

  script:
  """
  rsem-calculate-expression --bam --estimate-rspd --append-names --no-bam-output --seed 12345 --num-threads ${params.threads} ${input} ${indexGenome}/star_index ${rsemDir}/
  """
}

/* pass the BAM files to GATK to call the variants */

process gatk {
  publishDir "${params.outputDir}/gatk", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.fastq' from inputDir

  output:
  file '*.vcf' into gatkDir

  script:
  """
  gatk HaplotypeCaller -R ${genome} -I ${starAlignDir}/Aligned.sortedByCoord.out.bam -O ${gatkDir}/output.vcf
  """
}

/* Pass the VCF files to bcftools to filter the variants */

process bcftools {
  publishDir "${params.outputDir}/bcftools", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.vcf' from gatkDir

  output:
  file '*.vcf' into bcftoolsDir

  script:
  """
  bcftools filter -i 'QUAL > 20' ${gatkDir}/output.vcf -o ${bcftoolsDir}/output.vcf
  """
}

/* Pass the VCF files to snpEff to annotate the variants */

process snpEff {
  publishDir "${params.outputDir}/snpEff", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.vcf' from bcftoolsDir

  output:
  file '*.vcf' into snpEffDir

  script:
  """
  snpEff -v -stats ${snpEffDir}/snpEff_summary.html -c ${params.gtf} -s ${snpEffDir}/snpEff_summary.html -no-downstream -no-upstream -no-utr -no-intergenic \\
  -no-intron -no-utr -no-nc -no-ns -no-aa -no-ann -no-gene -no-effect -no-quiet -no-summary -no-stats -no-verbose -no-warnings -no-ver 
  """
} 

/* Annotate the variants with the gene names */

process annovar {
  publishDir "${params.outputDir}/annovar", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.vcf' from snpEffDir

  output:
  file '*.vcf' into annovarDir

  script:
  """
  table_annovar.pl ${snpEffDir}/output.vcf ${params.outputDir}/annovar -buildver hg38 -out ${annovarDir}/output -remove -protocol refGene -operation g -nastring . -vcfinput
  """
}

/* Use sift and polyphen to predict the effects of the variants */

process siftPolyphen {
  publishDir "${params.outputDir}/sift_polyphen", mode: 'copy'

  input:
  file genome from indexGenome
  file '*.vcf' from annovarDir

  output:
  file '*.vcf' into siftPolyphenDir

  script:
  """
  perl ${params.outputDir}/annovar/variants_function -vcfinput -buildver hg38 -out ${siftPolyphenDir}/output \\
  -protocol refGene -operation g -nastring . -vcfinput -polish -otherinfo -sift s -sift p -dbnsfp 3.5a \\
  -hgvs -exonicsplicing -genericdb ${params.outputDir}/annovar/hg38_genericdb.txt -genericdbfile ${params.outputDir}/annovar/hg38_genericdb.txt 
  """
}


/* run the pipeine using the workflow */

workflow{ 
    indexGenome
    fastqc
    multiqc
    starAlign
    picard
    rsem
    gatk
    bcftools
    snpEff
    annovar
    siftPolyphen   
}
