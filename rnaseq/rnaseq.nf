#! usr/bin/env nextflow

/*
 * pipeline input parameters
 */

params.reads = /path/to/reads/*.fastq.gz
params.transcriptome_index = /path/to/transcriptome/index/
params.transcriptome_file = /path/to/transcriptome/hg38.fa
params.transcriptome_gtf = /path/to/transcriptome/hg38.gtf
params.multiqc = /path/to/output/directory/for/multiqc
params.quantification = /path/to/salmon/quantification/output


process trimGaloreSingleEndReads {

	input:
	output:

	script:
	"""

	# write command here

	"""

}

process multiqc_reads {

	input:
        output:

        script:
        """

        # write command here

        """


}


process salmonQuant {

	input:
        output:

        script:
        """

        # write command here

        """




}

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    transcriptome: ${params.transcriptome_file}
    reads        : ${params.reads}
    multiqc_dir  : ${params.multiqc}
    salmon_quant : ${params.quantification}
    outdir       : ${params.outdir}
    """
    .stripIndent()

workflow {

	trimGaloreSingleEndReads(input_data)
	multiqc_reads(input_data)
	salmonQuant(input_data)

}
