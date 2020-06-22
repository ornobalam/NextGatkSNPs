#!/usr/bin/env nextflow

// ######################################
// #    Rafal Gutaker, NYU, 2018        #
// #    Pipeline that genotypes from    #
// #    GenomicsDB files		#
// #    				#
// #    Updated to GATK4 best practices #
// #	by Ornob Alam, NYU, 2020	#
// ######################################


params.help=""
params.exe="local"

if (params.help) {
        log.info " "
        log.info "This is Nextflow pipeline for genotyping SNPS in rice 3k dataset. Most parameters are stored in file call_3k.conf"
        log.info " "
        log.info "Usage: nextflow run call_3k.nf -c call_3k.conf --ref --chrom --(other options)"
        log.info "Options:"
        log.info "--help\t[BOOLEAN]\tShow this help message"
        log.info "--ref\t[STRING]\tPath to the indexed referece fasta file [OBLIGATORY]"
        log.info "--chrom\t[STRING]\tPath to the file with chromosome IDs [OBLIGATORY]"
        log.info "--exe\t[STRING]\tExecutor mode, -local- or -slurm- [DEFUALT: local]"
        log.info " "
        exit 1
}

// Initalize Input
CHROM = file("${params.chrom}")
REF = file("${params.ref}")

// Create Input Channel
SampleData = Channel.fromPath("${CHROM}").splitCsv(header: ['CHR'], skip: 0, by:1)


// ########### Genotype GVCF ###############
process genotype {

	errorStrategy 'retry'  
        maxRetries 2
        maxErrors 2

        executor = "${params.exe}"
        if ("${params.exe}" == "slurm")
        {
        	clusterOptions = "--cpus-per-task=${params.gv_cpu} --time=${params.gv_rt} --mem=${params.gv_vmem}"
	}

        input:
        val(CHROM) from SampleData

        output:
         file({ "${CALLD}" })

        script:
	CHR = "${CHROM.CHR}"
        CALLD = "${CHR}.vcf"
	DB = "../../../../${CHROM.CHR}_db"

         """
         ${params.gv_bin} ${params.gv_param} -V gendb://$DB -R ${REF} -L ${CHR} -O ${CALLD}
         """
}

