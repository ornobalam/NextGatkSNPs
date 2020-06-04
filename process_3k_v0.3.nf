#!/usr/bin/env nextflow

// ######################################
// #    Rafal Gutaker, NYU, 2017        #
// #    Pipeline that performs:         #
// #    - mapping to reference          #
// #    - sorting reads		        #
// #    - removal PCR duplicates        #
// #    - calling haplotypes            #
// #	- storing g.vcfs in GenomicsDB	#
// #					#
// #    Updated to GATK4 best practices #
// #	by Ornob Alam, NYU, 2020	#
// #                                    #
// ######################################


params.help=""
params.exe="local"

if (params.help) {
        log.info " "
        log.info "This is Nextflow pipeline for processing rice 3k dataset. Most parameters are stored in file process_3k.conf"
        log.info " "
        log.info "Usage: nextflow run process_3k.nf -c process_3k.conf --ref --list --chrom --(other options)"
        log.info "Options:"
        log.info "--help\t[BOOLEAN]\tShow this help message"
        log.info "--ref\t[STRING]\tPath to the indexed referece fasta file [OBLIGATORY]"
        log.info "--list\t[STRING]\tPath to the file with run IDs and fastq files to be processed [OBLIGATORY]"
        log.info "--exe\t[STRING]\tExecutor mode, -local- or -slurm- [DEFUALT: local]"
        log.info " "
        exit 1
}

// Initalize Input
DAT = file("${params.list}")
REF = file("${params.ref}")
CHR = file("${params.chrom}")

// Create Input Channel
SampleData = Channel.fromPath("${DAT}").splitCsv(header: ['SID','RID','P1','P2'], skip: 0, by:1, sep:",")
SampleData2 = Channel.fromPath("${CHR}").splitCsv(header: ['CHR'], skip: 0, by:1)

// ########### Run BWA #############
process bwa {

	errorStrategy 'finish'	

        executor = "${params.exe}"
	if ("${params.exe}" == "slurm")
	{
        	clusterOptions = "--cpus-per-task=${params.bwa_cpu} --time=${params.bwa_rt} --mem=${params.bwa_vmem}"
	}

        input:
         val(RAW) from SampleData

        output:
         set val(ID), val(SID), file({ "${MAP}" }) into aligned

        script:

	 SID = "${RAW.SID}"
	 ID = "${RAW.RID}"
         RG = "\'@RG\\tID:${RAW.P1}\\tSM:${RAW.SID}\\tPL:ILLUMINA\'"
         MAP = "${RAW.RID}.sam"

         P1 = file("${RAW.P1}")
         P2 = file("${RAW.P2}")

         """
         ${params.bwa_bin} mem ${params.bwa_param} -M -R ${RG} ${REF} ${P1} ${P2} > ${MAP}
	 """
}


// ########### Sort SAM ###############
process sort {

	errorStrategy 'finish'

        executor = "${params.exe}"
        if ("${params.exe}" == "slurm")
        {
        	clusterOptions = "--cpus-per-task=${params.srt_cpu} --time=${params.srt_rt} --mem=${params.srt_vmem}"
	}

        input:
         set val(ID), val(SID), file(SAM) from aligned

        output:
         set val(ID), val(SID), file({ "${SBAM}" }) into sorted

        script:

         SBAM = "${ID}.sorted.bam"
         del_sam = SAM.getName()

         """
         ${params.srt_bin} INPUT=${SAM} OUTPUT=${SBAM} SORT_ORDER=coordinate
         """
}

// ####### MERGE FILE CHANNELS #########


sorted
	.map {id, sid, file -> tuple(sid, file) }
	.groupTuple()
	.set { merg_sid }


// ###### Remove Read Duplicates #######
process rmdup {

	errorStrategy 'finish'

        executor = "${params.exe}"
        if ("${params.exe}" == "slurm")
        {
		clusterOptions = "--cpus-per-task=${params.ddp_cpu} --time=${params.ddp_rt} --mem=${params.ddp_vmem}"
	}

        input:
         set val(SID), file(SBAM) from merg_sid

        output:
         set val(SID), file({ "${RMDUP}" }) into rmdup

        script:

	 ALL_IN = SBAM.collect { "INPUT=$it" }.join(' ')
         RMDUP = "${SID}.dedup.bam"
         RMET = "${SID}.dedup.met"
//         del_sbam = SBAM.getName()

         """
         ${params.ddp_bin} ${ALL_IN} OUTPUT=${RMDUP} METRICS_FILE=${RMET}
         """
}


// ########### Check BAM ##########

process checkbam {

	errorStrategy 'ignore'

	executor = "${params.exe}"
	if ("${params.exe}" == "slurm")
	{
		clusterOptions = "--cpus-per-task=${params.vb_cpu} --time=${params.vb_rt} --mem=${params.vb_vmem}"
	}


	input:
	set val(SID), file(RMDUP) from rmdup

	output:
	set val(SID), file(RMDUP) into bchecked

	script:

	"""
	${params.vb_bin} INPUT=${RMDUP} ${params.vb_param}
	"""
}

// ########### Index BAM ##########
process index {

	errorStrategy 'finish'

        executor = "${params.exe}"
        if ("${params.exe}" == "slurm")
        {
	        clusterOptions = "--cpus-per-task=${params.idx_cpu} --time=${params.idx_rt} --mem=${params.idx_vmem}"
	}


        input:
         set val(SID), file(RMDUP) from bchecked

        output:
         set val(SID), file({ "${RMDUP}" }), file({ "${BAI}" }) into index

        script:

         BAI = "${SID}.dedup.bai"

         """
         ${params.idx_bin} INPUT=${RMDUP} OUTPUT=${BAI}
         """

}

// ######## HaplotypeCaller ########

process haplocall {

	errorStrategy 'finish'

	executor = "${params.exe}"
	if ("${params.exe}" == "slurm")
	{
		clusterOptions = "--cpus-per-task=${params.hc_cpu} --time=${params.hc_rt} --mem=${params.hc_vmem}"
	}

 	input:
 	set val(SID), file(RMDUP), file(BAI) from index
 
 	output:
 	file({ "${CALL}" }) into (called)
	file({ "${VCI}" }) into (vcfindex)
 
 	script:
 	
 	CALL = "${SID}.g.vcf"
	VCI = "${SID}.g.vcf.idx"
 
 	"""
 	${params.hc_bin} ${params.hc_param} -R ${REF} -I ${RMDUP} -O ${CALL}	
 	"""
 
}

// ########## Check GVCF ###########

process checkvcf {

	errorStrategy 'ignore'

	executor = "${params.exe}"
	if ("${params.exe}" == "slurm")
	{
		clusterOptions = "--cpus-per-task=${params.vv_cpu} --time=${params.vv_rt} --mem=${params.vv_vmem}"
	}

	input:
	file(CALL) from called
	
	output:
	file(CALL) into vchecked

	script:

	"""
	${params.vv_bin} ${params.vv_param} -R ${REF} -V ${CALL}
	"""
}


// ######## Convert into lists ##########

vchecked
	.toList()
	.set{ listed }

vcfindex
	.toList()
	.set{ vcflist }

// ######### CombineGVCFs ##########

process combine {
 
	errorStrategy 'finish'

 	executor = "${params.exe}"
 	if ("${params.exe}" == "slurm")
 	{
 		clusterOptions = "--cpus-per-task=${params.cg_cpu} --time=${params.cg_rt} --mem=${params.cg_vmem}"
 	}

	input:           
	val(CHROM) from SampleData2
	file(CALL) from listed
	file(VCI) from vcflist

	output:
	file( "${DB}" )                                                                                                                                 

	script:
	CHR = "${CHROM.CHR}"
	COMB_IN = CALL.collect { "--variant $it" }.join(' ')
	DB = "../../../../${CHROM.CHR}_db"

	"""
	${params.cg_bin} -R ${REF} ${COMB_IN} --genomicsdb-workspace-path $DB -L ${CHR}
	""" 
}
