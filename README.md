# NextGatkSNPs

[CURRENTLY UPDATING]

Nextflow pipelines for short read mapping and SNP calling - complies with GATK4 best practices. 

This is NextFlow pipeline that runs all processes to map short reads in fastq format to reference genome, 
filter and call SNPs. Branch princeV0.1 was optimized to work with slurm submission system. Please note 
that you will have to install bioinformatic packages (bioconda installations worked for me) and change 
configuration files (.conf) to reflect paths to your installations. In particular you will need BWA, 
PICARD,  and GATK. This pipeline complies with GATK4 best practices, and has been adapted from a previous
version created by Rafal Gutaker that complied with GATK3.8 best practices.

## PREPARE FILES FOR CALLING SNPS

The **process.nf** pipeline will run following processes:

  - map all fastq short reads to a reference genome independently,
  - sort all sam files and generate bam files,
  - merge mutliple bam files for the same indvidual,
  - remove PCR duplicates,
  - check if bam files are correct,
  - index merged bam file,
  - call haplotypes in g.vcf format,
  - check if vcf files are correct,
  - compare and index indvidual g.vcf files,
  - combine multiple g.vcf files in the GenomicsDB format.

This pipeline requires a configuration file, process.conf to find packages and modify cpu/time/memory usage. 
Additionally, it requires a comma separated "batch" file, which will include information about your fastq files. 
Finally, it will need a reference genome in fasta format indexed with bwa, and a text file with a single column
with all the chromosome IDs.

Batch files format is as follows: 

  - column 1: individual names,
  - column 2: library name,
  - column 3: fastq files with P1 read,
  - column 4: fastq files with P2 reads.
  

If all files are present, you can run this pipeline as follows:

```
nextflow run -c process_3k_v0.3.conf process_3k_v0.3.nf --ref ref.fa --list batch.csv --chrom chrom.txt
```

Note that you can add more individuals to the GenomicsDB databases by changing the GenomicsDBImport argument from --genomicsdb-workspace-path to --genomicsdb-update-workspace-path.

## CALL SNPS

Now that you have GenomicsDB files for each chromosome, you can use them to call SNPs

The **call.nf** pipeline will run following process in parallel for each chromosome:

  - genotype g.vcf haplotypes.
  
This pipeline requires a configuration file, call.conf to find packages and modify cpu/time/memory usage. Additionally, it requires a chrom.txt file with a single column with all the chromosome IDs. Finally, it will need a reference genome in fasta format indexed with bwa.


If all files are present, you can run this pipeline as follows:

```
nextflow run -c call.conf call.nf --ref ref.fa --chrom chrom.txt
```

