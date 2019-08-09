# NextGatkSNPs
Nextflow pipelines for short read mapping and SNP calling - complies with GATK3.8 best practices. 

This is NextFlow pipeline that runs all processes to map short reads in fastq format to reference genome, filter and call SNPs. Branch princeV0.1 was optimized to work with slurm submission system. Please note that you will have to install bioinformatic packages (bioconda installations worked for me) and change configuration files (*.conf) to reflect paths to your installations. In particular you will need BWA, PICARD, GATK, BGZIP and TABIX. This pipeline complies with GATK3.8 best practices and is unlikely to work with version 4.0+ without serious changes.

The process.nf pipeline will run following processes:
1 - map all fastq short reads to a reference genome independently
2 - sort all sam files and generate bam files
3 - merge mutliple bam files for the same indvidual
4 - remove PCR duplicates
5 - check if bam files are correct
6 - index merged bam file
7 - call haplotypes in g.vcf format
8 - check if vcf files are correct
9' - compres and index indvidual g.vcf files
9" - combine multiple g.vcf files

This pipeline requires a configuration file, process.conf to find packages and modify cpu/time/memory usage. Additionally, it requires a comma separated "batch" file, which will include information about your fastq files. Finally, it will need a reference genome in fasta format indexed with bwa. 

Batch files format is as follows: 
1st column - individual names
2nd column - library name
3rd column - fastq files with P1 reads
4th column - fastq files with P2 reads

If all files are present, you can run this pipeline as follows:

nextflow run -c process_3k_v0.3.conf process_3k_v0.3.nf --ref ref.fa --list batch.csv
