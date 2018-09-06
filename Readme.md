# GATK3 Human Exome Pipeline (github release v 1.4)

**Please read this important information before running the app.**

## What are typical use cases for this app?

Use this app when you have mapped exome reads to the human genome, and want to identify variants (SNPs and indels) using the GATK 3.x best practices method.



## What does this app do?

This app implements the GATK 3.x best practices pipeline. It takes a BAM file as input, refines the BAM (by deduplicating, realigning, and recalibrating) before it calls variants using GATK HaplotypeCaller. It outputs the refined mappings as well as the called variants. 
NOTE: This app does not perform variant recalibration.

This app has been modified to work with exome and custom panels.

The steps applied in the app are:
1. mark duplicates (using picard)
2. GATK indel realignment (uses known indels from 1000genomes phase 1 and Mills and 1000G gold standard indels
3. GATK base recalibration (using the known sites mentuoned above and from dbSNP (137).
4. Run GATK Haplotype caller. 
* For WES samples a BED file of the exome capture kit is used to focus variant calling (with padding of 100bp - unless stated otherwise in the parameters).
* Non-WES samples are not filtered by region and do not have any padding.

The reference sequence is obtained by reading the BAM header and is downloaded from the relevant public project within the app.


## What data are required for this app to run?

1. GATK JAR file
 * This app is only a wrapper for the GATK 3.x software, and requires that you appropriately license and obtain that software yourself.
 * After licensing GATK, you should have received a file with the `GenomeAnalysisTK` prefix and the `.jar` suffix, such as `GenomeAnalysisTK.jar` or `GenomeAnalysisTK-3.4-0.jar`. 

 * This file must be stated as an input.

2. BAM file 
This app requires a coordinate-sorted BAM file (`*.bam`). The app automatically detects the reference genome (hg19, GRCh37/b37, or GRCh37+decoy/hs37d5) based on the BAM file header, and uses the appropriate GATK resources (dbSNP and known indels).

3. BED file.
The exome capture bed file must be supplied for WES samples.


## What does this app output?
1. BAM (and index)
This app outputs the refined (deduplicated, realigned, and recalibrated) mappings in BAM format (`*.bam`), as well as the associated BAM index (`*.bai`).

2. Mark Duplicates Output Metrics 
The Mark duplicates output metrics file which is used to produce run-wide QC.

3. VCF (and index)
The app also outputs a _genotyped_ VCF file (`*.vcf.gz`) and its associated tabix index (`*.vcf.gz.tbi`), or an intermediate gVCF file (`*.g.vcf.gz`) and its associated tabix index (`*.g.vcf.gz.tbi`), or all of the above. This behavior depends on the "Output format" option. The option works as follows:

 * When set to `vcf`, the app runs GATK HaplotypeCaller in regular (_genotyped_ VCF) mode; this calls variants and outputs only the locations of variation.
 * When set to `gvcf`, the app runs GATK HaplotypeCaller in gVCF mode; this outputs information for all locations, including sections which lack variation. The gVCF is an intermediate file that can be later used as input to GATK GenotypeGVCFs, which can take multiple gVCF files (from multiple samples) and genotype them, creating a cohort-level genotyped VCF. For more information, consult [this GATK article](http://gatkforums.broadinstitute.org/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode).
 * When set to `both`, the app runs GATK HaplotypeCaller in gVCF mode, producing a gVCF file. Subsequently, it runs GATK GenotypeGVCFs to genotype the gVCF into a regular VCF.

