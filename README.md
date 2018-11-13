# GATK3 pre-processing v 1.0
**Please read this important information before running the app.**

## What are typical use cases for this app?
To perform GATK preprocessing steps as per GATK best practice. 
The refined BAM file output can be used for variant calling.

## What does this app do?
This app implements the GATK 3.x best practices pipeline. It takes a BAM file as input, refines the BAM (by deduplicating, realigning, and recalibrating) outputting the refined mappings. 

The steps applied in the app are:
1. Picard mark duplicates (if selected)
2. GATK indel realignment (uses known indels from 1000genomes phase 1 and Mills and 1000G gold standard indels
3. GATK base recalibration (using the known sites mentuoned above and from dbSNP (137).

The reference sequence is obtained by reading the BAM header and is downloaded from the relevant public project within the app.

## What data are required for this app to run?
1. BAM file 
This app requires a coordinate-sorted BAM file (`*.bam`). The app automatically detects the reference genome (hg19, GRCh37/b37, or GRCh37+decoy/hs37d5) based on the BAM file header, and uses the appropriate GATK resources (dbSNP and known indels).

2. A boolean flag if mark duplicates should be performed

## What does this app output?
1. BAM (and index)
This app outputs the refined (deduplicated, realigned, and recalibrated) mappings in BAM format (`*.bam`), as well as the associated BAM index (`*.bai`).

2. Mark Duplicates Output Metrics 
The Mark duplicates output metrics file which is used to produce run-wide QC.

## Custom modifications
* The app was modified by Viapath Genome Informatics
