#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

# check not all steps are being skipped otherwise we get a 'refined.bam' that hasn't had any preprocessing.
if [[ "$skip_markduplicates" == "true" ]] && [[ "$skip_BQSR" == "true" ]] && [[ "$skip_indelrealignment" == "true" ]]
then
exit
fi

#
# Fetch inputs
#
#mark-section "downloading inputs"

dx-download-all-inputs --parallel

# Show all the java versions installed on this worker
# Show the java version the worker is using
echo $(java -version)

#
# Use java7 as the default java version. 
# If java7 doesn't work with the GATK version (3.6 and above) then switch to java8 and try again.
#
update-alternatives --set java /usr/lib/jvm/java-7-openjdk-amd64/jre/bin/java
java -jar GenomeAnalysisTK.jar -version || (update-alternatives --set java /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java && java -jar GenomeAnalysisTK.jar -version)


#
# Calculate 80% of memory size, for java
#
head -n1 /proc/meminfo | awk '{print int($2*0.8/1024)}' >.mem_in_mb.txt
java="java -Xmx$(<.mem_in_mb.txt)m"


#
# Detect and download the appropriate human genome and related reference files
#
#mark-section "detecting reference genome"
samtools view -H "$sorted_bam_path" | grep ^@SQ | cut -f1-3 | md5sum | cut -c1-32 >.genome-fingerprint.txt
# parse the BAM header to look for all the reference sequences.
# use cut -f1-3 to reformat and calculate a checksum from this info. use cut -c1-32 to reformat and save this checksum into genome-fingerprint.txt
# this info is used to determine which reference genome was used for alignment and so the correct reference files can be downloaded below.
case "$(<.genome-fingerprint.txt)" in
  9220d59b0d7a55a43b22cad4a87f6797)
    genome=b37
    subgenome=b37
    ;;
  45340a8b2bb041655c6f6d4f9985944f)
    genome=b37
    subgenome=hs37d5
    ;;
  2f23b2f7c9731db07f0d1c8f9bc8c9d9)
    genome=hg19
    subgenome=hg19
    ;;
  53a8d91e94b765bd69c104611a2f796c)
    genome=grch38
    subgenome=grch38 # No alt analysis
    dx-jobutil-report-error "The GRCh38 reference genome is not supported by this app."
    ;;
  *)
    echo "Non-matching human genome. The input BAM contains the following chromosomes (names and sizes):"
    samtools view -H "$sorted_bam_path" | grep ^@SQ | cut -f1-3
    dx-jobutil-report-error "The reference genome of the input BAM file did not match any of the known human ones. Additional diagnostic information has been provided in the job log."
    ;;
esac

#mark-section "downloading reference genome"
APPDATA=project-B6JG85Z2J35vb6Z7pQ9Q02j8
dx cat "$APPDATA:/misc/gatk_resource_archives/${subgenome}.fasta-index.tar.gz" | tar zxf -


#echo $region_opts

# Fetch GATK resources (these have been prepared in archives)
#mark-section "fetching GATK resources"
dx cat "$APPDATA:/misc/gatk_resource_archives/gatk.resources.${genome}.tar" | tar xf -
known1="1000G_phase1.indels.${genome}.vcf.gz"
known2="Mills_and_1000G_gold_standard.indels.${genome}.vcf.gz"
dbsnp="dbsnp_137.${genome}.vcf.gz"


#
# Run Picard MarkDuplicates
#
# if mark duplicates is not skipped run the tool
if [[ "$skip_markduplicates" != "true" ]]
then
  mark-section "marking duplicates"
  $java -jar picard.jar MarkDuplicates I="$sorted_bam_path" O=deduplicated.bam M=output.metrics CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT $extra_md_options
  rm -f "$sorted_bam_path"
  realign_input="deduplicated.bam"
# if mark suplicates is skipped we need to set the inputs for the next tool to the files input to mark duplicates. This also requires an indexed BAM so need to create this index.
else
  mark-section "mark duplicates skipped"
  samtools index "$sorted_bam_path"
  realign_input="$sorted_bam_path"
fi

#
# Run GATK indel realignment
#
if [[ "$skip_indelrealignment" != "true" ]]
then
  mark-section "realigning indels"
  $java -jar GenomeAnalysisTK.jar -nct 1 -T RealignerTargetCreator -R genome.fa -I "$realign_input" -o realign.intervals -known $known1 -known $known2 $extra_rtc_options
  $java -jar GenomeAnalysisTK.jar -nct 1 -T IndelRealigner -R genome.fa -I "$realign_input" -targetIntervals realign.intervals -known $known1 -known $known2 -o realigned.bam $extra_ir_options
  realigned_bam="realigned.bam"
  samtools index "$realigned_bam"
  rm -f "$realign_input"
# if this step is skipped need to set the BQSR inputs to the inputs to indel realignment
else
  mark-section "not performing indel realignment"
  realigned_bam="$realign_input"
  # create path to index 
  realigned_bai="${realign_input/bam/"bai"}"
  fi

#
# Run GATK base recalibration
#
if [[ "$skip_BQSR" != "true" ]]
then
  mark-section "recalibrating base quality scores"
  $java -jar GenomeAnalysisTK.jar -nct 1 -T BaseRecalibrator -R genome.fa -I $realigned_bam -o recal.grp -knownSites $dbsnp -knownSites $known1 -knownSites $known2 $extra_br_options
  $java -jar GenomeAnalysisTK.jar -nct 1 -T PrintReads -R genome.fa -I $realigned_bam -BQSR recal.grp -o recal.bam $extra_pr_options
  rm -f realigned.bam
  recal_bam="recal.bam"
  recal_bai="recal.bai"
# if this step is skipped need to set the files to upload as the inputs to BQSR
else
  mark-section "not performing BQSR"
  recal_bam="$realigned_bam"
  recal_bai="$realigned_bai"
fi

mark-section "uploading results"
mkdir -p ~/out/bam/output/ ~/out/bai/output/ ~/out/outputmetrics/QC/
mv $recal_bam ~/out/bam/output/"$sorted_bam_prefix".refined.bam
mv $recal_bai ~/out/bai/output/"$sorted_bam_prefix".refined.bam.bai

if [[ "$skip_markduplicates" != "true" ]]
then
mark-section "uploading mark duplicates QC"
mv output.metrics ~/out/outputmetrics/QC/"$sorted_bam_prefix".output.metrics
fi


#
# Upload results
#
dx-upload-all-outputs --parallel
mark-success
