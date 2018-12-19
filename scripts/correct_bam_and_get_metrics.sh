#!/bin/bash

set -xueo pipefail

########################################################################################################
# PART 1: From raw basecalls to GATK-ready reads
########################################################################################################
## Clean BAM
#  - CleanSam using Picard tools
#  - FixMateInformation using Picard tools
#  - ValidateSamFile using Picard tools
#  - generate BAM index using samtools
## Metrics
#  - verify BAM id using verifyBamID
#  - get stats of BAM using Samtools flagstat
#  - wgs metrics using Picard tools CollectWgsMetrics
#  - insert size, alignment and GC bias metrics using Picard tools CollectMultipleMetrics
## File cleanup
#  - file cleanup
########################################################################################################

in_bam=$1
sanger_core_ref=$2
mem=$3

# untar ref
mkdir -p ref
tar xzf $sanger_core_ref -C ref --strip-components 1

path_ref="ref/genome.fa"
path_GATK="/opt/GenomeAnalysisTK.jar" # path to GATK jar
path_picard="/opt/picard.jar" # path to Picardtools jar
vcf_REF="/opt/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
tmp_bam_prefix=$(basename $in_bam)
tmp_bam_prefix=${tmp_bam_prefix%.*}  # remove the input bam extension

java_mem_tag="-Xmx${mem}g"

# expected outputs:

## bam and index
mateinfo_fixed_bam=${tmp_bam_prefix}.mateinfo_fixed.bam
# secondary file: ${tmp_bam_prefix}.mateinfo_fixed.bam.bai

## metrics_files
verifybamID_out=${tmp_bam_prefix}.verifybamID_out
# secondary file: ${tmp_bam_prefix}.verifybamID_out.depthSM
# secondary file: ${tmp_bam_prefix}.verifybamID_out.log
# secondary file: ${tmp_bam_prefix}.verifybamID_out.selfSM
flag_stats_out=${tmp_bam_prefix}.flag_stats.txt
wgs_metrics_out=${tmp_bam_prefix}.wgs_metrics.txt
multiple_metrics_out=${tmp_bam_prefix}.multiple_metrics
# secondary file: ${tmp_bam_prefix}.multiple_metrics.alignment_summary_metrics
# secondary file: ${tmp_bam_prefix}.multiple_metrics.gc_bias.detail_metrics
# secondary file: ${tmp_bam_prefix}.multiple_metrics.gc_bias.pdf
# secondary file: ${tmp_bam_prefix}.multiple_metrics.gc_bias.summary_metrics
# secondary file: ${tmp_bam_prefix}.multiple_metrics.insert_size_histogram.pdf
# secondary file: ${tmp_bam_prefix}.multiple_metrics.insert_size_metrics

## the log:
samplelog=${tmp_bam_prefix}.eze_gatk_part_1_correct.log
echo "$(date '+%d/%m/%y_%H:%M:%S'), Wake up to work" > "$samplelog"

####################################################
# Clean BAM
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Cleaning BAM---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard CleanSam \
TMP_DIR=/tmp \
I=$in_bam \
R=$path_ref \
O=${tmp_bam_prefix}.clean.bam >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished cleaning BAM---" >> "$samplelog"

####################################################
# FixMateInformation
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting FixMateInformation---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard FixMateInformation \
TMP_DIR=/tmp \
VALIDATION_STRINGENCY=LENIENT \
I=${tmp_bam_prefix}.clean.bam \
O=$mateinfo_fixed_bam \
TMP_DIR=/tmp >> "$samplelog"
rm ${tmp_bam_prefix}.clean.bam
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished Picard FixMateInformation" >> "$samplelog"

####################################################
# Generate BAM index
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Generating BAM index---" >> "$samplelog"
samtools index $mateinfo_fixed_bam
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished generating BAM index" >> "$samplelog"

####################################################
# verifyBamID
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting VerifyBamID---" >> "$samplelog"
verifyBamID \
--vcf $vcf_REF \
--bam $mateinfo_fixed_bam \
--out $verifybamID_out \
--ignoreRG \
--verbose >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished VerifyBamID" >> "$samplelog"

####################################################
# SAM FILE FLAG STATISTICS (samtools flagstat)
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Samtools flagstat---" >> "$samplelog"
samtools flagstat $mateinfo_fixed_bam >> $flag_stats_out
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished samtools flagstat" >> "$samplelog"

####################################################
# CollectWgsMetrics
####################################################

echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Picard CollectWgsMetrics---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard CollectWgsMetrics \
TMP_DIR=/tmp \
VALIDATION_STRINGENCY=LENIENT \
R=$path_ref \
I=$mateinfo_fixed_bam \
O=$wgs_metrics_out \
INCLUDE_BQ_HISTOGRAM=true >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished Picard CollectWgsMetrics" >> "$samplelog"

##################################################
# Collect multiple metrics
##################################################

# Collect metrics on insert size, GC bias, alignment summary
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Collecting multiple metrics---" >> "$samplelog"

java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard CollectMultipleMetrics \
TMP_DIR=/tmp \
R=$path_ref \
I=$mateinfo_fixed_bam \
O=$multiple_metrics_out \
PROGRAM=null \
PROGRAM=CollectAlignmentSummaryMetrics \
PROGRAM=CollectInsertSizeMetrics \
PROGRAM=CollectGcBiasMetrics \
METRIC_ACCUMULATION_LEVEL=null \
METRIC_ACCUMULATION_LEVEL=READ_GROUP \
METRIC_ACCUMULATION_LEVEL=SAMPLE \
VALIDATION_STRINGENCY=SILENT >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished collecting multiple metrics" >> "$samplelog"
