#!/bin/bash

set -xue

########################################################################################################
# PART 1: From raw basecalls to GATK-ready reads
########################################################################################################
## Clean BAM
#  - mark duplicates using Picard tools MarkDuplicates (removed)
#  - CleanSam using Picard tools
#  - FixMateInformation using Picard tools
#  - reformat BAM header (HiSeq -> illumina)
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
mkdir ref && tar xzf $sanger_core_ref -C ref --strip-components 1
path_ref="ref/genome.fa"

path_GATK="/opt/GenomeAnalysisTK.jar" # path to GATK jar
path_picard="/opt/picard.jar" # path to Picardtools jar
vcf_REF="/opt/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"
tmp_bam_prefix=$(basename $in_bam)

java_mem_tag="-Xmx${mem}g"

# expected outputs:
## bam and index
reheadered_bam=${tmp_bam_prefix}.reheader.bam
# ${tmp_bam_prefix}.reheader.bam.bai
## metrics_files
verifybamID_out=${tmp_bam_prefix}.verifybamID_out.txt
flag_stats_out=${tmp_bam_prefix}.flag_stats.txt
wgs_metrics_out=${tmp_bam_prefix}.wgs_metrics.txt
multiple_metrics_out=${tmp_bam_prefix}.multiple_metrics.txt
## the log:
samplelog=${tmp_bam_prefix}.eze_gatk_part_1_correct.log

###################################################
# Mark Duplicates
###################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Picard MarkDuplicates---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard MarkDuplicates \
I=$in_bam \
O=${tmp_bam_prefix}.md.bam \
M=${tmp_bam_prefix}.duplicate_metrics.txt \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
REMOVE_DUPLICATES=true >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'),Finished Picard MarkDuplicates" >> "$samplelog"

####################################################
# Clean BAM
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Cleaning BAM---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard CleanSam \
I=${tmp_bam_prefix}.md.bam \
R=$path_ref \
O=${tmp_bam_prefix}.clean.bam >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished cleaning BAM---" >> "$samplelog"

####################################################
# FixMateInformation
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting FixMateInformation---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard FixMateInformation \
I=${tmp_bam_prefix}.clean.bam \
O=${tmp_bam_prefix}.fixed.bam \
TMP_DIR=/tmp >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished Picard FixMateInformation" >> "$samplelog"

####################################################
# Reformat BAM header (PL: HiSeq -> illumina)
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Reformat header---" >> "$samplelog"

# NOTE:
# if sample belongs to "Pilot - Multicentric" or "Pilot - Asian - Shanghai", it doesn't have PL info in the header at all, we'll have to add it instead of replacing value, using the following command:
# samtools view -H $root_folder/input_bams/$sample.reheader.bam | sed -e 's/LB:/PL:illumina\tLB:/g' | samtools reheader - $root_folder/input_bams/$sample.reheader.bam > $root_folder/input_bams/$sample.reheader.with.PL.bam

# to add PL:illumina to BAM header
# samtools view -H $root_folder/input_bams/$sample.fixed.bam | sed -e 's/PU:/PL:illumina\tPU:/g' | samtools reheader - $root_folder/input_bams/$sample.fixed.bam > $root_folder/input_bams/$sample.reheader.bam \
# && echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished reformating header" >> "$samplelog" \
# && touch $root_folder/logs/part_1_Reformat_header_finished_$sample.txt

samtools view -H ${tmp_bam_prefix}.fixed.bam | sed -e 's/PL:HiSeq/PL:illumina/g' | samtools reheader - ${tmp_bam_prefix}.fixed.bam > $reheadered_bam
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished reformating header" >> "$samplelog"

####################################################
# Generate BAM index
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Generating BAM index---" >> "$samplelog"
samtools index $reheadered_bam
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished generating BAM index" >> "$samplelog"

####################################################
# Validate BAM
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Validating BAM---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard ValidateSamFile \
I=$reheadered_bam \
MODE=SUMMARY >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished BAM Validation---" >> "$samplelog"

####################################################
# verifyBamID
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting VerifyBamID---" >> "$samplelog"
verifyBamID \
--vcf $vcf_REF \
--bam $reheadered_bam \
--out $verifybamID_out \
--ignoreRG \
--verbose >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished VerifyBamID" >> "$samplelog"

####################################################
# SAM FILE FLAG STATISTICS (samtools flagstat)
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Samtools flagstat---" >> "$samplelog"
samtools flagstat $reheadered_bam >> $flag_stats_out
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished samtools flagstat" >> "$samplelog"

####################################################
# CollectWgsMetrics
####################################################

echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting Picard CollectWgsMetrics---" >> "$samplelog"
java $java_mem_tag -Djava.io.tmpdir=/tmp \
-jar $path_picard CollectWgsMetrics \
R=$path_ref \
I=$reheadered_bam \
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
R=$path_ref \
I=$reheadered_bam \
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