#!/bin/bash

set -xueo pipefail

########################################################################################################
# the script does the following
########################################################################################################
## Clean BAM
## Produce etrics
## Recalibration
## Variant Calling
########################################################################################################
########################################################################################################

in_bam=$1
# TODO check the index!!!
sanger_core_ref=$2
gatk_ref_bundle_gzipped_dbsnp=$3
gatk_ref_bundle_gzipped_mills_1000g_gold_indel=$4
gatk_ref_bundle_gzipped_1000g_indel=$5
mem=$6
cpu=$7

# untar ref
mkdir -p ref
tar xzf $sanger_core_ref -C ref --strip-components 1
path_ref=ref/genome.fa

this_gatk="java -Xmx${mem}g -Djava.io.tmpdir=/tmp -Djava.library.path=/tmp -jar /opt/GenomeAnalysisTK.jar -R $path_ref"
this_picard="java -Xmx${mem}g -Djava.io.tmpdir=/tmp -jar /opt/picard.jar"

tmp_bam_prefix=$(basename $in_bam)
tmp_bam_prefix=${tmp_bam_prefix%.*}  # remove the input bam extension

vcf_REF="/opt/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf"

# expected outputs:

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

samplelog=${tmp_bam_prefix}.bam2gvcf.log

## bam and index
recalibrated_bam=${tmp_bam_prefix}.recal.bam
# index ${tmp_bam_prefix}.recal.bai

## GVCF 
gvcf=${tmp_bam_prefix}.g.vcf.gz
# its index ${tmp_bam_prefix}.g.vcf.gz.tbi

echo "$(date '+%d/%m/%y_%H:%M:%S'), Wake up to work" > "$samplelog"

####################################################
# Clean BAM
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Cleaning BAM---" >> "$samplelog"
time ($this_picard CleanSam \
I=$in_bam \
R=$path_ref \
TMP_DIR=/tmp \
O=${tmp_bam_prefix}.clean.bam) >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished cleaning BAM---" >> "$samplelog"

####################################################
# FixMateInformation
####################################################
mateinfo_fixed_bam=${tmp_bam_prefix}.mateinfo_fixed.bam

echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting FixMateInformation---" >> "$samplelog"
time ($this_picard FixMateInformation \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/tmp \
I=${tmp_bam_prefix}.clean.bam \
O=$mateinfo_fixed_bam) >> "$samplelog"
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
time (verifyBamID \
--vcf $vcf_REF \
--bam $mateinfo_fixed_bam \
--out $verifybamID_out \
--ignoreRG \
--verbose) >> "$samplelog"
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
time ($this_picard CollectWgsMetrics \
VALIDATION_STRINGENCY=LENIENT \
TMP_DIR=/tmp \
R=$path_ref \
I=$mateinfo_fixed_bam \
O=$wgs_metrics_out \
INCLUDE_BQ_HISTOGRAM=true) >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished Picard CollectWgsMetrics" >> "$samplelog"

##################################################
# Collect multiple metrics
# Collect metrics on insert size, GC bias, alignment summary
##################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Collecting multiple metrics---" >> "$samplelog"

time ($this_picard CollectMultipleMetrics \
R=$path_ref \
TMP_DIR=/tmp \
I=$mateinfo_fixed_bam \
O=$multiple_metrics_out \
PROGRAM=null \
PROGRAM=CollectAlignmentSummaryMetrics \
PROGRAM=CollectInsertSizeMetrics \
PROGRAM=CollectGcBiasMetrics \
METRIC_ACCUMULATION_LEVEL=null \
METRIC_ACCUMULATION_LEVEL=READ_GROUP \
METRIC_ACCUMULATION_LEVEL=SAMPLE \
VALIDATION_STRINGENCY=SILENT) >> "$samplelog"
echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished collecting multiple metrics" >> "$samplelog"

########################################################################################################
# preparing ref files for GATK to use
########################################################################################################

# just because GATK does not think genome.fa.dict is the dict file name!!
cp ref/genome.fa.dict ref/genome.dict

echo "$(date '+%d/%m/%y_%H:%M:%S'), preparing GATK bundle files" > "$samplelog"
gatk_ref_bundle_dbsnp="ref/gatk_ref_bundle_dbsnp.vcf"
gatk_ref_bundle_mills_1000g_gold_indel="ref/gatk_ref_bundle_mills_1000g_gold_indel.vcf"
gatk_ref_bundle_1000g_indel="ref/gatk_ref_bundle_1000g_indel.vcf"
gunzip -c $gatk_ref_bundle_gzipped_dbsnp > $gatk_ref_bundle_dbsnp
gunzip -c $gatk_ref_bundle_gzipped_mills_1000g_gold_indel > $gatk_ref_bundle_mills_1000g_gold_indel
gunzip -c $gatk_ref_bundle_gzipped_1000g_indel > $gatk_ref_bundle_1000g_indel


########################################################################################################
# GATK BQSR (Base Quality Score Recalibration) 
# STEP 1: Creating covariates table masking known indel, SNP sites
########################################################################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK BaseRecalibrator Pre---" >> "$samplelog"

time ($this_gatk \
-nct $cpu \
-T BaseRecalibrator \
-I $mateinfo_fixed_bam \
-knownSites $gatk_ref_bundle_dbsnp \
-knownSites $gatk_ref_bundle_mills_1000g_gold_indel \
-knownSites $gatk_ref_bundle_1000g_indel \
-o ${tmp_bam_prefix}.recal_table) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished GATK BaseRecalibrator Pre" >> "$samplelog"

########################################################################################################
# GATK BQSR (Base Quality Score Recalibration) 
# STEP 2: Writing the new BAM with recalibrated Q scores
########################################################################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK PrintReads---" >> "$samplelog"

time ($this_gatk \
-nct $cpu \
-T PrintReads \
-I $mateinfo_fixed_bam \
-BQSR ${tmp_bam_prefix}.recal_table \
-o $recalibrated_bam) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK PrintReads finished---" >> "$samplelog"
rm $mateinfo_fixed_bam

####################################################
# Genotyping - HaplotypeCaller
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting HaplotypeCaller---" >> "$samplelog"
time ($this_gatk \
-nct $cpu \
-T HaplotypeCaller \
-I $recalibrated_bam \
--emitRefConfidence GVCF \
--dbsnp $gatk_ref_bundle_dbsnp \
-o $gvcf \
-pairHMM VECTOR_LOGLESS_CACHING ) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished HaplotypeCaller---" >> "$samplelog"

####################################################
# Check the validity of the g.vcf.gz file
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GVCF validation after HaplotypeCaller---" >> "$samplelog"
time ($this_gatk \
-T ValidateVariants \
-V $gvcf \
--validationTypeToExclude ALL ) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'),---GVCF validation after HaplotypeCaller COMPLETED---">> "$samplelog"