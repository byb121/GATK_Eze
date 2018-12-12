#!/bin/bash

set -xue

########################################################################################################
# PART 2: From raw basecalls to GATK-ready reads
########################################################################################################

## Recalibration
#  - GATK BQSR (Base Quality Score Recalibration) using BaseRecalibrator
#  - GATK PrintReads
#  - GATK BQSR (Base Quality Score Recalibration) using BaseRecalibrator
#  - GATK AnalyzeCovariates
## Variant Calling
#  - GATK HaplotypeCaller
#  - bgzip + tabix .g.vcf
#  - GATK ValidateVariants
## File cleanup
#  - file cleanup

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
# just because GATK does not think genome.fa.dict is the dict file name!!
cp ref/genome.fa.dict ref/genome.dict

this_gatk="java -Xmx${mem}g -Djava.io.tmpdir=/tmp -Djava.library.path=/tmp -jar /opt/GenomeAnalysisTK.jar -R ref/genome.fa"

tmp_bam_prefix=$(basename $in_bam)
tmp_bam_prefix=${tmp_bam_prefix%.*}  # remove the input bam extension

samplelog=${tmp_bam_prefix}.eze_gatk_part_2_bam2gvcf.log
echo "$(date '+%d/%m/%y_%H:%M:%S'), Wake up to work" > "$samplelog"

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
-I $in_bam \
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
-I $in_bam \
-BQSR ${tmp_bam_prefix}.recal_table \
-o ${tmp_bam_prefix}.recal.bam) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK PrintReads finished---" >> "$samplelog"

####################################################
# GATK BQSR post
# Count covariates in recal.BAM	to compare before/after recalibration
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK BaseRecalibrator Post---" >> "$samplelog"

time ($this_gatk \
-nct $cpu \
-T BaseRecalibrator \
-I $in_bam \
-knownSites $gatk_ref_bundle_dbsnp \
-knownSites $gatk_ref_bundle_mills_1000g_gold_indel \
-knownSites $gatk_ref_bundle_1000g_indel \
-BQSR ${tmp_bam_prefix}.recal_table \
-o ${tmp_bam_prefix}.recal_table_post) >> "$samplelog" 

echo "$(date '+%d/%m/%y_%H:%M:%S'), Finished GATK BaseRecalibrator Post" >> "$samplelog"

####################################################
# GATK AnalyzeCovariates
# Create a pdf plot to see the results of recalibration
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GATK AnalyzeCovariates---" >> "$samplelog"

time ($this_gatk \
-T AnalyzeCovariates \
-before ${tmp_bam_prefix}.recal_table \
-after ${tmp_bam_prefix}.recal_table_post \
-plots ${tmp_bam_prefix}.recal_plots.pdf \
-csv ${tmp_bam_prefix}.recal_plots.csv \
-l DEBUG ) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'),---GATK AnalyzeCovariates finished---" >> "$samplelog"

####################################################
# Genotyping - HaplotypeCaller
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting HaplotypeCaller---" >> "$samplelog"
time ($this_gatk \
-nct $cpu \
-T HaplotypeCaller \
-I ${tmp_bam_prefix}.recal.bam \
--emitRefConfidence GVCF \
--dbsnp $gatk_ref_bundle_dbsnp \
-o ${tmp_bam_prefix}.g.vcf.gz \
-pairHMM VECTOR_LOGLESS_CACHING ) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'),---Finished HaplotypeCaller---" >> "$samplelog"

####################################################
# Check the validity of the g.vcf.gz file
####################################################
echo "$(date '+%d/%m/%y_%H:%M:%S'),---Starting GVCF validation after HaplotypeCaller---" >> "$samplelog"
time ($this_gatk \
-T ValidateVariants \
-V ${tmp_bam_prefix}.g.vcf.gz \
--validationTypeToExclude ALL ) >> "$samplelog"

echo "$(date '+%d/%m/%y_%H:%M:%S'),---GVCF validation after HaplotypeCaller COMPLETED---">> "$samplelog"