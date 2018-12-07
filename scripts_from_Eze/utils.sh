#!/bin/sh

####################################################
# Universal paths
####################################################

root_folder="/path/to/project/folder" # eg "ICGC"
bundle2_8="/path/to/bundle/" # not needed for part 1, can be removed (it can be downloaded from ftp://ftp.broadinstitute.org/bundle//, username: gsapubftp-anonymous, password empty)
path_GATK="/apps/gatk/3.7-0/GenomeAnalysisTK.jar" # path to GATK jar
path_picard="/apps/picard-tools/2.0.1/picard.jar" # path to Picardtools jar
path_ref="$root_folder/human_37_ICGC_Sanger_genome.fa" # path to reference genome
vcf_REF="$root_folder/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf" # can be downloaded from http://www.google.com/url?q=http%3A%2F%2Fcsg.sph.umich.edu%2Fkang%2FverifyBamID%2Fdownload%2FOmni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf.gz&sa=D&sntz=1&usg=AFQjCNH_AcArf60EJcq6nl23dL8TOtaaGw
tmp="$root_folder/tmp"
path_rds_vcfs="/path/to/output/vcfs" # eg "$root_folder/output_vcfs", not needed for part 1, can be removed
logfile=$root_folder/logs/log.log # not needed for parts 1 and 2, can be removed
original_path="path/to/input/bams" # eg "$root_folder/input_bams"


