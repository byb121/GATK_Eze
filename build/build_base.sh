#!/bin/bash

set -uxe

# Basic preparation for apt package installation
apt-get -yq update
apt-get -yq install --no-install-recommends curl software-properties-common apt-transport-https ca-certificates

## install JAVA
add-apt-repository -y ppa:webupd8team/java
bash -c '/bin/echo debconf shared/accepted-oracle-license-v1-1 select true | /usr/bin/debconf-set-selections'
bash -c '/bin/echo debconf shared/accepted-oracle-license-v1-1 seen true | /usr/bin/debconf-set-selections'
apt-get -yq update
apt-get -yq install --no-install-recommends oracle-java8-installer oracle-java8-unlimited-jce-policy

# Install Samtools 1.7
apt-get -yq install --no-install-recommends samtools=1.7-1

# Install GATK 3.7
curl -sSL 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.7-0-gcfedb67' > GenomeAnalysisTK.tar.bz2
tar xjf GenomeAnalysisTK.tar.bz2
rm GenomeAnalysisTK.tar.bz2
rm -r resources

# Install Picard 2.0.1
DEBIAN_FRONTEND=noninteractive apt-get -yq install --no-install-recommends ant git r-base r-base-dev
curl -sSL https://github.com/broadinstitute/picard/archive/2.0.1.tar.gz > picard.tar.gz
tar xzf picard.tar.gz
rm picard.tar.gz
cd picard-2.0.1
curl -sSL https://github.com/samtools/htsjdk/archive/2.0.1.tar.gz > htsjdk.tar.gz
tar xzf htsjdk.tar.gz
rm htsjdk.tar.gz
mv htsjdk-2.0.1 htsjdk
ant clean all
mv dist/picard.jar ../
ant clean
cd ..
rm -r picard-2.0.1

# Install VerifyBamID v1.1.2
curl -sSL https://github.com/statgen/verifyBamID/releases/download/v1.1.2/verifyBamID.1.1.2 > verifyBamID
chmod +x verifyBamID
# download ref vcf
curl -sSL http://csg.sph.umich.edu/kang/verifyBamID/download/Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf.gz | gunzip -c > Omni25_genotypes_1525_samples_v2.b37.PASS.ALL.sites.vcf

apt-get clean autoclean
apt-get autoremove -y

Rscript -e 'install.packages(c("ggplot2", "gplots", "reshape", "gsalib"))'  # for GATK to use