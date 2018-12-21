# GATK_Eze

Dockerised version of Ezequiel Anokian's GATK pipeline for germline variants analysis. For pipeline's details please contact: Ezequiel.Anokian@icr.ac.uk

## Packaged tool versions

* **Java:** `1.8.0_191`
* **GATK:** `v3.7-0-gcfedb67`
* **Picardtools:** `2.0.1` (built with `htsjdk-2.0.1`)
* **VerifyBamID:** `1.1.2`
* **samtools:** `1.7.1`
* **R:** `3.4.4` with packages below:
  * ggplot2: `3.1.0`
  * gplots: `3.0.1`
  * gsalib: `2.1`
  * reshape `0.8.8`
  And their dependencies.

## Packaged CWL tool

### bam_to_vcf

It produces various quality metrics from the input BAM, clean the BAM file and fix mate info of it, then recalibrates the BAM and generates a gzip GVCF file from it.

Details of inputs and outputs are in `cwl/bam_to_gvcf.json`.

***Note:*** Only BAM recalibration and HaplotypCaller support multiple threading, 8 cores should be enough for most of cases. Set a sensible amount of RAM accordingly.