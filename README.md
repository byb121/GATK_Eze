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

## Packaged CWL tools

### correct_bam_and_get_metrics

This is part 1 of the whole pipeline, which mark duplicates, correct BAM header and produces various quality metrics.

Details of inputs and outputs are in `cwl/correct_bam_and_get_metrics.json`.

***Note:*** The tool does **NOT** support multiple threading. 16GB of RAM should be ok for most of BAMs.

### bam_to_vcf

This is part 2 of the whole pipeline, which recalibrate a BAM and generate gzip GVCF file from it.

Details of inputs and outputs are in `cwl/bam_to_gvcf.json`.

***Note:*** Most of the processed in this tool support multiple threading, 8 cores should be enough for most of cases. Set a sensible amount of RAM accordingly.