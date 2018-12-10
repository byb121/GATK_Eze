# GATK_Eze

Dockerised version of [Ezequiel Anokian](Ezequiel.Anokian@icr.ac.uk)'s GATK pipeline for germline variants analysis.

## Packaged tool versions

* **Java:** `1.8.0_191`
* **GATK:** `v3.7-0-gcfedb67`
* **Picardtools:** `2.0.1` (built with `htsjdk-2.0.1`)
* **VerifyBamID:** `1.1.2`
* **samtools:** `1.7.1`

## Packaged CWL tools

### correct_bam_and_get_metrics

This is part 1 of the whole pipeline, which mark duplicates, correct BAM header and produces various quality metrics.

Details of inputs and outputs are in `cwl/correct_bam_and_get_metrics.json`.