#!/usr/bin/env cwl-runner

$namespaces:
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/
  s: http://schema.org/

$schemas:
- https://schema.org/version/latest/schema.rdf
- http://dublincore.org/2012/06/14/dcterms.rdf
- http://xmlns.com/foaf/spec/20140114.rdf


class: CommandLineTool

id: "correct_bam_and_collect_metrics_for_gatk"

label: "Correct BMA and collect metrics for GATK pipeline"

cwlVersion: v1.0

doc: |
    ![build_status](https://quay.io/repository/wtsicgp/gatk_eze/status)
    See the [GATK_Eze](https://github.com/cancerit/GATK_Eze) website for more information.

dct:creator:
  "@id": "yaobo.xu@sanger.ac.uk"
  foaf:name: Yaobo Xu
  foaf:mbox: "yx2@sanger.ac.uk"

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/wtsicgp/gatk_eze:0.2.2"

hints:
  - class: ResourceRequirement
    coresMin: 1 # all processes are single threaded
    ramMin: 8000

inputs:

  in_bam: # Need the index file!!!
    type: File
    doc: "input BAM"
    inputBinding:
      position: 1
      shellQuote: true

  sanger_core_ref:
    type: File
    doc: "reference genome files in tar, as used in Sanger, GRCh37d5 can be downloaded here: http://ftp.sanger.ac.uk/pub/cancer/dockstore/human/core_ref_GRCh37d5.tar.gz"
    inputBinding:
      position: 2
      shellQuote: true

  gatk_ref_bundle_gzipped_dbsnp:
    type: File
    doc: "gzipped Miils and 1000G golden standard indels vcf file from GATK reference bundle, for Eze's pipeline, please use Mills_and_1000G_gold_standard.indels.b37.vcf.gz. It can be downloaded from ftp://ftp.broadinstitute.org/bundle/, username: gsapubftp-anonymous, no password"
    inputBinding:
      position: 3
      shellQuote: true

  gatk_ref_bundle_gzipped_mills_1000g_gold_indel:
    type: File
    doc: "gzipped 1000G phase1 indels vcf file from GATK reference bundle, for Eze's pipeline, please use 1000G_phase1.indels.b37.vcf.gz. It can be downloaded from ftp://ftp.broadinstitute.org/bundle/, username: gsapubftp-anonymous, no password"
    inputBinding:
      position: 4
      shellQuote: true

  gatk_ref_bundle_gzipped_1000g_indel:
    type: File
    doc: ""
    inputBinding:
      position: 5
      shellQuote: true

  mem:
    type: int
    doc: "number of GBs for max memory of Java processes to use"
    default: 12
    inputBinding:
      position: 6

  cpu:
    type: int
    doc: "number of GBs for max memory of Java processes to use"
    default: 1
    inputBinding:
      position: 7

outputs:

  g_vcf:
    type: File
    outputBinding:
      glob: $(inputs.in_bam.nameroot).g.vcf.gz
    secondaryFiles:
      - .tbi

  job_log:
    type: File
    outputBinding:
      glob: $(inputs.in_bam.nameroot).bam2gvcf.log

  flag_stats_out:
    type: File
    outputBinding:
      glob: $(inputs.in_bam.nameroot).flag_stats.txt

  wgs_metrics_out:
    type: File
    outputBinding:
      glob: $(inputs.in_bam.nameroot).wgs_metrics.txt

  verifybamID_out:
    type:
      type: array
      items: File
    outputBinding:
      glob: "$(inputs.in_bam.nameroot).verifybamID_out.*"
 
  multiple_metrics_out:
    type:
      type: array
      items: File
    outputBinding:
      glob: "$(inputs.in_bam.nameroot).multiple_metrics.*"

baseCommand: ["bam_to_gvcf.sh"]

s:codeRepository: https://github.com/cancerit/GATK_Eze
s:author:
  - class: s:Person
    s:email: mailto:cgphelp@sanger.ac.uk
    s:name: Yaobo Xu
