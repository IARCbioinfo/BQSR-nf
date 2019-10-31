# BQSR-nf
Nextflow script for base quality score recalibration of bam files using GATK

## Nextflow pipeline for base quality score recalibration with GATK processing
[![CircleCI](https://circleci.com/gh/IARCbioinfo/BQSR-nf/tree/master.svg?style=svg)](https://circleci.com/gh/IARCbioinfo/BQSR-nf/tree/master)
[![Docker Hub](https://img.shields.io/badge/docker-ready-blue.svg)](https://hub.docker.com/r/iarcbioinfo/bqsr-nf/)

## Decription

Nextflow pipeline for base quality score recalibration and quality control

## Dependencies

1. Nextflow: for common installation procedures see the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository.

2. [*multiQC*](http://multiqc.info/docs/)
3. [*GATK4*](https://software.broadinstitute.org/gatk/guide/quickstart) must be in the PATH variable
4. [GATK bundle](https://software.broadinstitute.org/gatk/download/bundle) VCF files with lists of indels and SNVs (recommended: 1000 genomes indels, dbsnp VCF)

You can provide a config file to customize the multiqc report (see https://multiqc.info/docs/#configuring-multiqc).

## Input 
 | Type      | Description     |
  |-----------|---------------|
  |--input_folder    | a folder with bam files |
 

## Parameters

* #### Mandatory
| Name | Example value | Description |
|-----------|--------------:|-------------| 
|--ref |    ref.fa | reference genome fasta file for GATK |

* #### Optional

| Name | Default value | Description |
|-----------|--------------|-------------| 
|--cpu          | 2 | number of CPUs |
|--mem         | 32 | memory for mapping|
|--output_folder   | . | output folder for aligned BAMs|
|--snp_vcf |  dbsnp.vcf | VCF file with known variants for GATK BQSR |
|--indel_vcf |  Mills_100G_indels.vcf | VCF file with known indels for GATK BQSR |
|--multiqc_config   |  null | config yaml file for multiqc | 

* #### Flags

| Name  | Description |
|-----------|-------------| 
|--help | print usage and optional parameters |

## Usage
To run the pipeline on a series of bam files in folder *bam*, a reference genome with indexes at *ref.fa*, and known snps and indels from the gatk bundle, one can type:
```bash
nextflow run iarcbioinfo/BQSR-nf --input_folder bam --ref ref.fa --snp_vcf GATK_bundle/dbsnp_146.hg38.vcf.gz --indel_vcf GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
``` 

## Output 
  | Type      | Description     |
  |-----------|---------------|
  | BAM/file.bam    | BAM files of alignments or realignments |
  | BAM/file.bam.bai    | BAI files of alignments or realignments |
  | QC/multiqc_BQSR_report.html      |     multiqc report  | 
  | QC/multiqc_BQSR_report_data      |  folder with data used to compute multiqc report |
  | QC/BAM/BQSR/file_recal.table | table of scores before recalibration   |
  | QC/BAM/BQSR/file_post_recal.table   | table of scores after recalibration |
  | QC/BAM/BQSR/file_recalibration_plots.pdf   |  before/after recalibration plots   |
          
The output_folder directory contains two subfolders: BAM and QC

## Directed Acyclic Graph

[![DAG BQSR](dag_BQSR.png)](http://htmlpreview.github.io/?https://github.com/IARCbioinfo/BQSR-nf/blob/dev/dag_BQSR.html)

## Contributions

  | Name      | Email | Description     |
  |-----------|---------------|-----------------| 
  | Nicolas Alcala*    | AlcalaN@fellows.iarc.fr    | Developer to contact for support |
