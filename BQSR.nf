#!/usr/bin/env nextflow

// Copyright (C) 2026 IARC/WHO
// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details <http://www.gnu.org/licenses/>.

nextflow.enable.dsl = 2

// --------------------------------------------------
// PARAMETERS
// --------------------------------------------------

params.help            = null
params.input_folder    = '.'
params.ref             = 'hg19.fasta'
params.cpu             = 2
params.mem             = 32
params.output_folder   = '.'
params.snp_vcf         = 'dbsnp.vcf'
params.indel_vcf       = 'Mills_1000G_indels.vcf'
params.multiqc_config  = 'NO_FILE'

//Header for the IARC tools - logo generated using the following page : http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}

// --------------------------------------------------
// FILE DEFINITION
// --------------------------------------------------

ref          = file(params.ref)
ref_fai      = file("${params.ref}.fai")
ref_dict     = file(params.ref.replaceFirst(/fasta|fa/, '') + 'dict')

known_snps         = params.snp_vcf ? file(params.snp_vcf) : null
known_snps_index   = file("${params.snp_vcf}.tbi")
known_indels       = params.indel_vcf ? file(params.indel_vcf) : null
known_indels_index = file("${params.indel_vcf}.tbi")

//ch_config_for_multiqc = Channel.value(file(params.multiqc_config))
multiqc = params.multiqc_config == 'NO_FILE'
    ? Channel.empty()
    : file(params.multiqc_config)

// --------------------------------------------------
// PROCESSES
// --------------------------------------------------

process BASE_QUALITY_SCORE_RECALIBRATION {

    tag { bam_tag }
    cpus params.cpu
    memory "${params.mem}G"

    input:
    tuple val(bam_tag), path(bam), path(bai)
    path known_snps
    path known_snps_index
    path known_indels
    path known_indels_index
    path ref
    path ref_fai
    path ref_dict

    output:
	path("${bam_tag}_BQSRecalibrated.bam"), emit: bam_out
	path("${bam_tag}_BQSRecalibrated.bai"), emit: bai_out
    path "*_recal.table",  emit: recal_tables
    path "*plots.pdf",     emit: recal_plots

    publishDir "${params.output_folder}/BQSR", mode: 'copy'

    script:
    """
	#!/bin/bash
    set -euo pipefail

	if [ ! -f "${bam}.bai" ]; then
		ln -s ${bai} ${bam}.bai
	fi

	//echo "[INFO] Running BQSR" - with AddOrReplaceReadGroups to avoid errors

	gatk AddOrReplaceReadGroups \
    	-I ${bam} \
    	-O ${bam_tag}_fixed.bam \
    	--RGID ${bam_tag} \
    	--RGLB lib1 \
    	--RGPL ILLUMINA \
    	--RGPU unit1 \
    	--RGSM ${bam_tag}

    gatk BaseRecalibrator \
        --java-options "-Xmx${params.mem}G" \
        -R ${ref} \
        -I ${bam}_fixed.bam \
        --known-sites ${known_snps} \
        --known-sites ${known_indels} \
        -O ${bam_tag}_recal.table

    gatk ApplyBQSR \
        --java-options "-Xmx${params.mem}G" \
        -R ${ref} \
        -I ${bam}_fixed.bam \
        --bqsr-recal-file ${bam_tag}_recal.table \
        -O ${bam_tag}_BQSRecalibrated.bam -- CREATE_INDEX true

    gatk BaseRecalibrator \
        --java-options "-Xmx${params.mem}G" \
        -R ${ref} \
        -I ${bam_tag}_BQSRecalibrated.bam \
        --known-sites ${known_snps} \
        --known-sites ${known_indels} \
        -O ${bam_tag}_recal.table

    gatk AnalyzeCovariates \
        --java-options "-Xmx${params.mem}G" \
        -before ${bam_tag}_recal.table \
        -after  ${bam_tag}_BQSRecalibrated_recal.table \
        -plots  ${bam_tag}_BQSRecalibrated_recalibration_plots.pdf

    mv ${bam_tag}_BQSRecalibrated.bai ${bam_tag}_BQSRecalibrated.bam.bai
    """
}

process MULTIQC_FINAL {
    cpus 2
    memory '1G'

    input:
    path recal_tables
    path recal_plots
    path multiqc_config //optional true

    output:
    path "*report.html", emit: report
    path "multiqc_BQSR_report_data", emit: report_data

    publishDir "${params.output_folder}/QC", mode: 'copy'

    script:
    """
    set -euo pipefail
	config_file='${multiqc_config}'
  	if [ "\$(basename "\$config_file")" = "NO_FILE" ]; then
        	opt=""
    	else
        	opt="--config \$config_file"
    	fi

    multiqc . \
        -n multiqc_BQSR_report.html \
        -m gatk \
        \$opt \
        --comment "GATK base quality score recalibration QC report"
    """
}

// --------------------------------------------------
// WORKFLOW
// --------------------------------------------------

workflow {

  		log.info IARC_Header()
// --------------------------------------------------
// INFO / HELP
// --------------------------------------------------

log.info ""
log.info "----------------------------------------------------------------------------------------------------------------"
log.info "BQSR-nf 1.1: BASE QUALITY SCORE RECALIBRATION"
log.info "----------------------------------------------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it under certain conditions; see LICENSE for details."
log.info "----------------------------------------------------------------------------------------------------------------"
log.info ""

if (params.help) {
    log.info ''
	log.info '-------------------------------------------------------------'
    log.info 'USAGE:'
    log.info 'nextflow run iarcbioinfo/BQSR-nf --input_folder input/ --ref hg19.fasta [--cpu 8] [--mem 32] [--output_folder output/]'
	log.info '-------------------------------------------------------------'
    log.info ''
	log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                 Folder containing BAM or fastq files to be aligned.'
    log.info '    --ref            FILE                   Reference fasta file (with index).'
    log.info 'Optional arguments:'
    log.info '    --cpu            INTEGER                Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem            INTEGER                Size of memory used by sambamba (in GB) (default: 32).'
    log.info '    --snp_vcf        STRING                 path to SNP VCF from GATK bundle (default : dbsnp.vcf)'
    log.info '    --indel_vcf      STRING                 path to indel VCF from GATK bundle (default : Mills_1000G_indels.vcf)'
    log.info '    --output_folder  STRING                Output folder (default: results_alignment).'
    log.info '    --multiqc_config STRING                 config yaml file for multiqc (default : none)'
    log.info ''
    exit 0
}

 else {
      /* Software information */
   log.info "input_folder = ${params.input_folder}"
   log.info "ref          = ${params.ref}"
   log.info "cpu          = ${params.cpu}"
   log.info "mem          = ${params.mem}"
   log.info "output_folder= ${params.output_folder}"
   log.info "snp_vcf          = ${params.snp_vcf}"
   log.info "indel_vcf          = ${params.indel_vcf}"
   log.info "help=${params.help}"
 }

        log.info "Running Base Quality Score Recalibration"
        bams = Channel.fromPath("${params.input_folder}/*.bam")
			.map { f -> tuple(f.baseName, f) }
			.ifEmpty { error "No BAM files found in ${params.input_folder}" }
        bais = Channel.fromPath("${params.input_folder}/*.bai")
            .map { f -> tuple(f.baseName, f) }
			.ifEmpty { error "No BAI files found in ${params.input_folder}" }
        
		bam_bai = bams.join(bais)
					.map {tag, bam, bai -> tuple(tag, bam, bai)} // emit tag, bam, bai
        bam_bai.view { "BAM_BAI → $it" }

    BASE_QUALITY_SCORE_RECALIBRATION(
        bam_bai,
        known_snps,
		known_snps_index,
        known_indels,
        known_indels_index,
        ref,
        ref_fai,
        ref_dict
    )

    MULTIQC_FINAL(
        BASE_QUALITY_SCORE_RECALIBRATION.out.recal_tables.collect(),
        BASE_QUALITY_SCORE_RECALIBRATION.out.recal_plots.collect(),
        multiqc
    )
}
