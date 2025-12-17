#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// --------------------------------------------------
// DEFAULT PARAMETERS
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

// --------------------------------------------------
// INFO / HELP
// --------------------------------------------------

log.info ""
log.info "-----------------------------------------------------------------"
log.info "BQSR-nf 1.1: BASE QUALITY SCORE RECALIBRATION"
log.info "-----------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-----------------------------------------------------------------"
log.info ""

if (params.help) {
    log.info ''
	log.info '-------------------------------------------------------------'
    log.info 'Usage:'
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

// --------------------------------------------------
// REFERENCE FILES
// --------------------------------------------------

ref          = file(params.ref)
ref_fai      = file("${params.ref}.fai")
ref_dict     = file(params.ref.replaceFirst(/fasta|fa/, '') + 'dict')

known_snps         = file(params.snp_vcf)
known_snps_index   = file("${params.snp_vcf}.tbi")
known_indels       = file(params.indel_vcf)
known_indels_index = file("${params.indel_vcf}.tbi")

ch_config_for_multiqc = Channel.value(file(params.multiqc_config))

// --------------------------------------------------
// INPUT BAM + BAI
// --------------------------------------------------

Channel bam_bai_files

if (file(params.input_folder).listFiles().any { it.name.endsWith('.bam') }) {

    bam_bai_files = Channel
        .fromFilePairs("${params.input_folder}/*.{bam,bai}", size: 2)
        .map { id, files -> tuple(files[0], files[1]) }
	println "BAM files found, proceed with realignment";
} else {
    error "ERROR: input folder contains no BAM files"
}

// --------------------------------------------------
// PROCESSES
// --------------------------------------------------

process BASE_QUALITY_SCORE_RECALIBRATION {

    tag { bam.baseName }
    cpus params.cpu
    memory "${params.mem}G"

    publishDir "${params.output_folder}/BAM", mode: 'copy', pattern: "*bam*"
    publishDir "${params.output_folder}/QC/BAM/BQSR", mode: 'copy',
        saveAs: { f ->
            f.contains('table') || f.contains('plots') ? f : null
        }

    input:
    tuple path(bam), path(bai)
    path known_snps
    path known_snps_index
    path known_indels
    path known_indels_index
    path ref
    path ref_fai
    path ref_dict

    output:
    path "*_recal.table", emit: recal_tables
    path "*plots.pdf",     emit: recal_plots
    tuple val(file_tag_new),
          path("${file_tag_new}.bam"),
          path("${file_tag_new}.bam.bai"),
          emit: final_bams

    script:
    """
    file_tag=${bam.baseName}
    file_tag_new=\${file_tag}_BQSRecalibrated

    gatk BaseRecalibrator \
        --java-options "-Xmx${params.mem}G" \
        -R ${ref} \
        -I ${bam} \
        --known-sites ${known_snps} \
        --known-sites ${known_indels} \
        -O \${file_tag}_recal.table

    gatk ApplyBQSR \
        --java-options "-Xmx${params.mem}G" \
        -R ${ref} \
        -I ${bam} \
        --bqsr-recal-file \${file_tag}_recal.table \
        -O \${file_tag_new}.bam

    gatk BaseRecalibrator \
        --java-options "-Xmx${params.mem}G" \
        -R ${ref} \
        -I \${file_tag_new}.bam \
        --known-sites ${known_snps} \
        --known-sites ${known_indels} \
        -O \${file_tag_new}_recal.table

    gatk AnalyzeCovariates \
        --java-options "-Xmx${params.mem}G" \
        -before \${file_tag}_recal.table \
        -after  \${file_tag_new}_recal.table \
        -plots  \${file_tag_new}_recalibration_plots.pdf

    mv \${file_tag_new}.bai \${file_tag_new}.bam.bai
    """
}

process MULTIQC_FINAL {

    cpus 2
    memory '1G'

    publishDir "${params.output_folder}/QC/BAM", mode: 'copy'

    input:
    path recal_tables
    path recal_plots
    path multiqc_config

    output:
    path "*report.html", emit: report
    path "multiqc_BQSR_report_data", emit: report_data

    script:
    """
    if [ "${multiqc_config.name}" = "NO_FILE" ]; then
        opt=""
    else
        opt="--config ${multiqc_config}"
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

    BASE_QUALITY_SCORE_RECALIBRATION(
        bam_bai_files,
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
        ch_config_for_multiqc
    )
}
