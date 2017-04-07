#! /usr/bin/env nextflow
// usage : ./BQSR.nf --input_folder input/ --cpu 8 --mem 32 --fasta_ref hg19.fasta
/*
vim: syntax=groovy
-*- mode: groovy;-*- */

// requirement:
// - samtools
// - samblaster
// - sambamba

//default values
params.help         = null
params.input_folder = '.'
params.fasta_ref    = 'hg19.fasta'
params.cpu          = 8
params.mem          = 32
params.mem_sambamba = 1
params.out_folder   = "."
params.intervals    = ""
params.GATK_bundle  = "bundle"
params.GATK_folder  = "."

if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW BASE QUALITY SCORE RECALIBRATION SCRIPT'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run BQSR.nf --input_folder input/ --fasta_ref hg19.fasta [--cpu 8] [--mem 32] [--out_folder output/]'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --input_folder   FOLDER                  Folder containing BAM or fastq files to be aligned.'
    log.info '    --fasta_ref          FILE                    Reference fasta file (with index).'
    log.info 'Optional arguments:'
    log.info '    --cpu          INTEGER                 Number of cpu used by bwa mem and sambamba (default: 8).'
    log.info '    --mem          INTEGER                 Size of memory used by sambamba (in GB) (default: 32).'
    log.info '    --out_folder     STRING                Output folder (default: results_alignment).'
    log.info ''
    exit 1
}

//read files
fasta_ref = file(params.fasta_ref)
fasta_ref_fai = file( params.fasta_ref+'.fai' )
fasta_ref_sa = file( params.fasta_ref+'.sa' )
fasta_ref_bwt = file( params.fasta_ref+'.bwt' )
fasta_ref_ann = file( params.fasta_ref+'.ann' )
fasta_ref_amb = file( params.fasta_ref+'.amb' )
fasta_ref_pac = file( params.fasta_ref+'.pac' )
fasta_ref_alt = file( params.fasta_ref+'.alt' )

if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
       println "BAM files found, proceed with realignment";
       bam_files = Channel.fromPath( params.input_folder+'/*.bam');
       bai_files = Channel.fromPath( params.input_folder+'/*.bai')
}else{
       println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
}

// base quality score recalibration
process base_quality_score_recalibration {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }
        
    input:
    file("${file_tag}.bam") from bam_files
    file("${file_tag}.bam.bai") from bai_files
    output:
    file("${file_tag_new}_recal.table") into recal_table_files
    file("${file_tag_new}_post_recal.table") into recal_table_post_files
    file("${file_tag_new}_recalibration_plots.pdf") into recal_plots_files
    set val(file_tag_new), file("${file_tag_new}.bam") into recal_bam_files
    file("${file_tag_new}.bam.bai") into recal_bai_files
    publishDir params.out_folder, mode: 'move'

    shell:
    file_tag=infile.baseName
    file_tag_new=file_tag+'_BQSR'
    '''
    indelsvcf=(`ls !{params.GATK_bundle}/*indels*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    dbsnpvcfs=(`ls !{params.GATK_bundle}/*dbsnp*.vcf* | grep -v ".tbi" | grep -v ".idx"`)
    dbsnpvcf=${dbsnpvcfs[@]:(-1)}
    knownSitescom=''
    for ll in $indelsvcf; do knownSitescom=$knownSitescom' -knownSites '$ll; done
    knownSitescom=$knownSitescom' -knownSites '$dbsnpvcf
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam $knownSitescom -L !{params.intervals} -o !{file_tag_new}_recal.table
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T BaseRecalibrator -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam $knownSitescom -BQSR !{file_tag_new}_recal.table -L !{params.intervals} -o !{file_tag_new}_post_recal.table		
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T AnalyzeCovariates -R !{params.fasta_ref} -before !{file_tag_new}_recal.table -after !{file_tag_new}_post_recal.table -plots !{file_tag_new}_recalibration_plots.pdf	
    java -jar !{params.GATK_folder}/GenomeAnalysisTK.jar -T PrintReads -nct !{params.cpu} -R !{params.fasta_ref} -I !{file_tag}.bam -BQSR !{file_tag_new}_recal.table -L !{params.intervals} -o !{file_tag_new}.bam
    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    '''
}