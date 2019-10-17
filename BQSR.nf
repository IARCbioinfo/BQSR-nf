#! /usr/bin/env nextflow
// requirement:
// - gatk4
// - samblaster
// - sambamba

//default values
params.help         = null
params.input_folder = '.'
params.ref          = 'hg19.fasta'
params.cpu          = 2
params.mem          = 32
params.output_folder   = "."
params.snp_vcf      = "dbsnp.vcf"
params.indel_vcf    = "Mills_1000G_indels.vcf"
params.multiqc_config = 'NO_FILE'


if (params.help) {
    log.info ''
    log.info '-------------------------------------------------------------'
    log.info 'NEXTFLOW BASE QUALITY SCORE RECALIBRATION SCRIPT'
    log.info '-------------------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run BQSR.nf --input_folder input/ --ref hg19.fasta [--cpu 8] [--mem 32] [--output_folder output/]'
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
    log.info ''
    exit 1
}

//read files
ref = file(params.ref)
ref_fai = file( params.ref+'.fai' )
ref_dict= file( params.ref.replaceFirst(/fasta/, "").replaceFirst(/fa/, "") +'dict')
//ref_sa  = file( params.ref+'.sa' )
//ref_bwt = file( params.ref+'.bwt' )
//ref_ann = file( params.ref+'.ann' )
//ref_amb = file( params.ref+'.amb' )
//ref_pac = file( params.ref+'.pac' )
//ref_alt = file( params.ref+'.alt' )

//get know site VCFs from GATK bundle
known_snps         = file( params.snp_vcf )
known_snps_index   = file( params.snp_vcf+'.tbi' )
known_indels       = file( params.indel_vcf )
known_indels_index = file( params.indel_vcf+'.tbi' )

//multiqc config file
ch_config_for_multiqc = file(params.multiqc_config)

if (file(params.input_folder).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){
       println "BAM files found, proceed with realignment";
       //bam_files = Channel.fromPath( params.input_folder+'/*.bam');
       //bai_files = Channel.fromPath( params.input_folder+'/*.bai')
	bam_bai_files = Channel.fromFilePairs("${params.input_folder}/*{.bam,.bai}")
			   .map { row -> tuple(row[1][0], row[1][1]) }

}else{
       println "ERROR: input folder contains no fastq nor BAM files"; System.exit(0)
}

// base quality score recalibration
process base_quality_score_recalibration {
    cpus params.cpu
    memory params.mem+'G'
    tag { file_tag }
    
    publishDir "$params.output_folder/BAM/", mode: 'copy', pattern: "*bam*"
    publishDir "$params.output_folder/QC/BAM/BQSR/", mode: 'copy',
	saveAs: {filename -> 
		if (filename.indexOf("table") > 0) "$filename"
		else if (filename.indexOf("plots") > 0) "$filename"
		else null
	}

    input:
    set file(bam), file(bai) from bam_bai_files
    file known_snps
    file known_snps_index
    file known_indels
    file known_indels_index
    file ref
    file ref_fai
    file ref_dict

    output:
    file("*_recal.table") into recal_table_files
    file("*plots.pdf") into recal_plots_files
    set val(file_tag_new), file("${file_tag_new}.bam"), file("${file_tag_new}.bam.bai") into final_bam_bai_files

    shell:
    file_tag=bam.baseName
    file_tag_new=file_tag+'_BQSRecalibrated'
    '''
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem}G" -R !{ref} -I !{file_tag}.bam --known-sites !{known_snps} --known-sites !{known_indels} -O !{file_tag}_recal.table
    gatk ApplyBQSR --java-options "-Xmx!{params.mem}G" -R !{ref} -I !{file_tag}.bam --bqsr-recal-file !{file_tag}_recal.table -O !{file_tag_new}.bam
    gatk BaseRecalibrator --java-options "-Xmx!{params.mem}G" -R !{ref} -I !{file_tag_new}.bam --known-sites !{known_snps} --known-sites !{known_indels} -O !{file_tag_new}_recal.table		
    gatk AnalyzeCovariates --java-options "-Xmx!{params.mem}G" -before !{file_tag}_recal.table -after !{file_tag_new}_recal.table -plots !{file_tag_new}_recalibration_plots.pdf	
    mv !{file_tag_new}.bai !{file_tag_new}.bam.bai
    '''
}

process multiqc_final {
    cpus 2
    memory '1G'

    publishDir "${params.output_folder}/QC/BAM/", mode: 'copy'

    input:
    file BQSR_results from recal_table_files.collect()
    file BQSR_results_plots from recal_plots_files.collect()
    file multiqc_config from ch_config_for_multiqc    

    output:
    file("*report.html") into final_output
    file("multiqc_BQSR_report_data/") into final_output_data

    shell:
    if( multiqc_config.name=='NO_FILE' ){
        opt = ""
    }else{
        opt = "--config ${multiqc_config}"
    }
    '''
    multiqc . -n multiqc_BQSR_report.html -m gatk !{opt} --comment "GATK base quality score recalibration QC report"
    '''
}

// Display completion message
workflow.onComplete {
  log.info "N E X T F L O W  ~  version ${workflow.nextflow.version} ${workflow.nextflow.build}"
  //log.info "iarcbioinfo/BQSR-nf ~ " + this.grabRevision() + (workflow.commitId ? " [${workflow.commitId}]" : "")
  log.info "Completed at: " + workflow.complete
  log.info "Duration    : " + workflow.duration
  log.info "Success     : " + workflow.success
  log.info "Exit status : " + workflow.exitStatus
  log.info "Error report: " + (workflow.errorReport ?: '-')
}



