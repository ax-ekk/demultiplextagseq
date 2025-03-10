/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UMITOOLS_EXTRACT       } from '../modules/nf-core/umitools/extract/main'
include { FQTK                   } from '../modules/nf-core/fqtk/main'
include { FASTQC; FASTQC as FASTQC_POST } from '../modules/nf-core/fastqc/main' 
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_demultiplextagseq_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEMULTIPLEXTAGSEQ {

    take:
    ch_samplesheet // channel: samplesheet read in from --input,
    
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ meta, zip -> return tuple(meta.id,zip[0],zip[1]) })    
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //
    // MODULE: Extract UMI
    //
    UMITOOLS_EXTRACT (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.join(UMITOOLS_EXTRACT.out.log.map { meta, log -> return tuple(meta.id, log)})
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())
    //
    // MODULE: demultiplexing
    //
    FQTK (
        UMITOOLS_EXTRACT.out.reads
    ) 
    ch_multiqc_files = ch_multiqc_files.join(FQTK.out.metrics.map { meta, metrics -> return tuple(meta.id, metrics)} )
    ch_versions = ch_versions.mix(FQTK.out.versions.first())  
    // 
    // MODULE: Run FastQC on demultiplexed reads
    // 
    FASTQC_POST (
        FQTK.out.sample_fastq.transpose().map { meta, fastq 
        -> def baseName = fastq.name.split('\\.')[0] 
        def namedTuple = [id: meta.id+"_"+baseName, single_end: true ] 
        return tuple (namedTuple,fastq)}
    )
    // name using multiplex id (sample) to match with rest of multiqc files (one multiqc report per multiplex)
    post_fastqcs =  FASTQC_POST.out.zip.map { meta, zip 
    -> def sample = meta.id.split('_')[0]+"_PAIRED_END"
    return tuple(sample, zip)}

    ch_multiqc_files = post_fastqcs.groupTuple().map{ id  -> return id.flatten()}.combine(ch_multiqc_files,by:0).map{id -> return tuple(id[0],id[1..-1])}
    
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_extra_files      = 
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_extra_files          = ch_extra_files.mix(ch_collated_versions)
    ch_extra_files          = ch_extra_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        ))

    MULTIQC (
        ch_multiqc_files,
        ch_extra_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
       []
        
    )

    emit: multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
