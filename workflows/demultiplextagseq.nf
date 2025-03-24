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

// my process to combibe the umi extracted fastq files
process combine_fastq {
    input :
    tuple val(meta), path(read)

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: reads

    script:
    """
    cat *extract_1.fastq.gz > ${meta.id}_1.fastq.gz
    cat *extract_2.fastq.gz > ${meta.id}_2.fastq.gz
    """
}


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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    //
    // Split the fastq files if parameter is set
    //
    if (params.split_fastq) {
        fq_split = ch_samplesheet.map { meta, reads -> 
                return tuple(meta, reads[0], reads[1]) 
        }
            .splitFastq(pe: true, by: 5000, file: true)
                .map { meta, read1, read2->
                    return tuple([id: meta.id + "_" + read1.baseName, single_end: meta.single_end],
                        [read1, read2])
                }
    } else {
      fq_split = ch_samplesheet
    }  
    //
    // MODULE: Extract UMI
    //
    UMITOOLS_EXTRACT (
        fq_split
    )
    ch_multiqc_files = ch_multiqc_files.mix(UMITOOLS_EXTRACT.out.log.collect{it[1]})
    ch_versions = ch_versions.mix(UMITOOLS_EXTRACT.out.versions.first())
    //
    // Combine the umi extracted fastq files
    //
    if (params.split_fastq) {
        UMITOOLS_EXTRACT.out.reads.map{ meta, file -> 
            return tuple( [id: meta.id.take(meta.id.lastIndexOf('.')), single_en: meta.single_end], file)}
            .set{ reads }

        reads = reads.groupTuple().map { meta, reads -> return tuple(meta, reads.flatten())}
        combine_fastq(
        reads
        )
        combine_fastq = combine_fastq.out.reads
    
    } else {
         combine_fastq = UMITOOLS_EXTRACT.out.reads
    }
    //
    // MODULE: demultiplexing
    //
    FQTK (
        combine_fastq
    ) 
    ch_multiqc_files = ch_multiqc_files.mix(FQTK.out.metrics.map { meta, metrics -> return metrics} )
    ch_versions = ch_versions.mix(FQTK.out.versions.first())  
    // 
    // MODULE: Run FastQC on demultiplexed reads
    // 
    FASTQC_POST (
        FQTK.out.sample_fastq.transpose().map { meta, fastq 
        -> def baseName = fastq.name.split('\\.')[0] 
        def namedTuple = [id: baseName, single_end: true ] 
        return tuple (namedTuple,fastq)}
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST.out.zip.collect{it[1]})
    //
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
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
        
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
