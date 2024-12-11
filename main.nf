nextflow.enable.dsl=2
import groovy.json.JsonOutput

include { GINKO_WF } from './nf-bioskryb-utils/subworkflows/cnv_ginko/main.nf' addParams( timestamp: params.timestamp )
include { printHeader; helpMessage } from './help' params ( params )

if ( params.help ) {
    helpMessage()
    exit 0
}

params.reference = params.genomes [ params.genome ] [ 'reference' ]
params.ginko_ref_dir = params.genomes [ params.genome ] [ "${params.bin_size}" ] [ "${params.ginko_readlen}" ] [ 'ginko_ref_dir' ]

workflow {
    
    if(params.bam_dir != ""){
        ch_bam = Channel.fromFilePairs(params.bam_dir + "/*" + params.sample_name + "*.{bam,bam.bai}")   
    }
    else if(params.input_csv != ""){
        ch_bam = Channel.fromPath( params.input_csv ).splitCsv( header:true )
                        .map { row -> [ row.biosampleName, row.bam, row.bam + ".bai" ] }
    }
    
    ch_reference = Channel.fromPath(params.reference)
    ch_binref = Channel.fromPath(params.ginko_ref_dir + "variable_" + params.bin_size + "_" + params.ginko_readlen + "_bwa")
    ch_gcref = Channel.fromPath(params.ginko_ref_dir + "GC_variable_" + params.bin_size + "_" + params.ginko_readlen + "_bwa")
    ch_boundsref_file = Channel.fromPath(params.ginko_ref_dir + "bounds_variable_" + params.bin_size + "_" + params.ginko_readlen + "_bwa")
    
    GINKO_WF(
            ch_bam,
            params.bin_size,
            ch_binref,
            ch_gcref,
            ch_boundsref_file,
            params.min_ploidy,
            params.max_ploidy,
            params.min_bin_width,
            params.is_haplotype,
            params.publish_dir,
            params.enable_publish,
            params.disable_publish
    )
    
    GINKO_WF.out.graph.flatten()
          .collectFile(name: "graph_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}")
          { it.getSimpleName().split('_')[3] + "\t" + "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginko/" +  it.getName() }
          
    GINKO_WF.out.CNVS.flatten()
          .collectFile(name: "cnvs_dots_files.txt", newLine: true, sort: { it[0] }, storeDir: "${params.tmp_dir}")
          { it.getSimpleName().split('_')[0] + "\t" + "${params.publish_dir}_${params.timestamp}/tertiary_analyses/cnv_ginko/"  + it.getName() }
    
    
}

workflow.onComplete {
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["graph"]       = [:]
    output["output"]["CNV1"]       = [:]
    output["output"]["CNV2"]       = [:]
    output["output"]["dots"]       = [:]

    
    graph_outfile = file("$params.tmp_dir/graph_files.txt")
    graph_outfile_lines = graph_outfile.readLines()
    for ( graph_line : graph_outfile_lines ) {
        def (sample_name, graph_path) = graph_line.split('\t')
        output["output"]["graph"][sample_name] = [:]
        output["output"]["graph"][sample_name]["jpeg"] = graph_path
    }
    
    tsvs_outfile = file("$params.tmp_dir/cnvs_dots_files.txt")
    tsvs_outfile_lines = tsvs_outfile.readLines()
    for ( tsvs_line : tsvs_outfile_lines ) {
        def (sample_name, tsvs_path) = tsvs_line.split('\t')
        if(tsvs_path.contains("CNV1")){
            output["output"]["CNV1"][sample_name] = [:]
            output["output"]["CNV1"][sample_name]["tsv"] = tsvs_path
        }
        else if(tsvs_path.contains("CNV2")){
            output["output"]["CNV2"][sample_name] = [:]
            output["output"]["CNV2"][sample_name]["tsv"] = tsvs_path
        }
        else if(tsvs_path.contains("dots")){
            output["output"]["dots"][sample_name] = [:]
            output["output"]["dots"][sample_name]["tsv"] = tsvs_path
        }
        
    }
    
    def output_json = JsonOutput.toJson(output)
    def output_json_pretty = JsonOutput.prettyPrint(output_json)
    File outputfile = new File("$params.tmp_dir/output.json")
    outputfile.write(output_json_pretty)
    // println(output_json_pretty)
}


workflow.onError {
    output                          = [:]
    output["pipeline_run_name"]     = workflow.runName
    output["pipeline_name"]         = workflow.manifest.name
    output["pipeline_version"]      = workflow.manifest.version
    output["pipeline_session_id"]   = workflow.sessionId
    output["output"]                = [:]
    output["output"]["fastq"]       = [:] 
    
    

    def subject = """\
        [nf-SV-pipeline] FAILED: ${workflow.runName}
        """

    def msg = """\
        Pipeline execution summary 
        --------------------------------
        Script name       : ${workflow.scriptName ?: '-'}
        Script ID         : ${workflow.scriptId ?: '-'}
        Workflow session  : ${workflow.sessionId}
        Workflow repo     : ${workflow.repository ?: '-' }
        Workflow revision : ${workflow.repository ? "$workflow.revision ($workflow.commitId)" : '-'}
        Workflow profile  : ${workflow.profile ?: '-'}
        Workflow cmdline  : ${workflow.commandLine ?: '-'}
        Nextflow version  : ${workflow.nextflow.version}, build ${workflow.nextflow.build} (${workflow.nextflow.timestamp})
        Error Report      : ${workflow.errorReport}
        """
        .stripIndent()
    
    log.info ( msg )
    
    if ( "${params.email_on_fail}" && workflow.exitStatus != 0 ) {
        sendMail(to: "${params.email_on_fail}", subject: subject, body: msg)
    }
}