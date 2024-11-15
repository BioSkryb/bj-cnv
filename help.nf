nextflow.enable.dsl=2

def printHeader() {
  
  log.info """\
  BJ-CNV   P I P E L I N E
  ===================================
  input     : ${ params.input_csv }
  publish_dir         : ${ params.publish_dir }
  timestamp           : ${ params.timestamp }
  genome              : ${ params.genome }
  \n
  """

}

def helpMessage() {

  yellow = "\033[0;33m"
  blue = "\033[0;34m"
  white = "\033[0m"
  red = "\033[0;31m"

  log.info """\
${blue}
    bj-cnv pipeline

    Usage:
        nextflow run main.nf [options]

    Script Options: see nextflow.config

${red}
        [required]
        --input_csv         FILE    Path to input csv file

        --genome            STR     Reference genome to use. Available options - GRCh38, GRCm39
                                    DEFAULT: ${params.genome}

${yellow}
        [optional]
        
        --genomes_base      STR     Path to the genomes
                                    DEFAULT: ${params.genomes_base}

        --publish_dir       DIR     Path to run output directory
                                    DEFAULT: ${params.publish_dir}
                                    
        --timestamp         STR     User can specify timestamp otherwise uses runtime generated timestamp
                                    
        --help              BOOL    Display help message
                                    
${white}
    """.stripIndent()
}

workflow{
  printHeader()
  helpMessage()
}
