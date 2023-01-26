process R_output_bed_files {
    // debug true
    // container 'emosyne/r_docker:1.92'
    container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "$EP_list"
    cache "lenient"
    

    input:
    tuple  val(EP_list), path(EP_file_hg19)

    output:
    tuple  val(EP_list), path("*EP_file_hg19.tsv"), emit: formattedfiles

    script:
    """

    R_output_bed_files.r ${EP_list} ${EP_file_hg19}
    

   
    """
}
